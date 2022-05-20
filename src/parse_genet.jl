function parse_ref(ref_file::String, chroms::UnitRange)
    df = CSV.File(ref_file; types=Dict(:A1=>Char,:A2=>Char)) |> DataFrame
    df.A1 = tochar.(df.A1)
    df.A2 = tochar.(df.A2)
    @assert df.CHR isa Vector{Int}
    @assert df.BP isa Vector{Int}
    @assert df.MAF isa Vector{T} where T<:Real
    filter!(row->row.CHR in chroms, df)
    return df
end
parse_ref(ref_file::String, chrom::Integer) =
    parse_ref(ref_file, chrom:chrom)
function tochar(x)::Union{Char,Missing}
    if x isa String
        if length(x) == 1
            return uppercase(first.(x))
        else
            return missing
        end
    elseif x isa Char
        uppercase(x)
    else
        x
    end
end

function parse_bim(bim_file::String, chroms::UnitRange)
    header = [:CHR, :SNP, :POS, :BP, :A1, :A2]
    df = CSV.File(bim_file*".bim"; header=header, types=Dict(:A1=>Char,:A2=>Char)) |> DataFrame
    df.A1 = tochar.(df.A1)
    df.A2 = tochar.(df.A2)
    @assert df.CHR isa Vector{Int}
    filter!(row->row.CHR in chroms, df)
    return df
end
parse_bim(bim_file::String, chrom::Integer) =
    parse_bim(bim_file, chrom:chrom)

nuc_map(char::Char) = nuc_map(Val(char))
nuc_map(::Val{'A'}) = 'T'
nuc_map(::Val{'T'}) = 'A'
nuc_map(::Val{'C'}) = 'G'
nuc_map(::Val{'G'}) = 'C'
function permute_snps(df)
    unique(vcat(
        df[:,[:SNP,:A1,:A2]],
        DataFrame(SNP=df.SNP,
                  A1=df.A2,
                  A2=df.A1),
        DataFrame(SNP=df.SNP,
                  A1=nuc_map.(first.(df.A1)),
                  A2=nuc_map.(first.(df.A2))),
        DataFrame(SNP=df.SNP,
                  A1=nuc_map.(first.(df.A2)),
                  A2=nuc_map.(first.(df.A1)))
    ))
end
function join_snps(ref_df, vld_df, sst_df; verbose=false)
    # TODO: Be more efficient, don't allocate all this memory
    vld_snps = vld_df[:,[:SNP,:A1,:A2]]
    ref_snps = permute_snps(ref_df)
    sst_snps = permute_snps(sst_df)
    snps = innerjoin(vld_snps, ref_snps, sst_snps, on=[:SNP,:A1,:A2], makeunique=true)
    verbose && @info "$(nrow(snps)) common SNPs"
    return snps
end
norm_ppf(x) = quantile(Normal(), x)
function parse_sumstats(ref_df, vld_df, sst_file, n_subj; verbose=false, missingstring="")
    sst_df = CSV.File(sst_file; missingstring=missingstring#=, types=Dict(:A1=>Char,:A2=>Char)=#) |> DataFrame
    sst_df.A1 = tochar.(sst_df.A1)
    sst_df.A2 = tochar.(sst_df.A2)
    nucs = Set(['A','C','T','G'])
    filter!(row->(row.A1 in nucs) && (row.A2 in nucs), sst_df)
    if hasproperty(sst_df, :BETA)
        filter!(row->!(row.P isa Missing) && !(row.BETA isa Missing), sst_df)
        sst_df.BETA = convert(Vector{Float64}, sst_df.BETA)
    elseif hasproperty(sst_df, :OR)
        filter!(row->!(row.P isa Missing) && !(row.OR isa Missing), sst_df)
        sst_df.BETA = log.(convert(Vector{Float64}, sst_df.OR))
        select!(sst_df, Not(:OR))
    end
    sst_df.P = convert(Vector{Float64}, sst_df.P)
    snps = join_snps(ref_df, vld_df, sst_df; verbose=verbose)
    sort!(snps, [:SNP, :A1, :A2])

    n_sqrt = sqrt(n_subj)
    sst_eff = Dict{String,Float64}()
    for row in Tables.namedtupleiterator(sst_df)
        if hassnp(snps, (row.SNP,row.A1,row.A2)) ||
           hassnp(snps, (row.SNP,nuc_map.(row.A1),nuc_map.(row.A2)))
            effect_sign = 1
        elseif hassnp(snps, (row.SNP,row.A2,row.A1)) ||
               hassnp(snps, (row.SNP,nuc_map.(row.A2),nuc_map.(row.A1)))
            effect_sign = -1
        else
            continue
        end
        #if hasproperty(row, :BETA)
        #    beta = row.BETA
        #elseif hasproperty(row, :OR)
        #    beta = log(row.OR)
        #end
        p = max(row.P, 1e-323)
        beta_std = effect_sign*sign(row.BETA)*abs(norm_ppf(p/2))/n_sqrt
        sst_eff[row.SNP] = beta_std
    end
    _sst_df = DataFrame(SNP=String[],CHR=Int[],BP=Int[],BETA=Float64[],A1=Char[],A2=Char[],MAF=Float64[],FLP=Int[])
    for (idx,row) in enumerate(Tables.namedtupleiterator(ref_df))
        haskey(sst_eff, row.SNP) || continue

        SNP = row.SNP
        CHR = row.CHR
        BP = row.BP
        BETA = sst_eff[row.SNP]
        A1,A2 = row.A1,row.A2
        if hassnp(snps, (SNP,A1,A2))
            MAF = row.MAF
            FLP = 1
        elseif hassnp(snps, (SNP,A2,A1))
            A1, A2 = A2, A1
            MAF = 1-row.MAF
            FLP = -1
        elseif hassnp(snps, (SNP,nuc_map(A1),nuc_map(A2)))
            A1, A2 = nuc_map(A1), nuc_map(A2)
            MAF = row.MAF
            FLP = 1
        elseif hassnp(snps, (SNP,nuc_map(A2),nuc_map(A1)))
            A1, A2 = nuc_map(A2), nuc_map(A1)
            MAF = 1-row.MAF
            FLP = -1
        else
            verbose && @warn "(Chromosome $CHR) Didn't find ($SNP,$A1,$A2) in snps"
            # FIXME: Skip?
        end
        push!(_sst_df, (SNP=SNP,CHR=CHR,BP=BP,BETA=BETA,A1=A1,A2=A2,MAF=MAF,FLP=FLP))
    end
    return _sst_df
end

function findsnp(snps, (snp,a1,a2))
    SNP_range = binary_range_search(snps, snp, :SNP)
    SNP_range === nothing && return nothing
    SNP_L, SNP_R = SNP_range
    SNP_sub = snps[SNP_L:SNP_R,:]

    A1_range = binary_range_search(SNP_sub, a1, :A1)
    A1_range === nothing && return nothing
    A1_L, A1_R = A1_range
    A1_sub = SNP_sub[A1_L:A1_R,:]

    A2_range = binary_range_search(A1_sub, a2, :A2)
    A2_range === nothing && return nothing
    A2_L, A2_R = A2_range
    @assert A2_L == A2_R
    return SNP_L + (A1_L-1) + (A2_L-1)
end
hassnp(snps, row) = findsnp(snps, row) !== nothing
# TODO: Allow warm-restarts
function binary_range_search(snps, x, col)
    _snps = snps[!,col]
    L = 1
    R = nrow(snps)
    while true
        (L > R) && return nothing
        M = floor(Int, (L+R)/2)
        _x = _snps[M]
        if _x == x
            L,R = M,M
            snps_rows = nrow(snps)
            while (L > 1) && (_snps[L - 1] == x)
                L -= 1
            end
            while R < (snps_rows) && (_snps[R + 1] == x)
                R += 1
            end
            return L,R
        elseif _x < x
            L = M+1
        elseif _x > x
            R = M-1
        end
    end
end

function parse_ldblk(ldblk_dir, sst_df, chrom)
    chr_name = ldblk_dir * "/ldblk_1kg_chr" * string(chrom) * ".hdf5"
    hdf_chr = h5open(chr_name, "r")
    n_blk = length(hdf_chr)
    ld_blk = [read(hdf_chr["blk_"*string(blk)]["ldblk"]) for blk in 1:n_blk]

    snp_blk = Vector{String}[]
    for blk in 1:n_blk
        push!(snp_blk, read(hdf_chr["blk_"*string(blk)]["snplist"]))
    end

    blk_size = Int[]
    _ld_blk = Matrix{Float64}[]
    mm = 1
    for blk in 1:n_blk
        idx = [(ii,findfirst(s->s==snp, sst_df.SNP)) for (ii, snp) in enumerate(snp_blk[blk]) if snp in sst_df.SNP]
        push!(blk_size, length(idx))
        if !isempty(idx)
            idx_blk = mm:(mm+length(snp_blk[blk])-1)
            flip = [sst_df.FLP[jj] for jj in last.(idx)]
            flipM = flip' .* flip
            _blk = Matrix{Float64}(undef, length(idx), length(idx))
            P = collect(Iterators.product(first.(idx),first.(idx)))
            for icol in 1:size(P,2)
                for irow in 1:size(P,1)
                    row,col = P[irow,icol]
                    _blk[irow,icol] = ld_blk[blk][row,col] * flipM[irow,icol]
                end
            end
            push!(_ld_blk, _blk)
            mm += length(snp_blk[blk])
        else
            push!(_ld_blk,  Matrix{Float64}(undef, 0, 0))
        end
    end

    return _ld_blk, blk_size
end

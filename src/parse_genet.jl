function parse_ref(ref_file, chrom)
    println("Parsing reference file: $ref_file")
    df = CSV.File(ref_file) |> DataFrame
    # TODO: Clean up header
    @assert df.BP isa Vector{Int}
    @assert df.MAF isa Vector{Int}
    # TODO: Filter on chrom
    println("$(nrow(df)) SNPs on chromosome $chrom in reference file")
    return df
end
function parse_bim(bim_file, chrom)
    println("Parsing BIM file: $(bim_file*".bim")")
    df = CSV.File(bim_file*".bim") |> DataFrame
    # TODO: Clean up header
    # TODO: Filter on chrom
    println("$(nrow(df)) SNPs in BIM file")
    return df
end
const NUC_MAPPING = Dict('A'=>'T','T'=>'A','C'=>'G','G'=>'C')
function permute_snps(df)
    vcat(
        df[!,[:SNP,:A1,:A2]],
        df[!,[:SNP,:A2,:A1]],
        @with(df[!,[:SNP,:A1,:A2]],
            DataFrame(SNP=:SNP,
                      A1=getindex.(Ref(NUC_MAPPING), :A1),
                      A2=getindex.(Ref(NUC_MAPPING), :A2))),
        @with(df[!,[:SNP,:A1,:A2]],
            DataFrame(SNP=:SNP,
                      A1=getindex.(Ref(NUC_MAPPING), :A1),
                      A2=getindex.(Ref(NUC_MAPPING), :A2)))
    )
end
function join_snps(ref_df, vld_df, sst_df)
    # TODO: Be more efficient, don't allocate all this memory
    vld_snps = vld_df[!,[:SNP,:A1,:A2]]
    ref_snps = permute_snps(ref_df)
    sst_snps = permute_snps(sst_df)
    snps = join(vld_snps, ref_snps, sst_snps; kind=:inner)
    println("$(nrow(snps)) common SNPs")
    return snps
end
function findsnp(df, snp)
    for (idx,row) in enumerate(Tables.rows(df))
        if Tuple(row) == SNP
            return idx
        elseif (row.SNP,NUC_MAPPING[row.A1],NUC_MAPPING[row.A1]) == SNP
            return idx
        end
    end
    return nothing
end
function parse_sumstats(ref_df, vld_df, sst_file, n_subj)
    println("Parsing summary statistics file: $sst_file")
    sst_df = CSV.File(sst_file) |> DataFrame
    # TODO: Clean up header
    @assert df.BETA isa Vector{Int}
    snps = join_snps(ref_df, vld_df, sst_df)

    n_sqrt = sqrt(n_subj)
    sst_eff = similar(sst_df, 0)
    for row in Tables.rows(sst_df)
        snp_rowidx = findsnp(snps, Tuple(row))
        if snp_rowidx !== nothing
            effect_sign = 1
        else
            snp_rowidx_flip = findsnp(snps, (row.SNP,row.A2,row.A1))
            snp_rowidx_flip === nothing && continue
            effect_sign = -1
        end
        if hasproperty(row, :BETA)
            beta = row.BETA
        elseif hasproperty(row, :OR)
            beta = log(row.OR)
        end
        p = max(row.P, 1e-323)
        # FIXME: norm_ppf
        beta_std = effect_sign*sign(beta)*abs(norm_ppf(p/2))/n_sqrt
        push!(sst_eff, merge(collect(row), (BETA_STD=beta_std,)))
    end
    _sst_df = similar(sst_df, 0)
    _sst_df.FLP = Int[]
    for (idx,row) in enumerate(Tables.rows(ref_df))
        snp_rowidx = findfirst(snp->snp==row.SNP, sst_eff.SNP)
        snp_rowidx === nothing && continue

        SNP = row.SNP
        CHR = row.CHR
        BP = row.BP
        BETA = sst_eff[snp_rowidx,:BETA_STD]
        A1,A2 = row.A1,row.A2
        if hassnp(snps, (SNP,A1,A2))
            MAF = row.MAF
            FLP = 1
        elseif hassnp(snps, (SNP,A2,A1))
            A1, A2 = A2, A1
            MAF = 1-row.MAF
            FLP = -1
        elseif hassnp(snps, (SNP,NUC_MAPPING[A1],NUC_MAPPING[A2]))
            A1 = NUC_MAPPING[A1]
            A2 = NUC_MAPPING[A2]
            MAF = row.MAF
            FLP = 1
        elseif hassnp(snps, (SNP,NUC_MAPPING[A2],NUC_MAPPING[A1]))
            A1 = NUC_MAPPING[A2]
            A2 = NUC_MAPPING[A1]
            MAF = 1-row.MAF
            FLP = -11
        end
        push!(_sst_df, (SNP=SNP,CHR=CHR,BP=BP,BETA=BETA,A1=A1,A2=A2,MAF=MAF,FLP=FLP))
    end
    println("$(nrow(_sst_df)) SNPs in summary statistics file")
    return _sst_df
end
function parse_ldblk(ldblk_dir, sst_df, chrom)
    println("Parsing reference LD on chromosome $chrom")
    error("Not yet implemented")
end

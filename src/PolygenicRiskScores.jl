module PolygenicRiskScores

## CSV parsing

using CSV
using DataFrames
using Dates, Distributions, Statistics, Random, LinearAlgebra, Printf
using HDF5
using Profile

include("parse_genet.jl")
include("gigrnd.jl")
include("mcmc_gtb.jl")

## argument parsing

using ArgParse

settings = ArgParseSettings()

@add_arg_table! settings begin
    "--ref_dir"
        help = "Path to the reference panel directory"
        required = true
    "--bim_prefix"
        help = "Directory and prefix of the bim file for the validation set"
        required = true
    "--sst_file"
        help = "Path to summary statistics file"
        required = true
    "--sst_missing"
        help = "Indicator for missing data in sumstats file (eg 'NA')"
        default = ""
    "--a"
        arg_type = Float64
        default = 1.0
    "--b"
        arg_type = Float64
        default = 0.5
    "--phi"
        arg_type = Float64
    "--n_gwas"
        help = "Sample size of the GWAS"
        #arg_type = Int
        required = true
    "--pop"
        help = "Population of the GWAS Sample"
        required = false
        default = nothing
    "--n_iter"
        help = "Number of MCMC iterations to perform"
        arg_type = Int
        default = 1000
    "--n_burnin"
        help = "Number of MCMC burn-in iterations"
        arg_type = Int
        default = 500
    "--thin"
        arg_type = Int
        default = 5
    "--out_dir"
        help = "Output file directory"
        required = true
    "--out_header"
        help = "Write header to output file"
        action = :store_true
    "--out_delim"
        help = "Output file delimiter"
        default = '\t'
    "--out_path"
        help = "Output file path (overrides --out_dir)"
        default = nothing
    "--out_name"
        help = "Output file prefix"
        default = nothing
    "--chrom"
        help = "Chromosomes to process"
        default = "1:22"
    "--beta_std"
        action = :store_true
    "--meta"
        help = "If true, return combined SNP effect sizes across populations using an inverse-variance-weighted meta-analysis of the population-specific posterior effect size estimates."
        default = false
        arg_type = Bool
    "--seed"
        help = "RNG seed for MCMC"
        arg_type = Int
    "--quiet"
        help = "Disable all unnecessary printing"
        action = :store_true
    "--hostsfile"
        help = "Hostsfile to use for parallel processing"
        default = nothing
    "--profile"
        help = "Enables output of profiling data of the MCMC"
        action = :store_true
end

function main()
    opts = parse_args(ARGS, settings)
    verbose = !opts["quiet"]

    chroms = eval(Meta.parse(opts["chrom"]))
    verbose && @info "Selecting chromosomes $chroms"

    ref_dir = opts["ref_dir"]
    verbose && @info "Parsing reference file: $ref_dir/snpinfo_1kg_hm3"
    t = now()
    if opts["pop"] === nothing
        ref_df = parse_ref(joinpath(ref_dir, "snpinfo_1kg_hm3"), chroms)
    else
        ref_df = parse_ref(joinpath(ref_dir, "snpinfo_mult_1kg_hm3"), chroms; multi=true)
    end
    verbose && @info "$(nrow(ref_df)) SNPs in reference file ($(round(now()-t, Dates.Second)))"

    bim_prefix = opts["bim_prefix"]
    verbose && @info "Parsing BIM file: $(bim_prefix*".bim")"
    t = now()
    vld_df = parse_bim(bim_prefix, chroms)
    verbose && @info "$(nrow(vld_df)) SNPs in BIM file ($(round(now()-t, Dates.Second)))"

    for chrom in chroms
        _main(chrom, ref_df, vld_df, opts; verbose=verbose)
    end
end
function _main(chrom, ref_df, vld_df, opts; verbose=false)
    sst_files = split(opts["sst_file"], ',')
    n_gwass = parse.(Int, split(opts["n_gwas"], ','))
    if opts["pop"] != nothing
        pops = split(opts["pop"], ',')
        n_pop = length(pops)
    else
        pops = [nothing]
        n_pop = 1
    end
    meta = opts["meta"]

    sst_dfs = Vector{DataFrame}(undef, length(n_gwass))
    ld_blks = Vector{Vector{Matrix{Float64}}}(undef, length(n_gwass))
    blk_sizes = Vector{Vector{Int}}(undef, length(n_gwass))
    for i in 1:length(n_gwass)
        sst_file = sst_files[i]
        pop = pops[i]
        verbose && @info "(Chromosome $chrom) (Population $pop) Parsing summary statistics file: $sst_file"
        t = now()
        sst_df = parse_sumstats(sst_file; verbose=verbose, missingstring=opts["sst_missing"])
        sst_eff = extract_effect_alleles(sst_df, ngwass[i])
        _ref_df = ref_df[ref_df.CHR .== chrom,:]
        _vld_df = vld_df[vld_df.CHR .== chrom,:]
        snps = join_snps(_ref_df, _vld_df, sst_df; verbose=verbose)
        sst_df = compute_frequencies(_ref_df,snps, sst_eff, pop; verbose=verbose)
        sst_dfs[i] = sst_df
        verbose && @info "(Chromosome $chrom) (Population $pop) $(nrow(sst_df)) SNPs in summary statistics file ($(round(now()-t, Dates.Second)))"

        verbose && @info "(Chromosome $chrom) (Population $pop) Parsing reference LD"
        t = now()
        snp_blk, n_blk, ld_blk = parse_ldblk(opts["ref_dir"], sst_df, chrom, pop)
        ld_blks[i], blk_sizes[i] = shrink_ldblk(snp_blk, n_blk, ld_blk, sst_df)
        verbose && @info "(Chromosome $chrom) (Population $pop) Completed parsing reference LD ($(round(now()-t, Dates.Second)))"
    end

    verbose && @info "(Chromosome $chrom) (Population $pop) Aligning LD blocks"
    t = now()
    snp_df, beta_vecs, frq_vecs, idx_vecs = align_ldblk(ref_df, vld_df, sst_dfs, length(n_gwass), chrom)
    verbose && @info "(Chromosome $chrom) (Population $pop) Aligned LD blocks ($(round(now()-t, Dates.Second)))"

    verbose && @info "(Chromosome $chrom) Initiating MCMC"
    t = now()
    beta_est, extra = mcmc(a=opts["a"], b=opts["b"], phi=opts["phi"], snp_df=snp_df, beta_vecs=beta_vecs, frq_vecs=frq_vecs, idx_vecs=idx_vecs, sst_df=sst_dfs, n=n_gwass, ld_blk=ld_blks, blk_size=blk_sizes, n_iter=opts["n_iter"], n_burnin=opts["n_burnin"], thin=opts["thin"], chrom=chrom, beta_std=opts["beta_std"], meta=meta, seed=opts["seed"], verbose=verbose, profile=opts["profile"])
    verbose && @info "(Chromosome $chrom) Completed MCMC ($(round(now()-t, Dates.Second)))"

    phi = opts["phi"]
    phi_str = phi === nothing ? "auto" : @sprintf("%1.0e", phi)

    out_dir = opts["out_path"] !== nothing ? dirname(opts["out_path"]) : opts["out_dir"]
    out_prefix = opts["out_name"] !== nothing ? (opts["out_name"] * "_") : ""

    for pp in 1:n_pop
        pop = pops[pp]
        verbose && @info "(Chromosome $chrom) (Population $pop) Writing posterior effect sizes"
        t = now()

        pop_str = n_pop == 1 ? "" : (pops[pp] * "_")
        eff_file = if opts["out_path"] === nothing
            if !isdir(out_dir)
                @warn "--out_dir does not exist; creating it at $(opts["out_dir"])\nNote: The old behavior of treating --out_dir as a prefix has been removed"
                mkdir(out_dir)
            end
            joinpath(out_dir, @sprintf("%s%spst_eff_a%d_b%.1f_phi%s_chr%d.txt", out_prefix, pop_str, opts["a"], opts["b"], phi_str, chrom))
        else
            @warn "Prefixing output file with $pop_str"
            joinpath(out_dir, pop_str * basename(opts["out_path"]))
        end

        out_df = snp_df[idx_vecs[pp], [:SNP, :BP, :A1, :A2]]
        out_df[!, :CHR] .= chrom
        out_df.BETA = map(b->@sprintf("%.6e", b), beta_est[pp])
        out_df = select(out_df, [:CHR, :SNP, :BP, :A1, :A2, :BETA])
        CSV.write(eff_file, out_df; header=opts["out_header"], delim=opts["out_delim"])
        verbose && @info "(Chromosome $chrom) (Population $pop) Finished writing posterior effect sizes ($(round(now()-t, Dates.Second)))"
    end

    if meta
        verbose && @info "(Chromosome $chrom) Writing meta posterior effect sizes"
        t = now()
        meta_eff_file = if opts["out_path"] === nothing
            joinpath(out_dir, @sprintf("_META_pst_eff_a%d_b%.1f_phi%s_chr%d.txt", opts["a"], opts["b"], phi_str, chrom))
        else
            @warn "Prefixing meta output file with META_"
            joinpath(out_dir, "META_" * basename(opts["out_path"]))
        end

        mu = extra.mu
        @assert mu !== nothing
        meta_out_df = snp_df[:, [:SNP, :BP, :A1, :A2]]
        meta_out_df[!, :CHR] .= chrom
        meta_out_df.BETA = map(b->@sprintf("%.6e", b), mu)
        meta_out_df = select(meta_out_df, [:CHR, :SNP, :BP, :A1, :A2, :BETA])
        CSV.write(meta_eff_file, meta_out_df; header=opts["out_header"], delim=opts["out_delim"])
        verbose && @info "(Chromosome $chrom) Finished writing meta posterior effect sizes ($(round(now()-t, Dates.Second)))"
    end

    if opts["profile"]
        prof_path = joinpath(out_dir, out_prefix * "chr$chrom.cpuprofile")
        open(prof_path, "w") do io
            Profile.print(io; C=true, mincount=100, maxdepth=100)
        end
    end

    if opts["phi"] === nothing && verbose
        @info @sprintf("Estimated global shrinkage parameter: %1.2e", extra.phi_est)
    end
end

end # module

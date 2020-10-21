module PRS

## CSV parsing

using CSV
using DataFrames, DataFramesMeta
using Dates, Distributions, Statistics, Random, LinearAlgebra, Printf
using HDF5

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
        arg_type = Int
        required = true
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
        help = "Output file directory and prefix"
        required = true
    "--out_header"
        help = "Write header to output file"
        action = :store_true
    "--out_delim"
        help = "Output file delimiter"
        default = '\t'
    "--out_path"
        help = "Output file path (overrides --out_dir)"
    "--chrom"
        help = "Chromosomes to process"
        default = "1:23"
    "--beta_std"
        action = :store_true
    "--seed"
        help = "RNG seed for MCMC"
        arg_type = Int
    "--quiet"
        help = "Disable all unnecessary printing"
        action = :store_true
    "--hostsfile"
        help = "Hostsfile to use for parallel processing"
        default = nothing
end

function main()
    opts = parse_args(ARGS, settings)
    verbose = !opts["quiet"]

    chroms = eval(Meta.parse(opts["chrom"]))
    verbose && @info "Selecting chromosomes $chroms"

    ref_dir = opts["ref_dir"]
    verbose && @info "Parsing reference file: $ref_dir/snpinfo_1kg_hm3"
    t = now()
    ref_df = parse_ref(ref_dir * "/snpinfo_1kg_hm3", chroms)
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
    sst_file = opts["sst_file"]
    verbose && @info "(Chromosome $chrom) Parsing summary statistics file: $sst_file"
    t = now()
    sst_df = parse_sumstats(ref_df[ref_df.CHR .== chrom,:], vld_df[vld_df.CHR .== chrom,:], sst_file, opts["n_gwas"]; verbose=verbose, missingstring=opts["sst_missing"])
    verbose && @info "(Chromosome $chrom) $(nrow(sst_df)) SNPs in summary statistics file ($(round(now()-t, Dates.Second)))"

    verbose && @info "(Chromosome $chrom) Parsing reference LD"
    t = now()
    ld_blk, blk_size = parse_ldblk(opts["ref_dir"], sst_df, chrom)
    verbose && @info "(Chromosome $chrom) Completed parsing reference LD ($(round(now()-t, Dates.Second)))"

    verbose && @info "(Chromosome $chrom) Initiating MCMC"
    t = now()
    beta_est = mcmc(opts["a"], opts["b"], opts["phi"], sst_df, opts["n_gwas"], ld_blk, blk_size, opts["n_iter"], opts["n_burnin"], opts["thin"], chrom, opts["beta_std"], opts["seed"]; verbose=verbose)
    verbose && @info "(Chromosome $chrom) Completed MCMC ($(round(now()-t, Dates.Second)))"

    verbose && @info "(Chromosome $chrom) Writing posterior effect sizes"
    eff_file = if opts["out_path"] === nothing
        out_path = opts["out_dir"]
        phi = opts["phi"]
        phi_str = phi === nothing ? "auto" : @sprintf("%1.0e", phi)
        out_path * @sprintf("_pst_eff_a%d_b%.1f_phi%s_chr%d.txt", opts["a"], opts["b"], phi_str, chrom)
    else
        opts["out_path"]
    end
    t = now()
    out_df = sst_df[:, [:SNP, :BP, :A1, :A2]]
    out_df[!, :CHR] .= chrom
    out_df.BETA = map(b->@sprintf("%.6e", b), beta_est)
    out_df = select(out_df, [:CHR, :SNP, :BP, :A1, :A2, :BETA])
    CSV.write(eff_file, out_df; header=opts["out_header"], delim=opts["out_delim"])
    verbose && @info "(Chromosome $chrom) finished writing posterior effect sizes ($(round(now()-t, Dates.Second)))"
end

end # module

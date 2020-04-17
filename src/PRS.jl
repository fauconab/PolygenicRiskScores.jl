module PRS

## CSV parsing

using CSV

function parse_ref(ref_file, chrom)
    println("Parsing reference file: $ref_file")
    df = CSV.File(ref_file) |> DataFrame
    # TODO: Assert header
    @assert df.BP isa Vector{Int}
    @assert df.MAF isa Vector{Int}
    println("$(nrow(df)) SNPs on chromosome $chrom in reference file")
    return df
end
function parse_bim(bim_file, chrom)
    println("Parsing BIM file: $(bim_file*".bim")")
    df = CSV.File(bim_file*".bim") |> DataFrame
    # TODO: Assert header
    println("$(nrow(df)) SNPs in BIM file")
    return df
end
function parse_sumstats(ref_df, vld_df, sst_file, n_subj)
    println("Parsing summary statistics file: $sst_file")
    df = CSV.File(sst_file) |> DataFrame
    # TODO: Assert header
    @assert df.BETA isa Vector{Int}
    println("$(nrow(df)) SNPs in summary statistics file")
    return df
end

## argument parsing

using ArgParse

settings = ArgParseSettings()

@add_arg_table! settings begin
    "--ref_dir"
        help="Path to the reference panel directory"
        required=true
    "--bim_prefix"
        help="Directory and prefix of the bim file for the validation set"
        required=true
    "--sst_file"
        help="Path to summary statistics file"
        required=true
    "--a"
        arg_type=Float64
        default=1.0
    "--b"
        arg_type=Float64
        default=0.5
    "--phi"
        arg_type=Float64
    "--n_gwas"
        help="Sample size of the GWAS"
        arg_type=Int
        required=true
    "--n_iter"
        arg_type=Int
        default=1000
    "--n_burnin"
        arg_type=Int
        default=500
    "--thin"
        arg_type=Int
        default=5
    "--out_dir"
        help="Output directory path"
        required=true
    "--chrom"
        default=1:23
    "--beta_std"
        action = :store_true
        default=false
    "--seed"
        arg_type=Int
end

function main()
    opts = parse_args(ARGS, settings)
    for chrom in opts["chrom"]
        ref_dict = parse_ref(opts["ref_dir"] * "/snpinfo_1kg_hm3", chrom)
        vld_dict = parse_bim(opts["bim_prefix"], chrom)
        sst_dict = parse_sumstats(ref_dict, vld_dict, opts["sst_file"], opts["n_gwas"])
        ld_blk, blk_size = parse_ldblk(opts["ref_dir"], sst_dict, chrom)
        mcmc(opts["a"], opts["b"], opts["phi"], sst_dict, opts["n_gwas"], ld_blk, blk_size, opts["n_iter"], opts["n_burnin"], opts["thin"], chrom, opts["out_dir"], opts["beta_std"], opts["seed"])
    end
end

end # module

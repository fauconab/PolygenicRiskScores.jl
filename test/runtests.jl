using Test
using Printf
using CSV, DataFrames, Tables
using HypothesisTests

PRS_DATA_PATH = get(ENV, "JULIA_PRS_TEST_DATA_PATH", nothing)
if PRS_DATA_PATH === nothing
    PRS_DATA_PATH = mktempdir()
end
cd(PRS_DATA_PATH) do
    isdir("PRScs") || run(`git clone https://github.com/getian107/PRScs`)
    isdir("PRScsx") || run(`git clone https://github.com/getian107/PRScsx`)
    isfile("test.bim") || run(`ln -s PRScs/test_data/test.bim .`)
    isfile("multi_test.bim") || run(`ln -s PRScsx/test_data/test.bim ./multi_test.bim`)
    isfile("sumstats.txt") || run(`ln -s PRScs/test_data/sumstats.txt .`)
    isfile("EUR_sumstats.txt") || run(`ln -s PRScsx/test_data/EUR_sumstats.txt .`)
    isfile("EAS_sumstats.txt") || run(`ln -s PRScsx/test_data/EAS_sumstats.txt .`)
    isfile("snpinfo_mult_1kg_hm3") || run(`wget -q https://www.dropbox.com/s/rhi806sstvppzzz/snpinfo_mult_1kg_hm3`)
    if !isdir("ldblk_1kg_eur")
        run(`wget -q https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz`)
        run(`tar xf ldblk_1kg_eur.tar.gz`)
    end
    if !isdir("ldblk_1kg_eas")
        run(`wget -q https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eas.tar.gz`)
        run(`tar xf ldblk_1kg_eas.tar.gz`)
    end
end

#These tests check the following:
#Script runs using phi auto version
#Script runs using specified fi
#Script successfully saves results using out_path and with out_dir

function test_harness(;a=1,b=0.5,phi=1e-02,chr=22,n_iter=1000,n_burnin=500,out_header=false,multi=false,meta=false,out_name="")
    a = repr(a)
    b = repr(b)
    phi = @sprintf("%1.0e", phi)
    chr = repr(chr)
    n_iter = repr(n_iter)
    bim_prefix = multi ? "multi_test" : "test"

    cmd_jl = `julia --project -e "using PolygenicRiskScores; PolygenicRiskScores.main()" -- --bim_prefix=$PRS_DATA_PATH/$bim_prefix --chrom=$chr --phi=$phi --n_iter=$n_iter --n_burnin=$n_burnin --out_dir=$PRS_DATA_PATH/output_jl`
    py_pkg = multi ? "PRScsx" : "PRScs"
    cmd_py = `python3 $PRS_DATA_PATH/$py_pkg/$py_pkg.py --bim_prefix=$PRS_DATA_PATH/$bim_prefix --chrom=$chr --phi=$phi --n_iter=$n_iter --n_burnin=$n_burnin --out_dir=$PRS_DATA_PATH/output_py`

    for (kind, cmd) in ((:jl, cmd_jl), (:py, cmd_py))
        if out_header && kind == :jl
            push!(cmd.exec, "--out_header")
        end
        if multi
            push!(cmd.exec, "--meta=$meta")
            push!(cmd.exec, "--sst_file=$PRS_DATA_PATH/EUR_sumstats.txt,$PRS_DATA_PATH/EAS_sumstats.txt")
            push!(cmd.exec, "--n_gwas=200000,100000")
            push!(cmd.exec, "--pop=EUR,EAS")
            push!(cmd.exec, "--ref_dir=$PRS_DATA_PATH")
            if kind == :py
                out_name = "py_test"
                mkpath(joinpath(PRS_DATA_PATH,"output_py"))
            end
            if out_name != ""
                push!(cmd.exec, "--out_name=$out_name")
            end
        else
            push!(cmd.exec, "--sst_file=$PRS_DATA_PATH/sumstats.txt")
            push!(cmd.exec, "--n_gwas=200000")
            push!(cmd.exec, "--ref_dir=$PRS_DATA_PATH/ldblk_1kg_eur")
        end
        run(cmd)
    end

    pops = multi ? ["EUR","EAS"] : [nothing]
    for pop in pops
        pop_str=multi ? pop*"_" : ""
        output_jl = joinpath(PRS_DATA_PATH, "output_jl", pop_str * "pst_eff_a$(a)_b$(b)_phi$(phi)_chr$(chr).txt")
        if multi
            out_prefix = out_name != "" ? out_name : "py_test"
            output_py = joinpath(PRS_DATA_PATH, "output_py", out_prefix * "_" * pop_str * "pst_eff_a$(a)_b$(b)_phi$(phi)_chr$(chr).txt")
        else
            output_py = joinpath(PRS_DATA_PATH, "output_py_$(pop_str)pst_eff_a$(a)_b$(b)_phi$(phi)_chr$(chr).txt")
        end

        beta_jl = nothing
        beta_py = nothing

        for (kind, output_path) in ((:jl, output_jl), (:py, output_py))
            @show output_path
            @test isfile(output_path)
            @test stat(output_path).size > 30_000
            header = [:CHR, :SNP, :BP, :A1, :A2, :BETA]
            if out_header
                f = CSV.File(output_path)
            else
                f = CSV.File(output_path; header=header)
            end
            for (idx,col) in enumerate(header)
                if out_header && kind == :jl
                    @test Tables.columnnames[idx] == col
                end
                if !multi
                    @test length(getproperty(f, col)) == 1000
                end
            end
            if kind == :jl
                beta_jl = f.BETA
            else
                beta_py = f.BETA
            end
        end

        # FIXME: Add accuracy/comparison tests
        # Do a t-test to make sure the beta results are not significantly different
        @test pvalue(EqualVarianceTTest(beta_jl,beta_py)) > 0.01

        rm(joinpath(PRS_DATA_PATH, output_jl))
        rm(joinpath(PRS_DATA_PATH, output_py))
    end
end

@testset "Single Ancestry" begin
    # Default run
    test_harness()
end
@testset "Multi Ancestry" begin
    # Default run
    test_harness(;multi=true, n_iter=10000, n_burnin=500)
end

# run julia once with seed X
# get hash of output file
# run julia again with seed X
# test that hash of new output file == old output file

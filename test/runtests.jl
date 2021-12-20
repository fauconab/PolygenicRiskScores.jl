using Test
using Printf
using CSV, DataFrames, Tables
using HypothesisTests

PRS_DATA_PATH = get(ENV, "JULIA_PRS_TEST_DATA_PATH", nothing)
if PRS_DATA_PATH === nothing
    PRS_DATA_PATH = mktempdir()
    cd(PRS_DATA_PATH) do
        isdir("PRScs") || run(`git clone https://github.com/getian107/PRScs`)
        isfile("test.bim") || run(`ln -s PRScs/test_data/test.bim .`)
        isfile("sumstats.txt") || run(`ln -s PRScs/test_data/sumstats.txt .`)
        if !isdir("ldblk_1kg_eur")
            run(`wget https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz`)
            run(`tar xvf ldblk_1kg_eur.tar.gz`)
        end
    end
end

#These tests check the following:
#Script runs using phi auto version
#Script runs using specified fi
#Script successfully saves results using out_path and with out_dir

function test_harness(;a=1,b=0.5,phi=1e-02,chr=22,niter=1000,out_header=false)
    a = repr(a)
    b = repr(b)
    phi = @sprintf("%1.0e", phi)
    chr = repr(chr)
    niter = repr(niter)

    cmd_jl = `julia --project -e "using PolygenicRiskScores; PolygenicRiskScores.main()" -- --ref_dir=$PRS_DATA_PATH/ldblk_1kg_eur --bim_prefix=$PRS_DATA_PATH/test --sst_file=$PRS_DATA_PATH/sumstats.txt --n_gwas=200000 --chrom=$chr --phi=$phi --n_iter=$niter --out_dir=$PRS_DATA_PATH/output_jl`
    if out_header
        push!(cmd_jl.exec, "--out_header")
    end
    run(cmd_jl)

    run(`python3 $PRS_DATA_PATH/PRScs/PRScs.py --ref_dir=$PRS_DATA_PATH/ldblk_1kg_eur --bim_prefix=$PRS_DATA_PATH/test --sst_file=$PRS_DATA_PATH/sumstats.txt --n_gwas=200000 --chrom=$chr --phi=$phi --n_iter=$niter --out_dir=$PRS_DATA_PATH/output_py`)

    output_jl = "output_jl_pst_eff_a$(a)_b$(b)_phi$(phi)_chr$(chr).txt"
    output_py = "output_py_pst_eff_a$(a)_b$(b)_phi$(phi)_chr$(chr).txt"

    beta_jl=nothing
    beta_py=nothing

    for (kind, output_path) in ((:jl, output_jl), (:py, output_py))
        output_path = joinpath(PRS_DATA_PATH, output_path)
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
            @test length(getproperty(f, col)) == 1000
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

@testset "Basic Options" begin
    # Default run
    test_harness()
end

# run julia once with seed X
# get hash of output file
# run julia again with seed X
# test that hash of new output file == old output file

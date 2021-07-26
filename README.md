# PolygenicRiskScores.jl (previously PRS.jl)

PolygenicRiskScores.jl is a port of [PRS-CS](https://github.com/getian107/PRScs) to Julia.

## Usage

Using the test data from PRS-CS at ~/prs-data, the following invocation should work when run
in the root directory of PolygenicRiskScores.jl:

```
julia --project -e "using PolygenicRiskScores; PolygenicRiskScores.main()" -- --ref_dir=~/prs-data/ldblk_1kg_eur --bim_prefix=~/prs-data/test --sst_file=~/prs-data/sumstats.txt --n_gwas=200000 --chrom=22 --phi=1e-2 --n_iter=1000 --out_dir=~/prs-data/output_jl
```

## Multi-threaded CSV reading

Julia's CSV reader supports multi-threaded CSV reading. In order to enable this, we need to load Julia with more than one thread. This is done through the `JULIA_NUM_THREADS` environment variable: 

```
JULIA_NUM_THREADS=8 julia --project -e "using PolygenicRiskScores; PolygenicRiskScores.main()" ...
```

The above code would run PolygenicRiskScores.jl with 8 threads.

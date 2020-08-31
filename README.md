# PRS.jl

PRS.jl is a port of [PRS-CS](https://github.com/getian107/PRScs) to Julia.

## Usage

Using the test data from PRS-CS at ~/prs-data, the following invocation should work when run
in the root directory of PRS.jl:

```
julia --project -e "using PRS; PRS.main()" -- --ref_dir=~/prs-data/ldblk_1kg_eur --bim_prefix=~/prs-data/test --sst_file=~/prs-data/sumstats.txt --n_gwas=200000 --chrom=22 --phi=1e-2 --n_iter=1000 --out_dir=~/prs-data/output_jl
```

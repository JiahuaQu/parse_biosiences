https://www.parsebiosciences.com/datasets/10-million-human-pbmcs-in-a-single-experiment/

```bash


bsub -P jhq -J parse -n 1 -q interactive -R "rusage[mem=4001]" -W 10:00 -Is "bash"

cd /research/sharedresources/immunoinformatics/common/jqu/parse_biosciences/datasets

nohup wget --continue --tries=0 -O Parse_10M_PBMC_cytokines.h5ad https://parse-wget.s3.us-west-2.amazonaws.com/10m/Parse_10M_PBMC_cytokines.h5ad > download.log 2>&1 &

md5sum Parse_10M_PBMC_cytokines.h5ad > md5sum.txt
# f850f0447147190c5e5216e0d2e2d162  Parse_10M_PBMC_cytokines.h5ad

# The MD5 value from the company:
# f850f0447147190c5e5216e0d2e2d162  Parse_10M_PBMC_cytokines.h5ad

```




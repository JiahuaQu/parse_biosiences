#!/bin/bash
#BSUB -P anndata
#BSUB -J transferh5ad
#BSUB -q priority
#BSUB -n 20
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=50000]"
#BSUB -oo run_%J.log
#BSUB -eo run_%J.error
#BSUB -W 24:00
#BSUB -B jiahua.qu@stjude.org
#BSUB -N jiahua.qu@stjude.org
#BSUB -notify start done exit suspend

module load mamba/1.4.2
mamba activate anndata_env

python DEG_post_pca.py
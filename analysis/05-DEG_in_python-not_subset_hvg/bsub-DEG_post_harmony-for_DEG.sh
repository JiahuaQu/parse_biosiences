#!/bin/bash
#BSUB -P anndata
#BSUB -J parse_DEG
#BSUB -q superdome
#BSUB -n 30
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=50000]"
#BSUB -oo run_%J.log
#BSUB -eo run_%J.error
#BSUB -W 100:00
#BSUB -B jiahua.qu@stjude.org
#BSUB -N jiahua.qu@stjude.org
#BSUB -notify "start done exit suspend"

module load mamba/1.4.2
mamba activate anndata_env

python DEG_post_harmony-for_DEG.py
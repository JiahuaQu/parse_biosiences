#!/bin/bash
#BSUB -P anndata
#BSUB -J parse
#BSUB -q superdome
#BSUB -n 20
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=50000]"
#BSUB -oo run_%J.log
#BSUB -eo run_%J.error
#BSUB -W 24:00
#BSUB -B jiahua.qu@stjude.org
#BSUB -N jiahua.qu@stjude.org
#BSUB -notify "start done exit suspend"

module load conda3/202311
conda activate anndata_env

python DEG_post_harmony-for_process_and_UMAP.py
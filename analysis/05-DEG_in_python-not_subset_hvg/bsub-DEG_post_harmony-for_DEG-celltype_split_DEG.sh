#!/bin/bash
#BSUB -P anndata
#BSUB -J split_DEG
#BSUB -q superdome
#BSUB -n 30
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=50000]"
#BSUB -oo run_%J.log
#BSUB -eo run_%J.error
#BSUB -W 100:00
#BSUB -B 
#BSUB -N 

# 使用环境变量传递的参数
echo "Cell type: ${CELL_TYPE}"

module load mamba/1.4.2
mamba activate anndata_env

/home/jqu/.conda/envs/anndata_env/bin/python DEG_post_harmony-for_DEG-celltype_split_DEG.py --cell_types "${CELL_TYPE}"
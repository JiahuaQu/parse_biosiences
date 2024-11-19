#!/bin/bash
#BSUB -P anndata
#BSUB -J cell_type
#BSUB -q superdome
#BSUB -n 5
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=50000]"
#BSUB -oo run_%J.log
#BSUB -eo run_%J.error
#BSUB -W 1:00
#BSUB -B 
#BSUB -N 

module load mamba/1.4.2
mamba activate anndata_env

/home/jqu/.conda/envs/anndata_env/bin/python cell_type_count_2_2.py
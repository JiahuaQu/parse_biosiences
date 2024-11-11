```bash
bsub -P jupyter -J scanpy -n 5 -q gpu -gpu "num=1/host" -R "span[hosts=1]" -R "rusage[mem=100000]" -W 24:00 -Is "bash"

cd /research/sharedresources/immunoinformatics/common/jqu/parse_biosciences/

module load mamba/1.4.2

#mamba create -n anndata_env anndata h5py python=3.9
#mamba env list
# anndata_env     /home/jqu/.conda/envs/anndata_env

mamba activate anndata_env

mamba list | grep scipy
# scipy                     1.13.1           py39haf93ffa_0    conda-forge
mamba update scipy
mamba list | grep scipy
# scipy                     1.13.1           py39haf93ffa_0    conda-forge

#mamba install anndata
#mamba install scanpy scipy numpy pandas seaborn matplotlib
#mamba install harmonypy
mamba install leidenalg
mamba deactivate

launch_jupyterlab.sh
# /home/jqu/.conda/envs/anndata_env
```

```bash
To access the server, open this file in a browser:
        file:///scratch_local/jqu/jupyter_runtime_dir_1/jpserver-1423578-open.html
    Or copy and paste one of these URLs:
        http://10.220.17.184:62162/lab?token=62f464ad2e25f6e5891e3e3dcac8f33f5490cb29db069a67
        http://127.0.0.1:62162/lab?token=62f464ad2e25f6e5891e3e3dcac8f33f5490cb29db069a67


To access the server, open this file in a browser:
        file:///scratch_local/jqu/jupyter_runtime_dir_0/jpserver-1057330-open.html
    Or copy and paste one of these URLs:
        http://10.220.19.182:57075/lab?token=1f1fb7e3f0e997fb39c3079cbe4fb0bb91d69116b2b41697
        http://127.0.0.1:57075/lab?token=1f1fb7e3f0e997fb39c3079cbe4fb0bb91d69116b2b41697

```

```bash
# For qsub
bsub -P qsub -J anndata -n 1 -q interactive -R "rusage[mem=4001]" -W 2:00 -Is "bash"

# Submit the job in this path.
cd /home/jqu/project/parse_biosiences/03-python_transfer_h5seurat
# Specify the read in and output absolute path whose prefix is
# "/research/sharedresources/immunoinformatics/common/jqu/parse_biosciences/datasets/".
bsub < bsub-transfer.sh
```

Process and generate GED in python by use of the script template:

[Parse_10M_PBMC_cytokine_clustering_tutorial](file:///Z:/ResearchHome/ClusterHome/jqu/project/parse_biosiences/01-download/Parse_10M_PBMC_cytokine_clustering_tutorial.html)

```bash
# For qsub
bsub -P qsub -J anndata -n 1 -q interactive -R "rusage[mem=4001]" -W 6:00 -Is "bash"

# Submit the job in this path.
cd /home/jqu/project/parse_biosiences/04-DEG_in_python
bsub < bsub-DEG_post_harmony_post_umap.sh
```


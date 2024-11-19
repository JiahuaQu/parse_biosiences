# 1. Jupyter Lab

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
        file:///scratch_local/jqu/jupyter_runtime_dir_0/jpserver-2468498-open.html
    Or copy and paste one of these URLs:
        http://10.220.17.172:52796/lab?token=dcdf584f0f5ba986a4b52f54013c7356fc9dc2f3a7ad08d6
        http://127.0.0.1:52796/lab?token=dcdf584f0f5ba986a4b52f54013c7356fc9dc2f3a7ad08d6
```

# 2. Bsub

## 2.1. Single python script

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
cd /home/jqu/project/parse_biosiences/05-DEG_in_python-not_subset_hvg
bsub < bsub-DEG_post_harmony_post_umap.sh
```

## 2.2. Split into multiple cell types and run them individually

To expedite the DEG calculation, I split the job into multiple small jobs ran in parallel:

(1) The python script (DEG_post_harmony-for_DEG-celltype_split_DEG.py) to run python. It receives the parameter from the bsub.sh file:

(2) The bsub.sh file (bsub-DEG_post_harmony-for_DEG-celltype_split_DEG.sh) submit a small job associated with one cell_type. It receives the environment variate from another bsub file () and inputs the cell_type as a parameter into the python script (DEG_post_harmony-for_DEG-celltype_split_DEG.py).

```bash
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
# Specify: "/home/jqu/.conda/envs/anndata_env/bin/python"
# Input parameter: --cell_types "${CELL_TYPE}"  
```

(3) Another bsub.sh file (bsub-DEG_post_harmony-for_DEG-celltype_split_DEG-bsub.sh). It reads in cell_type in loop and inputs it into the above bsub.sh file.

There are 18 different cell types.

```bash
#!/bin/bash
#BSUB -P anndata
#BSUB -J bsub_bsub
#BSUB -q superdome
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=5000]"
#BSUB -oo run_%J.log
#BSUB -eo run_%J.error
#BSUB -W 1:00
#BSUB -B 
#BSUB -N 

# 路径设置
analysis_save_path="./"
txt_file="${analysis_save_path}/cell_types.txt"
job_script="bsub-DEG_post_harmony-for_DEG-celltype_split_DEG.sh" # 你的具体 job 脚本路径

# 检查 txt 文件是否存在
if [[ ! -f "${txt_file}" ]]; then
    echo "Error: File '${txt_file}' does not exist or cannot be accessed!"
    exit 1
fi

# 初始化行号计数器
line_num=1

# 逐行读取文件内容
while IFS= read -r cell_type || [[ -n "${cell_type}" ]]; do
    echo "Processing line $line_num: $cell_type"
    
    # 确保读取的 cell_type 不为空
    if [[ -n "${cell_type}" ]]; then
        echo "Submitting job for cell type: ${cell_type}"
        export CELL_TYPE="${cell_type}"
        bsub -env "CELL_TYPE=${cell_type}" < "$job_script"
    fi
    ((line_num++))
done < "$txt_file"
```

# 3. Test script in Jupyter Lab again

```bash
bsub -P jupyter -J scanpy -n 5 -q gpu -gpu "num=1/host" -R "span[hosts=1]" -R "rusage[mem=100000]" -W 24:00 -Is "bash"

cd /home/jqu/project/parse_biosiences

launch_jupyterlab.sh
# /home/jqu/.conda/envs/anndata_env
```

```markdown
http://10.220.17.172:54363/lab?token=b7165f654c58c2b908a2a9b5f23fd5ecab2f529fc5d1f5a2
```

# 4. Submit each cell type individually

```python
import os
import sys
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
import warnings
import pickle
import argparse

# 设置命令行参数解析器
parser = argparse.ArgumentParser(description="Process cell types for DEG analysis.")
parser.add_argument('--cell_type', type=str, help='The cell type to analyze')

# 解析参数
args = parser.parse_args()

# 获取参数值
cell_type = args.cell_type

# 打印参数进行调试
print(f"Cell type: {cell_type}")

sc.settings.verbosity = 2

analysis_save_path = './'
adata = sc.read_h5ad(os.path.join(analysis_save_path, '../05-DEG_in_python-not_subset_hvg/adata_post_harmony.h5ad'))

### DEG
# 用于存储每个 cell_type 的差异基因信息
diff_gene_info = {}

# 子集化数据以包含当前 cell_type 的细胞
adata_sub = adata[adata.obs['cell_type'] == cell_type]
    
# 确保 cytokine 列为分类类型
if 'cytokine' not in adata_sub.obs:
    print(f"'cytokine' column not found for cell type {cell_type}. Skipping...")

# Check if 'cytokine' is already categorical
if not isinstance(adata_sub.obs['cytokine'].dtype, pd.CategoricalDtype):
    # Convert to categorical only if it's not already categorical
    adata_sub.obs['cytokine'] = adata_sub.obs['cytokine'].astype("category")

# 检查 'PBS' 是否在 cytokine 的类别中
if 'PBS' not in adata_sub.obs['cytokine'].unique():
    print(f"'PBS' not found in 'cytokine' for cell type {cell_type}. Skipping...")

# 确保 adata_sub 是一个副本
adata_sub = adata_sub.copy()

# 将 'PBS' 设置为参考组
if 'PBS' in adata_sub.obs['cytokine'].cat.categories:
    adata_sub.obs['cytokine'] = adata_sub.obs['cytokine'].cat.set_categories(
        ['PBS'] + [cat for cat in adata_sub.obs['cytokine'].cat.categories if cat != 'PBS'],
        ordered=True
    )

# 差异基因分析，PBS 作为参考
sc.tl.rank_genes_groups(adata_sub, groupby='cytokine', reference='PBS', method='wilcoxon')

# 提取差异基因信息并存储
result = adata_sub.uns['rank_genes_groups']
groups = result['names'].dtype.names  # 获取 cytokine 组名

# 将所有信息存储到 DataFrame 中，并转换为便于读取的格式
gene_data = {}
for group in groups:
    gene_data[group] = pd.DataFrame({
        'gene': result['names'][group],
        'logfoldchanges': result['logfoldchanges'][group],
        'pvals': result['pvals'][group],
         'pvals_adj': result['pvals_adj'][group]
    })
    
# 保存到主字典中
diff_gene_info[cell_type] = gene_data


### Save and export
## csv format
# 创建保存目录
output_dir = "./results"  # 输出目录
os.makedirs(output_dir, exist_ok=True)

# Replace space with underline
cell_type2 = cell_type.replace(" ", "_")

## Pickle format
output_file = os.path.join(output_dir, f"{cell_type2}_diff_gene_info.pkl")
# 保存文件
with open(output_file, 'wb') as f:
    pickle.dump(diff_gene_info, f)

# 创建一个空的列表，用于存储提取的 'OLAH' 行
combined_olaf_rows = []

# 遍历字典，保存每个 DataFrame 为 CSV 文件
for cell_type, cytokine_dict in diff_gene_info.items():
    for cytokine, df in cytokine_dict.items():
        # 为 DataFrame 添加新列 'cell_cytokine'
        df['cell_cytokine'] = f"{cell_type2}_{cytokine}"
        
        # 创建文件名，使用 {cell_type}_{cytokine}.csv 格式
        file_name = f"{cell_type2}_{cytokine}.csv"
        
        # 将 DataFrame 保存为 CSV 文件
        df.to_csv(os.path.join(output_dir, file_name), index=True)

        # 提取 gene 列为 "OLAH" 的行
        olah_rows = df[df['gene'] == 'OLAH']

        # 将提取的行添加到 combined_olaf_rows 列表中
        combined_olaf_rows.append(olah_rows)

# 合并所有提取的行
final_combined_df = pd.concat(combined_olaf_rows, ignore_index=True)

# 保存合并后的结果为一个 CSV 文件
# 使用 cell_type 动态命名文件名
final_combined_df.to_csv(os.path.join(output_dir, f"{cell_type}_combined_OLAH_rows.csv"), index=False)
```

```bash
#!/bin/bash
#BSUB -P anndata
#BSUB -J cell_type
#BSUB -q superdome
#BSUB -n 30
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=50000]"
#BSUB -oo run_%J.log
#BSUB -eo run_%J.error
#BSUB -W 100:00
#BSUB -B 
#BSUB -N 

module load mamba/1.4.2
mamba activate anndata_env

/home/jqu/.conda/envs/anndata_env/bin/python cell_type.py --cell_type "Plasmablast"
```




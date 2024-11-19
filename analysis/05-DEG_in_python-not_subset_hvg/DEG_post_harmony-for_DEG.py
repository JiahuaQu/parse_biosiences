import os
import sys
#import scipy
#import numpy as np
import pandas as pd
import scanpy as sc
#import seaborn as sns
#import scipy.io as sio
#import scanpy.external as sce
import matplotlib as mpl
import matplotlib.pyplot as plt
#import scipy.sparse as sparse
import warnings
#from numpy.fft import rfft, irfft
#import numpy.linalg as nla


sc.settings.verbosity = 1

# Adjust Scanpy figure defaults
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=400,
    facecolor = 'white', figsize=(6,6), format='png')

analysis_save_path = './'

adata = sc.read_h5ad(os.path.join(analysis_save_path, 'adata_post_harmony.h5ad'))


### DEG
# 用于存储每个 cell_type 的差异基因信息
diff_gene_info = {}

# 遍历每种 cell_type
for cell_type in adata.obs['cell_type'].unique():
    # 子集化数据以包含当前 cell_type 的细胞
    adata_sub = adata[adata.obs['cell_type'] == cell_type]
    
    # 确保 cytokine 列为分类类型
    print(pd.__version__)
    # Check if 'cytokine' is already categorical
    if not pd.api.types.is_categorical_dtype(adata_sub.obs['cytokine']):
        # Convert to categorical only if it's not already categorical
        adata_sub.obs['cytokine'] = adata_sub.obs['cytokine'].astype("category")

    # 将 'PBS' 设置为参考组，并保持其他类别不变
    if 'PBS' in adata_sub.obs['cytokine'].cat.categories:
        # Set categories without using inplace
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

## Pickle format
import pickle
with open('diff_gene_info.pkl', 'wb') as f:
    pickle.dump(diff_gene_info, f)

## csv format
# 创建保存目录
csv_path = './DEG_csv/'
os.makedirs(csv_path, exist_ok=True)

# 遍历字典，保存每个 DataFrame 为 CSV 文件
for cell_type, cytokine_dict in diff_gene_info.items():
    for cytokine, df in cytokine_dict.items():
        # 创建文件名，使用 {cell_type}_{cytokine}.csv 格式
        file_name = f"{cell_type}_{cytokine}.csv"
        # 将 DataFrame 保存为 CSV 文件
        df.to_csv(os.path.join('./DEG_csv/', file_name), index=True)

print("CSV files saved successfully.")

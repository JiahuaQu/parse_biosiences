import os
import sys
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import scipy.io as sio
import scanpy.external as sce
import matplotlib.pyplot as plt

sc.settings.verbosity = 1

# Adjust Scanpy figure defaults
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=400,
    facecolor = 'white', figsize=(6,6), format='png')

import scipy.sparse as sparse
import warnings
from numpy.fft import rfft, irfft
import numpy.linalg as nla

dataset_path = '/research/sharedresources/immunoinformatics/common/jqu/parse_biosciences/datasets/'
analysis_save_path = '/research/sharedresources/immunoinformatics/common/jqu/parse_biosciences/analysis/01-process/'

adata = sc.read_h5ad(os.path.join(analysis_save_path, 'adata_post_pca.h5ad'))

sce.pp.harmony_integrate(adata, 'donor')
adata.write(os.path.join(analysis_save_path, 'adata_post_harmony.h5ad'))

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, use_rep="X_pca_harmony")
sc.tl.umap(adata)
adata.write(os.path.join(analysis_save_path, 'adata_post_harmony_post_umap.h5ad'))

sc.tl.leiden(adata, resolution=1.0, flavor="igraph", n_iterations=2)
adata.write(os.path.join(analysis_save_path, 'adata_post_harmony_post_umap_post_leiden.h5ad'))

UMAP_path = '/research/sharedresources/immunoinformatics/common/jqu/parse_biosciences/analysis/01-process/UMAP/'
os.makedirs(UMAP_path, exist_ok=True)  # 如果目录不存在，则创建

# 绘制和保存 UMAP 图片
sc.pl.umap(adata, color='donor', legend_fontsize=8, save=os.path.join(UMAP_path, "umap_donor.png"))
sc.pl.umap(adata, color='donor', legend_fontsize=8, save=os.path.join(UMAP_path, "umap_donor.pdf"))

sc.pl.umap(adata, color='leiden', legend_fontsize=8, save=os.path.join(UMAP_path, "umap_leiden.png"))
sc.pl.umap(adata, color='leiden', legend_fontsize=8, save=os.path.join(UMAP_path, "umap_leiden.pdf"))

sc.pl.umap(adata, color='leiden', legend_fontsize=8, legend_loc="on data", save=os.path.join(UMAP_path, "umap_leiden_ondata.png"))
sc.pl.umap(adata, color='leiden', legend_fontsize=8, legend_loc="on data", save=os.path.join(UMAP_path, "umap_leiden_ondata.pdf"))

sc.pl.umap(adata, color='cell_type', legend_fontsize=8, save=os.path.join(UMAP_path, "umap_cell_type.png"))
sc.pl.umap(adata, color='cell_type', legend_fontsize=8, save=os.path.join(UMAP_path, "umap_cell_type.pdf"))

sc.pl.umap(adata, color='cell_type', legend_fontsize=8, legend_loc="on data", save=os.path.join(UMAP_path, "umap_cell_type_ondata.png"))
sc.pl.umap(adata, color='cell_type', legend_fontsize=8, legend_loc="on data", save=os.path.join(UMAP_path, "umap_cell_type_ondata.pdf"))

sc.pl.umap(adata, color='cytokine', legend_fontsize=8, save=os.path.join(UMAP_path, "umap_cytokine.png"))
sc.pl.umap(adata, color='cytokine', legend_fontsize=8, save=os.path.join(UMAP_path, "umap_cytokine.pdf"))

sc.pl.umap(adata, color='cytokine', legend_fontsize=8, legend_loc="on data", save=os.path.join(UMAP_path, "umap_cytokine_ondata.png"))
sc.pl.umap(adata, color='cytokine', legend_fontsize=8, legend_loc="on data", save=os.path.join(UMAP_path, "umap_cytokine_ondata.pdf"))

# 在 UMAP 图中显示 OLAH 基因的表达量
sc.pl.umap(adata, color='OLAH', cmap='viridis', legend_fontsize=8, save=os.path.join(UMAP_path, "umap_OLAH.png"))


### DEG
# 用于存储每个 cell_type 的差异基因信息
diff_gene_info = {}

# 遍历每种 cell_type
for cell_type in adata.obs['cell_type'].unique():
    # 子集化数据以包含当前 cell_type 的细胞
    adata_sub = adata[adata.obs['cell_type'] == cell_type]
    
    # 确保 cytokine 列为分类类型
    adata_sub.obs['cytokine'] = adata_sub.obs['cytokine'].astype("category")

    # 将 'PBS' 设置为参考组，并保持其他类别不变
    if 'PBS' in adata_sub.obs['cytokine'].cat.categories:
        adata_sub.obs['cytokine'].cat.set_categories(
            ['PBS'] + [cat for cat in adata_sub.obs['cytokine'].cat.categories if cat != 'PBS'],
            ordered=True,
            inplace=True
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
csv_path = '/research/sharedresources/immunoinformatics/common/jqu/parse_biosciences/analysis/01-process/DEG_csv/'
os.makedirs(csv_path, exist_ok=True)

# 遍历字典，保存每个 DataFrame 为 CSV 文件
for cell_type, cytokine_dict in diff_gene_info.items():
    for cytokine, df in cytokine_dict.items():
        # 创建文件名，使用 {cell_type}_{cytokine}.csv 格式
        file_name = f"{cell_type}_{cytokine}.csv"
        # 将 DataFrame 保存为 CSV 文件
        df.to_csv(os.path.join(csv_path, file_name), index=True)

print("CSV files saved successfully.")


### Draw a heat map or matrix
# Column: cell_type
# Row: cytokine
# Cell: everage of OLAH

#import scanpy as sc
#import pandas as pd
#import numpy as np

# 假设 adata 是您的 AnnData 对象，并且已经运行了 Harmony 整合

# 1. 提取 OLAH 基因的表达数据（使用 Harmony 之后的值）
ol_expression = adata[:, 'OLAH'].X  # 获取 OLAH 基因的表达值
ol_expression_df = pd.DataFrame(ol_expression, index=adata.obs.index, columns=['OLAH'])

# 2. 将表达数据与 obs 中的 cell_type 和 cytokine 合并
ol_expression_df['cell_type'] = adata.obs['cell_type'].values
ol_expression_df['cytokine'] = adata.obs['cytokine'].values

with open('OLAH-1-expression.pkl', 'wb') as f:
    pickle.dump(ol_expression_df, f)

# 3. 计算每个 cell_type 和 cytokine 的平均表达量
heatmap_data = ol_expression_df.groupby(['cytokine', 'cell_type']).mean().unstack()

# 4. 用 NaN 填充缺失的值，确保热图中的所有单元格都有值
heatmap_data = heatmap_data.fillna(0)
with open('OLAH-2-heatmap_data.pkl', 'wb') as f:
    pickle.dump(heatmap_data, f)

# 5. 绘制热图
plt.figure(figsize=(10, 20))  # 设置图形大小
sns.heatmap(heatmap_data, cmap='viridis', annot=False)  # annot=True 显示每个单元格的值
plt.title('Mean Expression of OLAH (Post-Harmony)')
plt.xlabel('Cell Type')
plt.ylabel('Cytokine')

# 6. 保存热图为 PNG 和 PDF 格式
plt.savefig(os.path.join(analysis_save_path, 'mean_expression_OLAH.png'), bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(analysis_save_path, 'mean_expression_OLAH.pdf'), bbox_inches='tight')  # 保存为 PDF


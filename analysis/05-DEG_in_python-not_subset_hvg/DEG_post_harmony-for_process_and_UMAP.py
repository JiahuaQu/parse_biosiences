import os
import sys
#import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
#import scipy.io as sio
import scanpy.external as sce
import matplotlib as mpl
import matplotlib.pyplot as plt
#import scipy.sparse as sparse
import warnings
from numpy.fft import rfft, irfft
import numpy.linalg as nla
import pickle


sc.settings.verbosity = 1

# Adjust Scanpy figure defaults
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=400,
    facecolor = 'white', figsize=(6,6), format='png')


analysis_save_path = './'
adata = sc.read_h5ad(os.path.join(analysis_save_path, 'adata_post_harmony.h5ad'))


### Histogram
# Check if 'OLAH' is in the variable names
if 'OLAH' not in adata.var_names:
    raise ValueError("The gene 'OLAH' was not found in the dataset's variable names.")

# Extract the OLAH expression values, ensuring dense and flattened format
ol_expression = adata[:, 'OLAH'].X

# If ol_expression is a sparse matrix, convert to a dense array
if hasattr(ol_expression, 'toarray'):
    ol_expression = ol_expression.toarray().flatten()  # Convert to dense and flatten
else:
    ol_expression = ol_expression.flatten()  # Just flatten if already dense


### Draw a heat map or matrix
# Column: cell_type
# Row: cytokine
# Cell: everage of OLAH

#import scanpy as sc
#import pandas as pd
#import numpy as np

# 假设 adata 是您的 AnnData 对象，并且已经运行了 Harmony 整合

# 1. 提取 OLAH 基因的表达数据（使用 Harmony 之后的值）
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
plt.savefig(os.path.join('mean_expression_OLAH.png'), bbox_inches='tight', dpi=300)
plt.savefig(os.path.join('mean_expression_OLAH.pdf'), bbox_inches='tight')  # 保存为 PDF

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, use_rep="X_pca_harmony")
sc.tl.umap(adata)
adata.write(os.path.join(analysis_save_path, 'adata_post_harmony_post_umap.h5ad'))

sc.tl.leiden(adata, resolution=1.0, n_iterations=2)
# delete flavor="igraph"

adata.write(os.path.join('adata_post_harmony_post_umap_post_leiden.h5ad'))

### UMAP

# 绘制和保存 UMAP 图片
sc.pl.umap(adata, color='donor', legend_fontsize=8, save=os.path.join("_donor.png"))
sc.pl.umap(adata, color='donor', legend_fontsize=8, save=os.path.join("_donor.pdf"))

sc.pl.umap(adata, color='leiden', legend_fontsize=8, save=os.path.join("_leiden.png"))
sc.pl.umap(adata, color='leiden', legend_fontsize=8, save=os.path.join("_leiden.pdf"))

sc.pl.umap(adata, color='leiden', legend_fontsize=8, legend_loc="on data", save=os.path.join("_leiden_ondata.png"))
sc.pl.umap(adata, color='leiden', legend_fontsize=8, legend_loc="on data", save=os.path.join("_leiden_ondata.pdf"))

sc.pl.umap(adata, color='cell_type', legend_fontsize=8, save=os.path.join("_cell_type.png"))
sc.pl.umap(adata, color='cell_type', legend_fontsize=8, save=os.path.join("_cell_type.pdf"))

sc.pl.umap(adata, color='cell_type', legend_fontsize=8, legend_loc="on data", save=os.path.join("_cell_type_ondata.png"))
sc.pl.umap(adata, color='cell_type', legend_fontsize=8, legend_loc="on data", save=os.path.join("_cell_type_ondata.pdf"))

sc.pl.umap(adata, color='cytokine', legend_fontsize=8, save=os.path.join("_cytokine.png"))
sc.pl.umap(adata, color='cytokine', legend_fontsize=8, save=os.path.join("_cytokine.pdf"))

sc.pl.umap(adata, color='cytokine', legend_fontsize=8, legend_loc="on data", save=os.path.join("_cytokine_ondata.png"))
sc.pl.umap(adata, color='cytokine', legend_fontsize=8, legend_loc="on data", save=os.path.join("_cytokine_ondata.pdf"))

# 在 UMAP 图中显示 OLAH 基因的表达量
sc.pl.umap(adata, color='OLAH', cmap='viridis', legend_fontsize=8, save=os.path.join("_OLAH.png"))


# Create a custom colormap: light yellow for low, dark red for high expression
cmap = mpl.colors.LinearSegmentedColormap.from_list('yellowred', ['yellow', 'darkred'])

# Plot UMAP with custom colormap
sc.pl.umap(adata, color='OLAH', cmap=cmap, legend_fontsize=8, save="_OLAH_yellowred.png")
sc.pl.umap(adata, color='OLAH', cmap=cmap, legend_fontsize=8, save="_OLAH_yellowred.pdf")
import os
import scanpy as sc
import pandas as pd

sc.settings.verbosity = 2

# 加载数据
analysis_save_path = './'
adata = sc.read_h5ad(os.path.join(analysis_save_path, '../05-DEG_in_python-not_subset_hvg/adata_post_harmony.h5ad'))

# 检查是否存在目标基因 OLAH
if 'OLAH' not in adata.var_names:
    raise ValueError("Gene 'OLAH' not found in the dataset.")

# 计算整个数据集中表达 OLAH 值大于 0 的细胞数
total_olah_positive = (adata[:, 'OLAH'].X > 0).sum()

# 计算 OLAH 是否表达的布尔值列
olah_positive_by_cell = (adata[:, 'OLAH'].X > 0).toarray().flatten()

# 添加布尔列到 obs
adata.obs['OLAH_positive'] = olah_positive_by_cell

# 计算每种 cell_type 的总细胞数
cell_type_counts = adata.obs['cell_type'].value_counts().reset_index()
cell_type_counts.columns = ['cell_type', 'total_cells']

# 计算每种 cell_type 中表达 OLAH 的细胞数
olah_positive_by_cell_type = (
    adata.obs.groupby('cell_type')['OLAH_positive'].sum().reset_index()
)
olah_positive_by_cell_type.columns = ['cell_type', 'olah_positive_cells']

# 合并两个数据框
summary_df = pd.merge(cell_type_counts, olah_positive_by_cell_type, on='cell_type', how='left')

# 计算总数行
total_row = pd.DataFrame([{
    'cell_type': 'Total',
    'total_cells': adata.n_obs,
    'olah_positive_cells': total_olah_positive
}])

# 添加总数行到结果数据框
final_summary = pd.concat([summary_df, total_row], ignore_index=True)

# 导出为 CSV 文件
output_file = 'num.csv'
final_summary.to_csv(output_file, index=False)
print(f"Summary saved to {output_file}")

# 打印结果用于检查
print("Summary:")
print(final_summary)

import os
import pandas as pd
import scanpy as sc

# Adjust Scanpy figure defaults
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=400,
    facecolor = 'white', figsize=(6,6), format='png')


analysis_save_path = './'
adata = sc.read_h5ad(os.path.join(analysis_save_path, '../05-DEG_in_python-not_subset_hvg/adata_post_harmony.h5ad'))

# 获取 cell_type 的类别及其计数
cell_type_counts = adata.obs['cell_type'].value_counts()

# 将结果转换为 DataFrame
cell_type_df = pd.DataFrame({
    'cell_type': cell_type_counts.index,  # 类别名称
    'count': cell_type_counts.values     # 对应的计数
})

# 导出为 CSV 文件
output_file = "./cell_type_counts.csv"  # 指定输出文件路径
cell_type_df.to_csv(output_file, index=False)

import matplotlib.pyplot as plt

# 设置图形大小
plt.figure(figsize=(10, 6))

# 绘制柱状图
plt.bar(cell_type_df['cell_type'], cell_type_df['count'], color='skyblue', edgecolor='black')

# 添加标题和标签
plt.title('Cell Type Distribution', fontsize=16)
plt.xlabel('Cell Type', fontsize=14)
plt.ylabel('Count', fontsize=14)

# 旋转 x 轴标签以避免重叠
plt.xticks(rotation=45, ha='right', fontsize=12)

# 添加网格线
plt.grid(axis='y', linestyle='--', alpha=0.7)

# 保存图像（可选）
plt.savefig('cell_type_distribution.png', dpi=300, bbox_inches='tight')
plt.savefig('cell_type_distribution.pdf', bbox_inches='tight')

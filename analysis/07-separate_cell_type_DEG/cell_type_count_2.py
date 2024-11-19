import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
import scanpy as sc

# 路径设置
analysis_save_path = './'
adata = sc.read_h5ad(os.path.join(analysis_save_path, '../05-DEG_in_python-not_subset_hvg/adata_post_harmony.h5ad'))

# 确保 adata.obs 中包含 'cell_type' 和 'cytokine' 列
# 统计组合类别的数量
combination_counts = adata.obs.groupby(['cell_type', 'cytokine']).size().reset_index(name='count')

# 转换为 1e4 为单位
combination_counts['scaled_count'] = combination_counts['count'] / 1e4

# 创建一个pivot表，便于热图绘制
pivot_table = combination_counts.pivot(index='cytokine', columns='cell_type', values='scaled_count').fillna(0)

# 导出统计结果为 CSV 文件
combination_counts.to_csv(os.path.join(analysis_save_path, 'combination_counts_2.csv'), index=False)
pivot_table.to_csv(os.path.join(analysis_save_path, 'pivot_table.csv'))

# 绘制气泡热图
plt.figure(figsize=(10, 35))

# 获取坐标和大小
x_labels = pivot_table.columns  # 横坐标（cell_type）
y_labels = pivot_table.index    # 纵坐标（cytokine）
sizes = pivot_table.values      # 气泡大小

# 绘制气泡
for i, y in enumerate(y_labels):
    for j, x in enumerate(x_labels):
        size = sizes[i, j]
        plt.scatter(j, i, s=size * 100, color='skyblue', edgecolor='black')  # s 控制气泡大小
        if size > 0:
            plt.text(j, i, f"{size:.1f}", ha='center', va='center', fontsize=10)  # 在气泡中标注数值

# 设置轴标签
plt.xticks(ticks=np.arange(len(x_labels)), labels=x_labels, rotation=45, ha='right', fontsize=12)
plt.yticks(ticks=np.arange(len(y_labels)), labels=y_labels, fontsize=12)
plt.xlabel('Cell Type', fontsize=14)
plt.ylabel('Cytokine', fontsize=14)
plt.title('Bubble Heatmap of Cell Type and Cytokine Combinations (Unit: 1e4)', fontsize=16)

# 添加图例（放置在图外右侧）
bubble_sizes = [1, 10, 30, 50]  # 气泡大小对应的示例数值
legend_handles = [
    plt.scatter([], [], s=size * 100, color='skyblue', edgecolor='black', label=f"{size} x 1e4")
    for size in bubble_sizes
]
plt.legend(
    handles=legend_handles, title='Counts (Unit: 1e4)', title_fontsize=12, fontsize=10, loc='center left',
    bbox_to_anchor=(1.05, 0.5), frameon=True
)

# 添加网格线
plt.grid(axis='both', linestyle='--', alpha=0.5)

# 保存图像
output_png = os.path.join(analysis_save_path, 'cell_type_cytokine_distribution.png')
output_pdf = os.path.join(analysis_save_path, 'cell_type_cytokine_distribution.pdf')
output_pickle = os.path.join(analysis_save_path, 'cell_type_cytokine_distribution.pkl')

plt.savefig(output_png, dpi=300, bbox_inches='tight')
plt.savefig(output_pdf, bbox_inches='tight')

# 保存 Matplotlib 图形对象为 Pickle 文件
fig = plt.gcf()  # 获取当前图形对象
#output_pickle = os.path.join(analysis_save_path, 'cell_type_cytokine_distribution.pkl')

with open(output_pickle, 'wb') as f:
    pickle.dump(fig, f)

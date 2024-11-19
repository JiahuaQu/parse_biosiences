import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
import scanpy as sc

# 设置文件路径
pivot_table_csv_path = './pivot_table.csv'

# 读取 pivot_table.csv
pivot_table = pd.read_csv(pivot_table_csv_path, index_col=0)

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

# 保存 Matplotlib 图形对象为 Pickle 文件
fig = plt.gcf()  # 获取当前图形对象
output_pickle = 'cell_type_cytokine_distribution.pkl'

with open(output_pickle, 'wb') as f:
    pickle.dump(fig, f)

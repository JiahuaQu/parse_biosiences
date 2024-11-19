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
parser.add_argument('--cell_types', type=str, help='The cell type to analyze')

# 解析参数
args = parser.parse_args()

# 获取参数值
cell_types = args.cell_types

# 打印参数进行调试
print(f"Cell type: {cell_types}")

sc.settings.verbosity = 1

analysis_save_path = './'
adata = sc.read_h5ad(os.path.join(analysis_save_path, 'adata_post_harmony.h5ad'))

### DEG
# 用于存储每个 cell_type 的差异基因信息
diff_gene_info = {}

# 遍历每种 cell_type
for cell_type in cell_types:
    # 子集化数据以包含当前 cell_type 的细胞
    adata_sub = adata[adata.obs['cell_type'] == cell_type]
    
	# 确保 cytokine 列为分类类型
    if 'cytokine' not in adata_sub.obs:
        print(f"'cytokine' column not found for cell type {cell_type}. Skipping...")
        continue

	# Check if 'cytokine' is already categorical
    if not isinstance(adata_sub.obs['cytokine'].dtype, pd.CategoricalDtype):
        # Convert to categorical only if it's not already categorical
        adata_sub.obs['cytokine'] = adata_sub.obs['cytokine'].astype("category")

	# 检查 'PBS' 是否在 cytokine 的类别中
    if 'PBS' not in adata_sub.obs['cytokine'].unique():
        print(f"'PBS' not found in 'cytokine' for cell type {cell_type}. Skipping...")
        continue

    # 确保 adata_sub 是一个副本
    adata_sub = adata_sub.copy()
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
## csv format
# 创建保存目录
output_dir = "./results"  # 输出目录
os.makedirs(output_dir, exist_ok=True)

# Replace space with underline
cell_types2 = cell_types.replace(" ", "_")

## Pickle format
output_file = os.path.join(output_dir, f"{cell_types2}_diff_gene_info.pkl")
# 保存文件
with open(output_file, 'wb') as f:
    pickle.dump(diff_gene_info, f)

# 遍历字典，保存每个 DataFrame 为 CSV 文件
for cell_type, cytokine_dict in diff_gene_info.items():
    for cytokine, df in cytokine_dict.items():
        # 创建文件名，使用 {cell_type}_{cytokine}.csv 格式
        file_name = f"{cell_type2}_{cytokine}.csv"
        # 将 DataFrame 保存为 CSV 文件
        df.to_csv(os.path.join(output_dir, file_name), index=True)


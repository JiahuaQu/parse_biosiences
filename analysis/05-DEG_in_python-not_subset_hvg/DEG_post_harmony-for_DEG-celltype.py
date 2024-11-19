import os
import scanpy as sc
import pickle

# 设置路径
analysis_save_path = './'
adata_file = os.path.join(analysis_save_path, 'adata_post_harmony.h5ad')
pickle_file = os.path.join(analysis_save_path, 'cell_types.pkl')
txt_file = os.path.join(analysis_save_path, 'cell_types.txt')

# 读取 AnnData 文件
adata = sc.read_h5ad(adata_file)

# 提取唯一的 cell_type 值
cell_types = adata.obs['cell_type'].unique()

# 保存为 pickle 文件
with open(pickle_file, 'wb') as f:
    pickle.dump(cell_types, f)

# 保存为 txt 文件
with open(txt_file, 'w') as f:
    for cell_type in cell_types:
        f.write(f"{cell_type}\n")

# 验证写入是否成功
print(f"Cell types saved to:\n- Pickle: {pickle_file}\n- TXT: {txt_file}")

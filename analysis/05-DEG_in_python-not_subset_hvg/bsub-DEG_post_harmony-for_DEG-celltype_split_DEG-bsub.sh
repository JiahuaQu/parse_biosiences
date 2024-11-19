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
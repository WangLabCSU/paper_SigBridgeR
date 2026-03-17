#!/usr/bin/bash

# ==========================================
# 脚本名称: setup_folders.sh
# 功能: 批量创建实验目录结构及 R 脚本占位文件
# ==========================================

# 1. 定义顶层目录列表
# 注意：你列出了 7 个目录，脚本将处理所有这 7 个目录
top_dirs=("DEGAS" "PIPET" "Scissor" "scPP" "LP_SGL" "scAB" "scPAS")

# 2. 定义子目录列表
sub_dirs=("survival" "binary")

# 3. 定义文件后缀列表
# 注意：这里严格按照你的要求保留了 "10l"
file_suffixes=("1k" "5k" "10k" "50k" "100k")

echo "正在开始创建目录结构和文件..."

# 开始循环处理
for top in "${top_dirs[@]}"; do
    for sub in "${sub_dirs[@]}"; do
        
        # 构建并创建路径：例如 DEGAS/survival
        target_path="./${top}/${sub}"
        mkdir -p "$target_path"
        
        # 在该路径下创建 5 个 R 脚本
        for suffix in "${file_suffixes[@]}"; do
            file_name="seurat_${suffix}.R"
            
            # 使用 touch 创建空文件
            touch "${target_path}/${file_name}"
        done
        
    done
    echo "已处理目录: ${top}"
done

echo "--------------------------------------"
echo "所有任务执行完毕！"
#!/bin/bash

# ==========================================
# 脚本名称: copy_and_replace.sh
# 功能: 复制模板代码并批量替换文本内容
# ==========================================

# 1. 定义顶层目录列表
top_dirs=("DEGAS" "PIPET" "Scissor" "scPP" "LP_SGL" "scAB" "scPAS")

# 2. 定义子目录列表
sub_dirs=("survival" "binary")

# 3. 定义目标尺寸列表 (源文件是 1k，所以这里列出要生成的目标)
# 注意：这里使用了 10k 而不是 10l
target_sizes=("5k" "10k" "50k" "100k")

echo "开始执行复制与替换任务..."

# 开始遍历目录
for top in "${top_dirs[@]}"; do
    for sub in "${sub_dirs[@]}"; do
        
        # 定义当前的工作路径
        work_path="./${top}/${sub}"
        
        # 定义源文件
        source_file="${work_path}/seurat_1k.R"
        
        # 检查源文件是否存在，如果不存在则跳过，防止报错
        if [ ! -f "$source_file" ]; then
            echo "[警告] 找不到源文件: $source_file，跳过该目录。"
            continue
        fi

        # 遍历需要生成的尺寸
        for size in "${target_sizes[@]}"; do
            target_file="${work_path}/seurat_${size}.R"
            
            # 使用 sed 命令读取源文件，替换字符，并写入新文件
            # s/1k/${size}/g 意思是：substitute(替换) 所有的 1k 为 对应的size，g代表全局替换
            sed "s/1k/${size}/g" "$source_file" > "$target_file"
            
            # (可选) 打印详细信息，如果觉得太吵可以注释掉下面这行
            echo "已生成: $target_file (内容中 1k 已替换为 $size)"
        done
        
    done
    echo "目录处理完成: ${top}"
done

echo "--------------------------------------"
echo "所有文件的复制与文本替换已完成！"  
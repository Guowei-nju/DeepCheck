#!/bin/bash
B_SCRIPT="/data/yan_code/deepcheck/data/random_protein_fragmentation.sh"

if [ ! -f "$B_SCRIPT" ]; then
  echo "B脚本不存在或路径不正确。请检查B_SCRIPT变量的值。"
  exit 1
fi

DATA_DIR="/data/yan_code/deepcheck/data/results_test/"

# 使用sem命令限制并发进程数为128
sem -j128 --id myjob

# 遍历目录下的文件并处理
for file in "$DATA_DIR"/*; do
  if [ -f "$file" ]; then
    sem --id myjob "$B_SCRIPT" "$file"
  fi
done

# 等待所有任务完成
sem --wait --id myjob
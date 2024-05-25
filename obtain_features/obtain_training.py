import os
import multiprocessing
import subprocess

def process_file(file_path, b_script):
    if os.path.isfile(file_path):
        subprocess.call(["/bin/bash", b_script, file_path])

def main():
    b_script = "./random_protein_fragmentation.sh"
    data_dir = "./training_proteins/"
    
    if not os.path.isfile(b_script):
        print("B脚本不存在或路径不正确。请检查b_script变量的值。")
        exit(1)

    pool = multiprocessing.Pool(processes=128)  # 使用所有可用的CPU核心
    jobs = []

    for root, dirs, files in os.walk(data_dir):
        for file in files:
            file_path = os.path.join(root, file)
            job = pool.apply_async(process_file, args=(file_path, b_script))
            jobs.append(job)

    # 等待所有任务完成
    for job in jobs:
        job.get()

    pool.close()
    pool.join()

if __name__ == "__main__":
    main()

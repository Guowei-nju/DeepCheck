import os
import subprocess

# 指定输入文件夹
""" input_folder = 'refseq'
genes_folder= 'traing_data/training_genes'
proteins_folder= 'traing_data/training_proteins' """

input_folder = './data/'
genes_folder= './Patescibacteria_genes'
proteins_folder= './Patescibacteria_proteins'

# 遍历输入文件夹中的所有.fna文件
for filename in os.listdir("no_rundent_refseq89"):
    if filename.endswith('.fna'):
        # 构造输出文件名
        base_name = os.path.splitext(filename)[0]
        genome_input = f'{base_name}.fna'
        genes_output = f'{base_name}.genes'
        proteins_output = f'{base_name}.faa'
        
        # 构造prodigal命令
        prodigal_command = f'prodigal -i {os.path.join(input_folder, genome_input)} -o {os.path.join(genes_folder, genes_output)} -a {os.path.join(proteins_folder, proteins_output)}'

        # 执行prodigal命令
        subprocess.run(prodigal_command, shell=True)

        print(f'Processed {filename} and generated {genes_output} and {proteins_output}')   
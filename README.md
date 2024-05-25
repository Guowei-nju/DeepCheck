![DeepCheck Logo](https://github.com/Guowei-nju/DeepCheck/blob/main/figures/log.png)<!-- 替换为实际的Logo URL -->


**Multi-task learning aids in assessing microbial genome quality**

![DeepCheck Banner](https://github.com/Guowei-nju/DeepCheck/blob/main/figures/frame.png) <!-- 替换为实际的Banner URL -->

Metagenomic analyses afford biologists the opportunity to delve into the microbial world, thereby enhancing our comprehension of the contributions of microbes to a myriad of ecological and biological processes. Evaluating the quality of metagenome-assembled genomes (MAGs) represents a pivotal step within metagenomic analysis and serves as a fundamental prerequisite for the derivation of accurate biological insights. However, the prevailing machine learning-based methods for quality assessment train two distinct models to predict the completeness and contamination of MAGs, respectively. This approach overlooks the intrinsic relationship between these two tasks, resulting in diminished generalization capability of the models.

Here, we propose a multi-tasking deep learning framework for simultaneous prediction of completeness and contamination of MAGs, named **DeepCheck**. In various experimental scenarios, DeepCheck outperforms other existing tools in accuracy with considerable prediction speed. Even for new lineages, DeepCheck still achieves high-quality prediction accuracy. Furthermore, we also employ interpretable machine learning techniques to identify specific genes and pathways upon which the model bases its predictions. This facilitates the independent investigation and assessment of these biological elements for further insight.

# Run without installing

For simplicity, you can just download CheckM2 from GitHub and run it directly without installing.

Retrieve the files:
```
git clone --recursive https://github.com/Guowei-nju/DeepCheck.git && cd checkm2

```
Create an appropriate Conda environment with prerequisites using the `DeepCheck.yml` file:

```
conda env create -n DeepCheck -f DeepCheck.yml
conda activate DeepCheck
```


# Training your data

If you want to train your own data, you'll need to use the [`obtain_training.py`](obtain_features/obtain_training.py) in the `contain_features` folder to simulate genome integrity and contamination, and then use [`generate_features.py`](generate_features.py) to generate the features that DeepCheck needs as input. After generating the features, you need to read the features and simulate the labels (completeness and contamination) using the `process_data` function in [`process_data.py`](process_data.py) as:
- `train_data_path`
- `train_comp_path`
- `train_y_cont`

Finally, call:
```bash
python train.py --train_x train_data_path --train_y_comp train_comp_path --train_y_cont train_y_cont
```
# Direct prediction using our model

We also provide the associated notebooks for use in making predictions directly:

[Evaluating Macro Genome Assemblies](https://github.com/Guowei-nju/DeepCheck/blob/main/evaluating%20macro%20genome%20assemblies!.ipynb)

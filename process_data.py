import pickle
from os import listdir
import pandas as pd
from scipy.sparse import csr_matrix, vstack
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix, vstack, hstack
from scipy import sparse
import scipy
import sys
import warnings

import os
#读取pkl文件
import pickle
#读取pkl文件
import pandas as pd
#读取pkl文件
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

def process_data(cm2_folders, cont_threshold=100000):
    # Initialize lists to hold vectors and load scaler
    comp_scaled_vectors, comp_scaled_labels_comp, comp_scaled_labels_cont = [], [], []
    comp_raw_vectors, comp_raw_labels = [], []
    cont_raw_vectors, cont_raw_labels = [], []
    import pickle
    scalerfile = "scaler.sav"
    scaler = pickle.load(open(scalerfile, 'rb'))
    count=0
    for folder in cm2_folders:
        print ("this is the "+str(count)+"th file")
        print (folder)
        exclude_from_cont = False
        pickled_location = f'data/traing_data/output/{folder}'
        pickles = listdir(pickled_location)
        pickles = sorted([f for f in pickles if f.endswith('.pkl')])
        interim_pickled_list = []
        for pickle in pickles:
            pickled_vector = pd.read_pickle(f'checkm2_output/{folder}/{pickle}') #读取checkm2输出的pkl文件作为特征
            pickled_vector['Completeness'] = pickled_vector['Name'].apply(lambda x: x.split('_pc_complete')[0].split('_')[-1]).astype(float)
            pickled_vector['Contamination'] = pickled_vector['Name'].apply(lambda x: x.split('pc_contaminated')[0].split('_pc_completeXXX')[-1]).astype(float)
            interim_pickled_list.append(pickled_vector)

        full_pickles = pd.concat(interim_pickled_list)
        if full_pickles[(full_pickles['Completeness'] == 100) & (full_pickles['Completeness'] == 100)]['AALength'].values[0] < cont_threshold:
            exclude_from_cont = True

        # Make raw comp vectors
        del full_pickles['Name']
        comp_raw_vectors.append(csr_matrix(full_pickles.iloc[:, :-2].values))
        comp_raw_labels.append(np.array(full_pickles['Completeness'].values))

        # Make raw cont vectors
        if not exclude_from_cont:
            cont_raw_vectors.append(csr_matrix(full_pickles.iloc[:, :-2].values))
            cont_raw_labels.append(np.array(full_pickles['Contamination'].values))

        # Make scaled vectors
        comp_scaled_vectors.append(csr_matrix(scaler.transform(full_pickles.iloc[:, :-2].values))[:, :20021])
        comp_scaled_labels_comp.append(np.array(full_pickles['Completeness'].values))
        comp_scaled_labels_cont.append(np.array(full_pickles['Contamination'].values))
        count=count+1
    # Concatenate all vectors
    comp_scaled_vectors = vstack(comp_scaled_vectors)
    comp_scaled_labels_comp = np.concatenate(comp_scaled_labels_comp)
    comp_scaled_labels_cont = np.concatenate(comp_scaled_labels_cont)
    comp_raw_vectors = vstack(comp_raw_vectors)
    comp_raw_labels = np.concatenate(comp_raw_labels)
    cont_raw_vectors = vstack(cont_raw_vectors)
    cont_raw_labels = np.concatenate(cont_raw_labels)

    return (
        comp_scaled_vectors, comp_scaled_labels_comp, comp_scaled_labels_cont,
        comp_raw_vectors, comp_raw_labels,
        cont_raw_vectors, cont_raw_labels
    )
def read_feature(fasta_names, cont_threshold=100000):
    # Initialize lists to hold vectors and load scaler
    comp_scaled_vectors, cont_raw_vectors = [], []
    comp_raw_vectors = []
    import pickle
    from os import listdir
    import pandas as pd
    import numpy as np
    from scipy.sparse import csr_matrix, vstack

    scalerfile = "scaler.sav"
    scaler = pickle.load(open(scalerfile, 'rb'))
    count = 0
    
    for folder in fasta_names:
        print(folder)
        exclude_from_cont = False
        pickled_location = f'/data/yan_code/deepcheck/data/traing_data/output/{folder}/'
        pickles = listdir(pickled_location)
        pickles = sorted([f for f in pickles if f.endswith('.pkl')])
        interim_pickled_list = []
        
        for pickle in pickles:
            pickled_vector = pd.read_pickle(f'/data/yan_code/deepcheck/data/traing_data/output/{folder}/{pickle}')

            interim_pickled_list.append(pickled_vector)
        
        full_pickles = pd.concat(interim_pickled_list)
        
        # Make raw comp vectors
        del full_pickles['Name']
        comp_raw_vectors.append(csr_matrix(full_pickles.iloc[:, :-2].values))
        
        # Make raw cont vectors
        if not exclude_from_cont:
            cont_raw_vectors.append(csr_matrix(full_pickles.iloc[:, :-2].values))
        
        # Make scaled vectors
        comp_scaled_vectors.append(csr_matrix(scaler.transform(full_pickles.iloc[:, :-2].values))[:, :20021])
        count += 1
    
    # Concatenate all vectors
    comp_scaled_vectors = vstack(comp_scaled_vectors)
    comp_raw_vectors = vstack(comp_raw_vectors)
    cont_raw_vectors = vstack(cont_raw_vectors)
    
    return comp_scaled_vectors
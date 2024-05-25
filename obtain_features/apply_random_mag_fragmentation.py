import os
import sys
import numpy as np
import pandas as pd
from random import shuffle


# import file
file_list = []
#this is a 1-contig complete genome
input_file = str(sys.argv[1])
#this is the complete genome with contig fragmentation applied
output_file = str(sys.argv[2])
#this is the "MAG_CONTIG_LENGTH_DISTRIBUTIONS" folder
fragment_dir = str(sys.argv[3])
for f in [f for f in os.listdir(fragment_dir) if f.endswith(".tsv")]:
    file_list.append(f)
file_list.sort()

#start by setting up our random frag data
frags_to_choose_from = len(file_list)
random_pick = int(np.random.randint(1,frags_to_choose_from, size=1))
chosen_file = file_list[random_pick]

#open our chosen file into pandas
fragment_size_list = pd.read_csv("{}{}".format(fragment_dir, chosen_file), sep='\t', index_col=None, names=['ContigSize'])
total_size = fragment_size_list['ContigSize'].sum()
#calculate the proportional size of every contig relative to genome size
fragment_size_list['RelativeSize'] = fragment_size_list['ContigSize'].apply(lambda x: x/total_size)

#lets now read in our file:
with open(input_file, 'r') as file:
    fasta_sequence = file.read().replace('\n', '')
    
input_genome_len = len(fasta_sequence)
#now let's calculate the absolute length (in nucleotides) of each fragment that needs to be generated
fragment_size_list['InputFileFragmentLength'] = fragment_size_list['RelativeSize'].apply(lambda x: int(round(x*input_genome_len,0)))
input_frag_size_list = fragment_size_list['InputFileFragmentLength'].tolist()
index=0
running_char_tally = 0
contig_list = []
shuffle(input_frag_size_list)
for entry in input_frag_size_list:
    if index is not len(input_frag_size_list):
        contig_list.append('>{}_{}'.format(index, input_file[:-4].split('/')[-1]))
        contig_list.append(fasta_sequence[running_char_tally:(running_char_tally + entry)])
        index = index + 1
        running_char_tally = running_char_tally + entry
    else:
        contig_list.append('>{}'.format(index))
        contig_list.append(fasta_sequence[running_char_tally:]) # if it's the last contig, continue until the end of fasta sequence
    
#write contigs to newly fragmented file: 
with open(output_file, "w") as fout:
    print(*contig_list, sep="\n", file=fout)

#!/usr/bin/env python

import json, os, gc, re
from Bio import SeqIO
import pandas as pd
import numpy as np
from tqdm import tqdm
from glob import glob
from argparse import ArgumentParser
from collections import OrderedDict

parser = ArgumentParser(description="Merge MethylDackel BedGraph files from from all the samples into methylated and coverage files.")
parser.add_argument('json_file', help="Json containing metadata")
parser.add_argument('-g',  '--genome_dir', help="Directory with the genome", required=True)
parser.add_argument('--prefix', help="Base path for data files", default='.')
args= parser.parse_args()

def replicate_iter(jsonfile, base_path ) : 
    data = json.load(open(jsonfile))
    experiments = data['experiments']

    for experiment in experiments :
        experiment_path = os.path.join(base_path, experiment['id'])
        samples = experiment['samples']

        for sample in samples : 
            sample_path = os.path.join(experiment_path, sample['id'])
            
            for r, replicate in enumerate(sample['accession']) :
                replicate_path = os.path.join(sample_path, str(r))

                assert os.path.isdir(replicate_path)

                yield (experiment['id'], sample['id'], r) , replicate_path


basename = json.load(open(args.json_file))['genome']['basename']
genome_path = os.path.join(args.genome_dir, f'{basename}.fna')

print (f'Found genome: {genome_path}')

genome = SeqIO.to_dict(SeqIO.parse(genome_path, 'fasta'))

#c_list = list()
c_index = OrderedDict()

for seqid in tqdm(genome) :
    #c_index[seqid] = dict()
    for match in re.finditer('[CcGgn]', str(genome[seqid].seq)) : 
        i = match.start()
        c_index[(seqid,i)]=len(c_index)
        #c_list.append((seqid,i))
        
        
        # char = match.group()

        # c_masked.append(char.islower())
        # if char in ('C','c') :
        #     c_strand.append(True)
        #     c_context.append(genome[seqid][i:i+3])
        # else : 
        #     c_strand.append(False)
        #     c_context.append(genome[seqid][i-2:i+1].reverse_complement())

#masked = list(map( lambda pos: genome[pos[0]][pos[1]].islower(), c_list))

index = pd.MultiIndex.from_tuples(c_index.keys(), names=('seqid', 'pos'))

replicate_paths = list(replicate_iter(args.json_file, args.prefix))
nsamples = len(replicate_paths)

print(f'{len(c_index)} cytosines in the genome.')
print(f'{nsamples} samples found.')

methylated = np.zeros((len(c_index),nsamples), dtype='uint16')
covered = np.zeros((len(c_index),nsamples), dtype='uint16')

sample_list = list()

for j, (label, replicate_path) in enumerate(tqdm(replicate_paths)) :

    files = glob(os.path.join(replicate_path, 'results','MethylDackel', '*.bedGraph'))

    sample_list.append(label)

    assert len(files) %3 == 0

    for file in files :
        #print(file)
        # with open(file) as f :
        #     assert f.readline().startswith('track type')
        #     for line in f : 
        #         seqid, pos, x , y , m, u = line.strip().split()
        #         i = c_index[seqid][int(pos)]
        #         methylated[i,j]+= int(m)
        #         covered[i,j]   += int(m)+int(u)
        df = pd.read_table(file, 
                           skiprows=1,
                           names=('seqid','pos','y','x','m','u'), 
                           usecols=('seqid','pos','m','u'), 
                           dtype={'m':'uint16', 'u':'uint16'})
        df['c'] = df.m+df.pop('u')
        for seqid, pos, m, c in df.itertuples(index=False) :
            #i = c_index[seqid][pos]
            i = c_index[(seqid,pos)]
            methylated[i,j] += m
            covered[i,j] +=c


columns = pd.MultiIndex.from_tuples(sample_list, names=('study','sample','replicate'))

methylated = pd.DataFrame(methylated, index=index, columns=columns)
covered    = pd.DataFrame(covered, index = index, columns=columns)

# Save master DataFrames
methylated.to_pickle('methylated.pkl')
covered.to_pickle('covered.pkl')


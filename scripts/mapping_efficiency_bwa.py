#!/usr/bin/env python

import json, os, gc, re
from Bio import SeqIO
import pandas as pd
import numpy as np
from tqdm import tqdm
from glob import glob
from argparse import ArgumentParser
from collections import OrderedDict

parser = ArgumentParser(description="Parses bwa-mem mapping efficiency.")
parser.add_argument('json_file', help="Json containing metadata")
parser.add_argument('--base_path', help="Base path for data files", default='.')
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


replicate_paths = list(replicate_iter(args.json_file, args.base_path))
nsamples = len(replicate_paths)

mapping_eff = dict()

for label, replicate_path in replicate_paths :

    files = glob(os.path.join(replicate_path, 'results','bwa-mem_alignments', 'logs','*stats_report.txt'))

    mefficiency = list()

    for file in files :
        with open(file) as f:
            for line in f :
                if line.startswith('SN\traw total sequences:'):
                    n_sequences = int(line.strip().split()[-1])
                elif line.startswith('SN\treads mapped:') :
                    n_mapped = int(line.strip().split()[-1])
                    break
            mefficiency.append(n_mapped/n_sequences)

    mapping_eff[label] = np.mean(mefficiency)

pd.Series(mapping_eff).sort_index().to_csv('mapping_efficiency.csv', sep=' ', header=None)





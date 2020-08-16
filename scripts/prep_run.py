#!/usr/bin/env python

from sys import argv
import json
import os 
from uuid import uuid4
from glob import glob
from argparse import ArgumentParser
from shutil import move

parser = ArgumentParser(description="Helper script to create run directory.")
parser.add_argument(dest='json_file', help="master json file")
parser.add_argument('-base_dir', default='.')
parser.add_argument('-fastq_dir', required=True)
args = parser.parse_args()

data = json.load(open(args.json_file))

experiments = data['experiments']

for experiment in experiments :
    if 'id' not in experiment : 
        experiment['id'] =  uuid4().hex[:4]
    
    experiment_dir = os.path.join(args.base_dir, experiment['id'])

    for sample in experiment['samples'] : 

        if 'id' not in sample : 
            sample['id'] = uuid4().hex[:4]

        sample_dir = os.path.join(experiment_dir, sample['id'])
        

        for replicate, replicate_accessions in enumerate (sample['accession']) :
            replicate_dir = os.path.join(sample_dir, str(replicate))
            os.makedirs(replicate_dir, exist_ok=True)
            for accession in replicate_accessions : 
                print (accession) 
                fastq_files = glob(os.path.join(os.path.abspath(args.fastq_dir), accession + '*.fastq.gz'))
                assert len(fastq_files) > 0 
                for fn in  fastq_files:
                    os.link(fn, os.path.join(replicate_dir, os.path.basename(fn)) )
                
        with open(os.path.join(sample_dir, 'sample.json'), 'w') as f :
            json.dump(sample, f)

    with open(os.path.join(experiment_dir, 'experiment.json'), 'w') as f :
        json.dump(experiment, f)

move(args.json_file, args.json_file+'.backup')
with open(args.json_file, 'w') as f:
    json.dump(data, f)


#!/usr/bin/env python
# coding: utf-8

# Boilerplate code for notebooks


import numpy as np
import pandas as pd
import json, pickle, gc
from sys import argv 
from os import path, makedirs

from utils import *

prefix = argv[1]

# Derive filenames from prefix

meta_file = path.join(prefix, 'meta.json')

metadata = json.load(open(meta_file, 'r'))

genome_file = f"{metadata['genome']['basename']}.fna"
annotation_file = f"{metadata['genome']['basename']}.gff"

organism = metadata['organism']
organism_prefix = path.basename(prefix)


print(f'Running boilerplate for {organism} based off {prefix}')


# Assembly
genome = load_genome(path.join(prefix, genome_file ))

print(
    f"Total {sum([len(sequence) for sequence in genome.values()])} bps in {len(genome)} contigs."
)

# Annotation
db = load_annotation(path.join(prefix, annotation_file))

print(f"Total {db.count_features_of_type('gene')} genes.")

# Parse accession numbers for genes 

gene_accession = dict()

for gene in db.all_features(featuretype='gene') :
    for label in gene.attributes['Dbxref'] :
        if label.startswith('GeneID') :
            gene_accession[gene.id] = label.split(':')[1] 


# Choosing one mrna per gene

print(f"Choosing a representative mrna per protein coding gene..")

# Note that protein coding genes should be all genes
# if we are using canon-gff parsed version of the annotation
protein_coding_genes = list(
    filter(
        lambda gene: gene.attributes["gene_biotype"][0] == "protein_coding",
        db.all_features(featuretype="gene"),
    )
)

representative_mrnas = list()

for gene in protein_coding_genes:
    mrnas = list(db.children(gene, featuretype="mRNA"))

    if len(mrnas) == 1:
        representative_mrnas.append(mrnas[0])
    else:
        # for multiple transcripts, pick longest translated product
        len_longest_mrna = 0
        for mrna in mrnas:
            mrna_len = db.children_bp(mrna, child_featuretype="CDS")
            if mrna_len > len_longest_mrna:
                longest_mrna = mrna
                len_longest_mrna = mrna_len
        representative_mrnas.append(longest_mrna)

print(f'{len(representative_mrnas)} representative mrnas.')

print("Loading coverage matrix..")
covered = pd.read_pickle(path.join(prefix, 'covered_filtered.pkl'))

print("Loading methylation matrix..")
methylated = pd.read_pickle(path.join(prefix, 'methylated_significant.pkl'))

assert methylated.columns.equals(covered.columns)

print('Loading index..')
index = pd.read_pickle(path.join(prefix, 'index.pkl'))
index.set_index(['seqid', 'pos'], inplace=True)
index.strand = index.strand.astype(bool)
assert covered.index.equals(index.index)

samples = parse_experiments(meta_file)
samples = samples.reindex(methylated.columns)
assert samples.index.equals(methylated.columns)
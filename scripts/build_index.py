#!/usr/bin/env python

import pandas as pd
import numpy as np
from tqdm import tqdm
from utils import *

from argparse import ArgumentParser
from os import path

parser = ArgumentParser(description="Significance testing for methylation.")
parser.add_argument('-c', '--covered', help="DataFrame dump containing total call counts", required=True)
parser.add_argument('-g', '--genome', help="Genome in fasta format.", required=True)

parser.add_argument('-gff', '--gff', help="Genome annotation in gff format.", required=True)
#parser.add_argument('-i', '--index', help="DataFrame with the strand and context data", required=True)

parser.add_argument('--chr_prefix', help="Seqid prefix for chromosomal sequences", default = 'NC_0376')
args= parser.parse_args()


covered = pd.read_pickle(args.covered)
genome = load_genome(args.genome, upper=False)

index = covered.index.to_frame()

strand  = np.empty(len(index), dtype=bool)
context = np.empty(len(index), dtype='U3')
masked  = np.empty(len(index), dtype=bool)

for i, (seqid, pos) in tqdm(enumerate(index.itertuples(index=False)), total= len(index)) : 
        masked[i] = genome[seqid][pos].islower() 
        if genome[seqid][pos] in ('C','c') : 
            strand[i]=True 
            context[i] = str(genome[seqid][pos:pos+3].upper()) 
        else : 
            strand[i]=False 
            context[i] = str(genome[seqid][pos-2:pos+1].reverse_complement().upper()) 

index['strand'] = strand
index['masked'] = masked
index['context'] = context

index['nuclear'] = index.seqid.apply(lambda seqid: seqid.startswith(args.chr_prefix)).astype(bool)


index['di_context'] = index.context.apply(lambda x: x[:2]).astype('category')

db = load_annotation(args.gff)

gene_accession = dict()

for gene in db.all_features(featuretype='gene') :
    for label in gene.attributes['Dbxref'] :
        if label.startswith('GeneID') :
            gene_accession[gene.id] = label.split(':')[1] 

# Choosing one mrna per gene

print(f"Choosing a representative mrna per gene..")

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


# Build dict of feature positions for fast search

feature_pos = dict()

for mrna in representative_mrnas :
    if mrna.seqid not in feature_pos : 
        feature_pos[mrna.seqid] = dict()
    
    GeneID = gene_accession[mrna.attributes['Parent'][0]]
    
    gene_pos = (mrna.start, mrna.end)
    exon_pos = [(exon.start, exon.end) for exon in db.children(mrna, featuretype='exon') ]
    
    feature_pos[mrna.seqid][gene_pos] = (GeneID, exon_pos)

# Function to  search in the feature pos index

def pos_infeature( postuple, feature_pos) :
    seqid, pos = postuple
    
    # offset from 0 based pos to 1 based annotation
    pos +=1 
    
    if seqid not in feature_pos : 
        return (None, None) 
    
    for gene_pos in feature_pos[seqid] : 
        if  gene_pos[0] <=  pos <= gene_pos[1] :
            GeneID, exon_list = feature_pos[seqid][gene_pos]
            # in gene
            for i, exon_pos in  enumerate(exon_list, start=1):
                if  exon_pos[0] <=  pos <= exon_pos[1] :
                    # in exon
                    return (GeneID, i)
            # in gene, but not in exon
            return (GeneID, None)
            
    return (None, None)

# Note that we are neglecting the case where the position could be inexon in ANOTHER gene
# There should not be too many of these.

df = pd.DataFrame([pos_infeature(pos, feature_pos) for pos in tqdm(index.index, total=len(index)) ],
                      index = index.index ,
                      columns=['ingene', 'inexon'], dtype='category') 

index  = index.join(df)

index.drop(['seqid', 'pos'], axis=1, inplace=True)
index.reset_index(inplace=True)

index.seqid = index.seqid.astype('category')

if index.pos.max() < 2**32:
    index.pos = index.pos.astype('uint32')

index.context = index.context.astype('category')
index.di_context = index.di_context.astype('category')
index.strand = index.strand.astype(bool)

index.to_pickle('index.pkl')
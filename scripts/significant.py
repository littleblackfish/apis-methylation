#!/usr/bin/env python

import pandas as pd
import numpy as np
from os import path
from scipy.stats import binom
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests
from argparse import ArgumentParser

test_methods = ['bonferroni', 'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky', 'holm', 'holm-sidak' ]

parser = ArgumentParser(description="Significance testing for methylation.")
parser.add_argument('-m', dest='methylated', help="DataFrame containing methylated call counts", required=True)
parser.add_argument('-c', dest='covered', help="DataFrame containing total call counts", required=True)

parser.add_argument('--conversion_rate', help="(Bisulfite) Conversion rate for significance testing", type=float, default=0.995)
parser.add_argument('--alpha', help="False discovery rate", type=float, default=0.001)
parser.add_argument('--min_covered', help="Minimum coverage", type=int, default=5)
parser.add_argument('--correction_method', help="Multiple testing correction method", default='fdr_bh', choices=test_methods)

args= parser.parse_args()

print (args)
print ('Loading files...')

methylated_fn = path.basename(args.methylated).split('.')[0]
covered_fn = path.basename(args.covered).split('.')[0]

methylated = pd.read_pickle(args.methylated)
covered = pd.read_pickle(args.covered)

assert methylated.index.equals(covered.index)
assert methylated.columns.equals(covered.columns)

print (f'Loaded {len(covered)} C sites.')

print('Filtering coverage..')

# We zero out coverage values below the threshold
covered = covered.where(covered >= args.min_covered, other=0)

# We zero out the 0.1% highest covered positions for each experiment
# (outliers, likely experimental artifacts)
#max_covered = covered.where(covered>0).quantile(0.999)
max_covered = pd.Series({ label:data[data>0].quantile(0.999) for label, data in covered.iteritems() }, 
                        dtype='uint16')


covered = covered.where(covered <= max_covered, other=0)

# Saving us some space
if covered.max().max()<256 :
    print ('Casting to uint8 for space savings.')
    covered = covered.astype('uint8')
    methylated = methylated.astype('uint8')

# This is the boolean mask for significance, same shape as methylated/covered matrices
# True means we reject null hypothesis
significant = pd.DataFrame(np.zeros(methylated.shape, dtype=bool), index=methylated.index, columns=methylated.columns)

failed_conversion_prob = 1 - args.conversion_rate

# Buffer for p values to avoid re-calculating binomial cdf
p_buffer = dict()

print ('Significance testing..')

# Significance testing per column
for (sample, met), (smp,cov) in tqdm(zip(methylated.iteritems(), covered.iteritems()), total=methylated.shape[1]) :	
    assert sample == smp
    
    # We only take positions with some methylated calls
    mask = (cov > 0) & (met > 0)
    p = np.empty(mask.sum())

    # We compare methylated vs. total call count per position 
    # through a binomial model to get p values
    for i , m_c in enumerate(zip(met[mask], cov[mask])) :
        try:
            p[i] = p_buffer[m_c]
        except :
            p[i] = p_buffer[m_c] = 1- binom.cdf(*m_c, failed_conversion_prob)

    # We correct p values for multiple testing
    reject, pvals, ignore, ignore = multipletests(p, alpha=args.alpha, method=args.correction_method)

    # Write results to relevant column
    significant[sample][mask]=reject

# We only take rows where that are significantly methylated in a sample
significant = significant[significant.any(axis=1)]

print(f'{len(significant)} significantly methylated sites.')

# We zero out insgnificant methylation calls (noise)
methylated = methylated.reindex_like(significant).where(significant, other=0)

methylated.to_pickle(f'{methylated_fn}_significant.pkl')
covered.to_pickle(f'{covered_fn}_filtered.pkl')

#covered.reindex_like(significant).to_pickle(f'{covered_fn}_filtered_significant.pkl')
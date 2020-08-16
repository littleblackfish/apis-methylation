####### NOTEBOOK UTILITIES ###########
import pandas as pd 
import numpy as np

# Return the matching F strand index for CpG sites in the R sites
def index_R2F(index_R) : 
    assert type(index_R) == pd.MultiIndex
    tmp = index_R.to_frame() 
    tmp.pos -=1
    return pd.MultiIndex.from_frame(tmp)

# Return the matching R strand index for CpG sites in the F sites
def index_F2R(index_F) : 
    assert type(index_F) == pd.MultiIndex
    tmp = index_F.to_frame() 
    tmp.pos +=1
    return pd.MultiIndex.from_frame(tmp)

# Calculate Spearman's correlation and Wilcoxon test for 
# methylation levels in complementary CpG sites in each sample

def correlation_CG(methylome_index, methylated, covered ) :
    from scipy.stats.stats import spearmanr
    from scipy.stats import wilcoxon
    
    assert methylome_index.index.equals(methylated.index)
    
    # Take only CpG context
    methylome_index_CG = methylome_index[methylome_index.di_context =='CG']
    methylated_CG = methylated[methylome_index.di_context =='CG']
    
    # Split into strands
    methylated_CG_F = methylated_CG[ methylome_index_CG.strand ]
    methylated_CG_R = methylated_CG[~methylome_index_CG.strand ]
    
    print(f'{methylated_CG_F.shape[0]} methylable CpG sites in the (F)orward strand.')
    print(f'{methylated_CG_R.shape[0]} methylable CpG sites in the (R)everse strand.')

    # We report some statistics while we're at it
    
    both_index_F = methylated_CG_F.index.intersection(index_R2F(methylated_CG_R.index))
    either_index_F = methylated_CG_F.index.union(index_R2F(methylated_CG_R.index))
    
    n_either = len(either_index_F)

    print(f'{n_either} destranded CpG sites methylable in either strand.')
    print(f'{len(both_index_F)} ({100*len(both_index_F)/n_either:.1f}%) destranded CpG sites methylable in both strands.')

    only_R_index = methylated_CG_R.index.difference(index_F2R(both_index_F))
    only_F_index = methylated_CG_F.index.difference(both_index_F)    
    n_only_F = len(only_F_index)
    n_only_R = len(only_R_index)
    
    tmp = (methylome_index.reindex(only_F_index).methylatedin == 1).sum()    
    
    print(f'{n_only_F} ({100*n_only_F/n_either:.1f}%) CpG sites only methylable in F ({100*tmp/n_only_F:.1f}% uniquely methylated).')

    tmp = (methylome_index.reindex(only_R_index).methylatedin == 1).sum()    

    print(f'{n_only_R} ({100*n_only_R/n_either:.1f}%) CpG sites only methylable in R ({100*tmp/n_only_R:.1f}% uniquely methylated).')
    
    assert n_only_F +n_only_R +len(both_index_F) == n_either
    
    # Back to business..
    # Reindex each strand to include complement of the other
    
    methylated_CG_F = methylated_CG_F.reindex(either_index_F, fill_value=0)
    methylated_CG_R = methylated_CG_R.reindex(index_F2R(either_index_F), fill_value=0)
    
    # Get coverages for each strand
    
    covered_CG_F = covered.reindex(either_index_F, fill_value=0)
    covered_CG_R = covered.reindex(index_F2R(either_index_F), fill_value=0)
    
    # From methylated read counts to methylation levels
    
    methylation_CG_F = methylated_CG_F / covered_CG_F
    methylation_CG_R = methylated_CG_R / covered_CG_R

    # Shift the R strand's index to match F
    
    methylation_CG_R.index = index_R2F(methylation_CG_R.index)
    
    assert methylation_CG_R.index.equals(methylation_CG_F.index)

    mstats = dict(covered_on_both=dict(), 
                  spearman_r=dict(), 
                  spearman_p=dict(), 
                  wilcoxon_w=dict(), 
                  wilcoxon_p=dict())

    for sample in methylated_CG.columns: 

        df = pd.concat([methylation_CG_F[sample], 
                        methylation_CG_R[sample]], axis=1, join='inner', names=['F','R']).dropna()

        df.columns=['F','R']

        mstats['covered_on_both'][sample]=len(df)

        mstats['spearman_r'][sample], mstats['spearman_p'][sample] = spearmanr(df)
        mstats['wilcoxon_w'][sample], mstats['wilcoxon_p'][sample] = wilcoxon(df.F,df.R)

    return pd.DataFrame(mstats, index=methylation_CG_F.columns)

def destrand_CG(methylome_index, methylated, index, covered):
    
    assert methylome_index.index.equals(methylated.index)
    assert index.index.equals(covered.index)
    
    # Take only CpG context
    methylome_index_CG = methylome_index[methylome_index.di_context =='CG']
    methylated_CG = methylated[methylome_index.di_context =='CG']
    
    # Split into strands
    methylated_CG_F = methylated_CG[ methylome_index_CG.strand ]
    methylated_CG_R = methylated_CG[~methylome_index_CG.strand ]
    
    print(f'{methylated_CG_F.shape[0]} methylable CpG sites in the (F)orward strand.')
    print(f'{methylated_CG_R.shape[0]} methylable CpG sites in the (R)everse strand.')

    either_index_F = methylated_CG_F.index.union(index_R2F(methylated_CG_R.index))
    
    n_either = len(either_index_F)

    print(f'{n_either} destranded CpG sites methylable in either strand.')    # We report some statistics while we're at it
    
    # Reindex each strand to include complement of the other
    
    methylated_CG_F = methylated_CG_F.reindex(either_index_F, fill_value=0)
    methylated_CG_R = methylated_CG_R.reindex(index_F2R(either_index_F), fill_value=0)
    
    # Get coverages for each strand
    
    covered_CG_F = covered.reindex(either_index_F, fill_value=0)
    covered_CG_R = covered.reindex(index_F2R(either_index_F), fill_value=0)

    nmissing_F = (covered_CG_F==0).all(axis=1).sum()
    nmissing_R = (covered_CG_R==0).all(axis=1).sum()

    print(f'{nmissing_F} ({100*nmissing_F/n_either:.1f}%) CpG sites uncovered in F.')
    print(f'{nmissing_R} ({100*nmissing_F/n_either:.1f}%) CpG sites uncovered in R.')

    # Shift the R strand's index to match F
    
    methylated_CG_R.index = index_R2F(methylated_CG_R.index)
    covered_CG_R.index = index_R2F(covered_CG_R.index)
    
    assert methylated_CG_F.index.equals(methylated_CG_R.index)
    assert covered_CG_F.index.equals(covered_CG_R.index)
    assert covered_CG_F.index.equals(methylated_CG_F.index)
    
    # Destrand 

    covered_CG_destranded = covered_CG_F + covered_CG_R
    methylated_CG_destranded = methylated_CG_F + methylated_CG_R
    
    # Update mean methylation levels after binning

    index_CG_destranded = index.reindex(either_index_F).copy()

    index_CG_destranded['coveredin'] = (covered_CG_destranded>0).sum(axis=1)

    index_CG_destranded['methylatedin'] = (methylated_CG_destranded>0).sum(axis=1)

    index_CG_destranded['methylated_ratio'] = index_CG_destranded.methylatedin / index_CG_destranded.coveredin

    index_CG_destranded['mean_methylation'] = methylated_CG_destranded.sum(axis=1)/covered_CG_destranded.sum(axis=1)

    filt_coverage = covered_CG_destranded.where(methylated_CG_destranded>0, other=0)

    index_CG_destranded['cond_mean_methylation'] =  methylated_CG_destranded.sum(axis=1) / filt_coverage.sum(axis=1)
    
    return index_CG_destranded, covered_CG_destranded, methylated_CG_destranded


# Shuffle the order of samples and report discovered methylome size at ith sample

def methylome_size(ismethylated, nshuffle=100) :
    from random import shuffle

    assert (ismethylated.dtypes==bool).all()

    samples = ismethylated.columns.to_list()

    methylome_size = np.zeros((nshuffle, len(samples)))

    for i in range(nshuffle) :
        shuffle(samples)
        
        # Discovered methylable so far
        discovered = np.zeros(len(ismethylated), dtype=bool)
    
        for j, sample in enumerate(samples) :
            discovered = discovered | ismethylated[sample]
            methylome_size[i,j] = discovered.sum()

    df = pd.DataFrame(dict(methylable=methylome_size.mean(axis=0), 
                           err=methylome_size.std(axis=0)),
                      index=np.arange(1,len(samples)+1))
        
    return df
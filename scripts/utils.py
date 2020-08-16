import json, os 
import pandas as pd
import numpy as np

# Loads the json file that holds all the relevant filepaths for a given analysis
def load_paths(filename):
    with open(filename, "r") as f:
        paths = json.load(f)

    prefix = paths["prefix"]

    for key, value in paths.items():
        # Do not prepend prefix if it is obviously absolute or relative path
        if not value.startswith("/") and not value.startswith(".."):
            paths[key] = os.path.join(prefix + value)

    return paths

# Parses master json file to a DataFrame
def parse_experiments(filename) :
    experiments = json.load(open(filename))['experiments']

    for experiment in experiments :
        experiment['samples'] = {sample.pop('id'): sample for sample in experiment['samples']}

    experiments = {experiment.pop('id'):experiment for experiment in experiments}

    tmp = dict()
    for e in experiments :
        for s in experiments[e]['samples'] :
            sample = experiments[e]['samples'][s]

            sample.pop('layout')
            for i,r in enumerate(sample.pop('accession')) :
                tmp[(e,s,i)]=dict(doi=experiments[e]['doi'],
                                  author=experiments[e]['author'],
                                  **sample)

    df = pd.DataFrame.from_dict(tmp, orient='index').sort_index()
    df.index.names=('study', 'sample', 'replicate')
    #assert methylated.columns.equals(df.index)
    return df

# Extracts the sequence for a given feature
def feature_sequence(genome, feature):
    sequence = genome[feature.seqid][feature.start - 1 : feature.end]
    if feature.strand == "+":
        return sequence
    else:
        return sequence.reverse_complement()


def feature_cpgs(feature, df):
    return df[(feature.seqid, feature.start) : (feature.seqid, feature.end)]


# Loads the genome into a dict
def load_genome(genomefile, upper = True):
    from Bio import SeqIO
    from Bio.Alphabet import generic_dna

    print(f"Loading genome assembly from {genomefile}")

    genome = dict()

    for contig in SeqIO.parse(genomefile, "fasta", alphabet=generic_dna):
        seqid = contig.name.strip("|").split("|")[-1]
        if upper :
            genome[seqid] = contig.seq.upper()
        else : 
            genome[seqid] = contig.seq

    return genome


# Loads the annotation db if it exists
# Parses gff3 into db if it doesn't.
def load_annotation(annotationfile):
    import gffutils

    print(f"Loading annotation from {annotationfile}")

    db_path = annotationfile + ".db"

    try:
        db = gffutils.FeatureDB(db_path)
    except:
        db = gffutils.create_db(
            annotationfile,
            # dbfn=':memory:',
            dbfn=db_path,
            force=True,
            force_gff=True,
            id_spec={"gene": "ID", "mRNA": "ID", "exon": "ID"},
        )
    #db.analyze()

    return db


# Merges replicates,
# effectively collapsing the last level of the column index
def merge_replicates(methylation, coverage):

    # Methylated read counts
    methylated_reads = (methylation / 100 * coverage.reindex_like(methylation)).round()

    # Sum all replicates for a sample
    merged_methylated_reads = methylated_reads.groupby(level=(0, 1), axis=1).sum()
    merged_coverage = coverage.groupby(level=(0, 1), axis=1).sum()

    # Go back to percent methylation rates in the merged samples
    merged_methylation = (
        merged_methylated_reads
        / merged_coverage.reindex_like(merged_methylated_reads)
        * 100
    )

    return merged_methylation, merged_coverage


# Iterates over a sequence with a k-width window and stride
def kmer_iter(sequence, k, offset=0, stride=1):
    for i in range(offset, len(sequence) - k + 1, stride):
        yield sequence[i : i + k]

# Returns all possible kmers of size k
def all_kmers(k, alphabet="ATCG"):
    from itertools import product

    kmers = ["".join(kmer) for kmer in product(alphabet, repeat=k)]
    return kmers


# Returns chromosome id for chromosome number
def chr_id(chromosome_no, prefix="NC_0376", offset=38, postfix=".1"):
    assert type(chromosome_no) == int
    n = chromosome_no + offset
    return prefix + str(n) + postfix

# Loads a GloVe embeddding into a dictionary keyed by kmers
def load_glove_dict(path):
    from numpy import asarray

    mapping = dict()
    with open(path) as f:
        for line in f:
            values = line.split()
            word = values[0]
            mapping[word] = asarray(values[1:], dtype="float32")
    return mapping


# Loads a GloVe embedding into 4^k x d matrix
def load_glove_mat(vectors_path):
    from numpy import log2, empty

    glove_dict = load_glove_dict(vectors_path)

    # infer k and d
    k = int(log2(len(glove_dict)) / 2)
    d = len(next(iter(glove_dict.values())))

    glove_mat = empty((4 ** k, d), dtype="float32")

    for i, kmer in enumerate(all_kmers(k)):
        glove_mat[i] = glove_dict[kmer]

    return k, glove_mat


# Parses fasta description lines (as specified in this project)
# id is seqid:start-end
# description is a dict in json format


def parse_fasta_description(record):
    id, data = record.description.split()

    seqid, pos = id.split(":")
    start, end = pos.split("-")

    data = json.loads(data)

    return (seqid, int(start), int(end)), data

# Plotting functions

def plot_cpgs(feature, relativePos=False):

    feature_coverage = feature_cpgs(feature)
    feature_methylation = methylation.reindex(feature_coverage.index, fill_value=0)

    methylation_mean = feature_methylation.mean(axis=1)
    methylation_std = feature_methylation.std(axis=1)

    cpg_pos = array([pos for seqid, pos in feature_coverage.index])
    if relativePos:
        cpg_pos -= db[feature].start

    data = [
        go.Scatter(
            x=cpg_pos,
            y=methylation_mean,
            #      text=list(methylation_mean.index[1]) ,
            error_y=dict(
                type="data",  # value of error bar given in data coordinate
                array=methylation_std,
                visible=True,
            ),
            mode="markers",
        )
    ]
    layout = go.Layout(
        xaxis=dict(title="position"),
        yaxis=dict(title="methylation"),
        title=f"{feature.id} on {feature.seqid}:{feature.start}-{feature.end} {feature.strand}",
    )

    return go.Figure(data=data, layout=layout)


def plot_mrna(feature, df_coverage, df_methylation):

    feature_coverage = feature_cpgs(feature, df_coverage)
    feature_methylation = df_methylation.reindex(feature_coverage.index, fill_value=0)

    cpg_pos = array(
        [pos for seqid, pos in feature_coverage.index] * len(feature_coverage.columns)
    )

    met = pd.Series(
        feature_methylation.values.flatten(order="F"), index=cpg_pos
    ).dropna()

    data = [
        go.Box(
            x=met.index,
            y=met,
            boxpoints="all",
            line_width=0,
            width=0.001,
            marker_size=3,
        )
    ]

    layout = go.Layout(
        xaxis=dict(title="position"),
        yaxis=dict(title="methylation"),
        title=f"{feature.id} on {feature.seqid}:{feature.start}-{feature.end} {feature.strand}",
    )

    layout.shapes = [
        go.layout.Shape(
            type="rect",
            fillcolor="lightgrey",
            layer="below",
            line_width=0,
            y0=0,
            y1=feature_methylation.max().max() + 0.5,
            x0=cds.start,
            x1=cds.end,
        )
        for cds in db.children(feature, featuretype="CDS")
    ]

    return go.Figure(data=data, layout=layout)

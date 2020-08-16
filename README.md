This repository contains code used to re-analyse WGBS data from the honey bee (*Apis mellifera*) and some other insects.
It exists primarily for reproducibility purposes;
while it is not polished for public use, it may provide some useful components for re-use.
Indeed, the structure herein was initially developed for a systematic review of DNA methylation patterns in the honey bee, then re-used to expand analysis into 5 more Hymenopterans, namely:
  * *Bombus terrestris* (buff-tailed bumblebee)
  * *Nasonia vitripennis* (jewel wasp)
  * *Ooceraea biroi* (clonal raider ant)
  * *Harpegnathos saltator* (Indian jumping ant)
  * *Camponotus floridanus* (Florida carpenter ant)


**Data** for each species is curated in json files such as [this one](data/Apis_mellifera/meta.json).

**Primary analysis** of WGBS data was done with [nf-core/methylseq](https://github.com/nf-core/methylseq/).

A number of [scripts](scripts/) exist to aid in heavy lifting. Mainly:
  * [prep_run.py](scripts/prep_run.py) will parse a meta.json file and generate the directory structure to run nf-core/methylseq on all samples. (note that this assumes all fastq files are accessible through filesystem)
  * [methylseq.pbs](scripts/pbs/methylseq.pbs) contains the parameters used for nf-core/methylseq execution within this structure.
  * [submit.sh](scripts/pbs/submit.sh) submits jobs for all samples within such a directory.
  * [merge_bedgraph.py](scripts/merge_bedgraph.py) parses nf-core/methylseq output for all samples into two master matrices that contain entirety of experimental data. Cytosines in rows, samples in columns; total calls in one, methylated calls in the other.
  * [significant.py](scripts/significant.py) filters those matrices for minimum coverage and tests for significant methylation against a null of incomplete bisulfite conversion using a binomial model.
  * [build_index.py](scripts/build_index.py) builds an index containing strand, context and annotation data for every cytosine in the genome to match the master matrices..

**Secondary analysis** is then contained in [notebooks](notebooks/), and is entirely based on the output of last 2 scripts.

A complete, although not minimal [conda environment](environment.yml) that supports the entirety of analysis is provided.

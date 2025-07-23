#!/usr/bin/env python3
# File: /projects/salmonella/work/20210212/source1.py

import os
import datetime
import pandas as pd
import pickle
from rotifer.genome.data import NeighborhoodDF

### Data loading and preparation
# Before execution, the merged2.tsv and ssg.tsv files must be downloaded
# from the Zenodo repository at https://zenodo.org/records/16358274
# and placed in ../data

# Loading genome annotations
if os.path.exists("genome.pkl"):
    genome = pickle.load(open("genome.pkl","rb"))
else:
    # Importing MMseqs clusters: (100% id, 100% cov, clusthash) and (80% cov, 0.001 e-value, easy-clust)
    # ssg: sequence similarity groups
    ssg = pd.read_csv("../data/ssg.tsv", sep="\t")

    # Add protein domain architecture to the genome dataframe
    genome = NeighborhoodDF(pd.read_csv("../data/merged2.tsv", sep="\t"))
    genome = genome.merge(ssg, on='pid', how='left')
    del(ssg)

    # Import T6SS components detected using both Pfam and Rocha's profiles
    t6ss = pd.read_csv("../data/t6ss.acc", sep="\t", names=['pid'])
    genome['t6ss'] = genome.c100i100.isin(genome[genome.pid.isin(t6ss.pid)].c100i100)

    # Save genome annotation
    genome = genome.query('type != "gene"')
    pickle.dump(genome, open("genome.pkl",'wb'))

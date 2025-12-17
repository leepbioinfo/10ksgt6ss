#!/usr/bin/env python3
# File: /projects/salmonella/work/20210212/source1.py

import os
import pandas
import pickle
from rotifer.genome import io as rgio
from rotifer.genome import utils as rgutils
from rotifer.genome.data import NeighborhoodDF

# Globals
MIN_CTERMINAL_LENGTH = 100 # Minimal for putative new C-terminal domains
selected = ['assembly', 'nucleotide', 'c80e3', 'c100i100', 'pid', 'locus', 'plen', 't6ss', 'pfam', 'aravind', 'rocha', 'cdd']

### Data loading and preparation

# Loading genome annotations
genome = pickle.load(open("genome.pkl","rb"))

# Table of protein domain coordinates (architecture2tabke + source column)
if os.path.exists("domain.pkl"):
    domain = pickle.load(open("domain.pkl","rb"))
else:
    # Series for selecting CDS rows: genome DF wont change its number of rows
    cds = (genome.type == "CDS")
    domain = pd.read_csv("/projects/salmonella/work/20210212/merged2.scan.arch.tsv", sep="\t").drop_duplicates()
    domain.rename({domain.columns[0]:'pid'}, axis=1, inplace=True)
    tmp = genome[cds].filter(['c100i100','pid']).drop_duplicates()
    domain = tmp.merge(domain, on="pid", how="right")        # Replace column pid with c100i100
    domain = domain.drop(['pid'], axis=1).drop_duplicates()  #
    tmp = genome[cds].filter(selected).drop_duplicates()
    domain = tmp.merge(domain, on="c100i100", how="left")
    pickle.dump(domain, open("domain.pkl",'wb'))
    del(tmp, cds)

### Analysis: identify toxins


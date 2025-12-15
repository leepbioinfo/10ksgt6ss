#!/usr/bin/env python3
# File: /projects/salmonella/work/20210212/source1.py

import os
import datetime
import pandas
import pickle
from rotifer.genome import io as rgio
from rotifer.genome import utils as rgutils
from rotifer.genome.data import NeighborhoodDF

# Loading genome annotations
if os.path.exists("ssg.pkl"):
    ssg = pickle.load(open("ssg.pkl","rb"))
else:
    # Importing MMseqs clusters: (100% id, 100% cov, clusthash) and (80% cov, 0.001 e-value, easy-clust)
    # Marking all T6SS components
    # ssg: sequence similarity groups
    ssg = pd.read_csv("/projects/salmonella/data/merged2/merged2.c80e3_cluster.tsv", sep="\t", names=['c80e3','c100i100']).drop_duplicates()
    c100i100 = pd.read_csv("/projects/salmonella/data/merged2/merged2.c100i100.tsv", sep="\t", names=['c100i100','pid']).drop_duplicates()
    ssg = ssg.merge(c100i100, on='c100i100', how='left')
    ssg = ssg[['c80e3','c100i100','pid']]
    del(c100i100)

    # Load domain architectures
    f = ['merged2.c100i100.aravind.scan.arch',
         'merged2.c100i100.cdd.hmmsearch.arch',
         'merged2.c100i100.pfam.hmmscan.arch',
         'merged2.c100i100.rocha.hmmsearch.arch']
    d = pd.DataFrame()
    for s in ['pfam','aravind','rocha','cdd']:
        tmp = pd.DataFrame()
        for p in [ x for x in f if s in x ]:
            p = pd.read_csv(p, sep="\t")
            p.rename({'ID':'pid', 'architecture':s}, axis=1, inplace=True)
            tmp = pd.concat([tmp,p])
        if s == "cdd":
            tmp['cdd'] = tmp.cdd.replace({'MIX_III':'MIX'})
        if d.empty:
            d = tmp.filter(['pid',s]).drop_duplicates() 
        else:
            d = d.merge(tmp.filter(['pid',s]).drop_duplicates(), on='pid', how='outer')
    d = ssg.filter(['c100i100','pid']).drop_duplicates().merge(d, on='pid', how='right')
    d = d.drop(['pid'], axis=1).drop_duplicates()
    ssg = ssg.merge(d, on="c100i100", how="left")
    del(tmp)
    del(d, f, p, s)

    # Save genome annotation
    pickle.dump(ssg, open("ssg.pkl",'wb'))

# Save
ssg.to_csv("../data/ssg.tsv", sep="\t", index=False)


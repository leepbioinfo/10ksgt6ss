#!/usr/bin/env python3
# File: /projects/salmonella/work/20210212/source1.py

import os
import pandas
import pickle
from rotifer.genome import io as rgio
from rotifer.genome import utils as rgutils
from rotifer.devel.beta import sequence as rdbs

# Globals
MIN_CTERMINAL_LENGTH = 100 # Minimal for putative new C-terminal domains
selected = ['assembly', 'nucleotide', 'c80e3', 'c100i100', 'pid', 'locus', 'plen', 't6ss', 'pfam', 'aravind', 'rocha', 'cdd']

### Data loading and preparation

# Loading genome annotations
genome = pickle.load(open("genome.pkl","rb"))
df = pd.read_pickle('./10k_T6SS_neighbor_df_jaccard.pk')

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

# Identifying evolved proteins with unlabeled C-terminal domain:
selected = df.query('pfam.fillna("-").str.contains("PAAR|HCP|VgrG|DcrB|DUF4150")', engine= 'python')
selected = selected.aravind.str.split('+').explode().reset_index().drop_duplicates(['index', 'aravind'])
selected = selected.sort_values(['index', 'aravind']).sort_index().drop_duplicates('index', keep='last')
selected = df.loc[selected['index']].c80e3.value_counts()
pt = vndf[vndf.c80e3.isin(selected.index.tolist())].copy()
ptdom = domains[domains.pid.isin(pt.pid)].copy()
cterm = ptdom.sort_values(['pid','end'], ascending=False).drop_duplicates(['pid'])
cterm = cterm.filter(['pid','domain','end','plen'])
cterm = cterm.eval('cterm = plen - end + 1')
cterm = cterm.query('cterm > 50')
cterm = cterm.astype({ x:'int' for x in ['end','plen','cterm'] })
cterm = cterm.rename({'end':'start', 'plen':'end'}, axis=1)

# Collect sequneces and slice c-termini
ctseq = rdbs.sequence(cterm.pid.tolist())
ctseq.df = ctseq.df.merge(cterm.rename({'pid':'id','length':'cterm'}, axis=1), on='id', how='left')
ctseq.df.sequence = ctseq.df.apply(lambda row: row['sequence'][row['start'] - 1:row['end']], axis=1)

# cluster sequences using MMseqs2, register representatives
# for 80% coverage and 0% identity (same as e-value <= 1e-3)
# and store alignments for each cluster
#
# These alignemnets were manually revised, and expanded to 
# generate the final datasets of each ST model. 
cts = dict()
ctseq.add_cluster(80,0,inplace=True)
for x in ctseq.df.c80i0.unique().tolist():
    cts[x] = ctseq.filter(f'c80i0 == "{x}"').align()
    cts[x].to_file(f'{x}.aln')

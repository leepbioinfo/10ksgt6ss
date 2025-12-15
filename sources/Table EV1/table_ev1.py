#!/usr/bin/env python3
# IPython log file

import os
import sys
import pickle
import numpy as np
import pandas as pd
from glob import glob
import rotifer
from rotifer.genome import data as rgd
from rotifer.interval import utils as riu
from rotifer.pandas import functions as rpf
from rotifer.devel.beta import sequence as rdbs
from rotifer.devel.alpha import gian_func as gf
import models as mparser
import sql
from other_functions import *

def find_project_root(target_folder='10ksgt6ss'):
    current = os.path.abspath(os.getcwd())
    while True:
        if os.path.basename(current) == target_folder:
            return current
        parent = os.path.dirname(current)
        if parent == current:
            raise FileNotFoundError(f"'{target_folder}' not found in path hierarchy.")
        current = parent

# Encontra a raiz do projeto
project_root = find_project_root('10ksgt6ss')

# Adiciona a pasta /sources ao sys.path
sources_path = os.path.join(project_root, 'sources')
data_path = os.path.join(project_root, 'data')
sys.path.insert(0, sources_path)

#get_ipython().run_line_magic('logstart', 'source4.py')
rotifer.logger.warning("Starting...")
today = "20240211"
reuse = True # Load some pre-calculated data, instead of recalculting
save = True # False # Do not save output files

# Connect to SQLite3 database: full genome annotation and c100i100 clusters
rotifer.logger.warning("Connecting to the Salmonella database...")

# Load annotated neighborhoods
rotifer.logger.warning("Loading Gian results...")
dfj = pd.read_csv(os.path.join(data_path,"10k_vizinho_novo_df_jaccard.tsv"), sep="\t")

# Load model annotations
rotifer.logger.warning("Loading model annotations...")
#cdd = pd.read_csv("/databases/cdd/cddid_all.tbl", sep="\t", names=["id","basename","name","description","length"])
models = pd.read_excel(os.path.join(data_path,"models.xlsx"))
xlsx = pd.ExcelFile(os.path.join(data_path,"control.xlsx"))
t6ssmap = pd.read_excel(xlsx, "t6ssmap")
toxann = pd.read_excel(xlsx, "toxann")
t6loci = pd.read_excel(xlsx, "t6loci")
modex = pd.read_csv(os.path.join(data_path,"modex.tsv"), names=['modex']).modex.tolist()

# Load architectures
rotifer.logger.warning("Loading cluster, architectures and coordinates...")
info = pd.read_pickle(os.path.join(data_path,"info.pkl"))
coord = pd.read_pickle(os.path.join(data_path,"coord.pkl"))

# Laod clusters
ssg = pd.read_csv(os.path.join(data_path,"ssg.tsv"))

# Load genomic locations
rotifer.logger.warning("Loading genomic locations...")
if reuse and os.path.exists("toxnei.pkl"):
    toxnei = pd.read_pickle("toxnei.pkl")
else:
    toxnei = genome[genome.c100i100.isin(coord.ID.to_list())]
    toxnei.to_pickle("toxnei.pkl")
toxnei['block_id'] = ((toxnei.nucleotide != toxnei.nucleotide.shift(1)) | (toxnei.feature_order - toxnei.feature_order.shift(1) > 10)).cumsum()

# Selecting toxins
rotifer.logger.warning("Identifying toxins...")
toxdoms = set(coord.query('function == "toxin"').basename)
toxins  = coord.ID.isin(coord[coord.basename.isin(toxdoms)].ID)
toxins  = toxins | ((coord.source == "st") & (~coord.basename.isin(["TM","SP"]))) # Make sure all ST models are included
#toxins  = coord[(~faketox) & toxins & (~coord.basename.isin(badmodels))].copy()
toxins  = coord[toxins].copy() # Filter by ID: all functions will remain
toxins.rename({'ID':'c100i100'}, axis=1, inplace=True)
map_column(toxins, toxann,  origin="basename", destination="Toxin", through=['basename','model','Toxin'])
map_column(toxins, t6ssmap, origin="basename", destination="Toxin", through=['basename','model','Toxin'])

# Filtering: excluding models
rotifer.logger.warning("Filtering toxin matches...")
toxins = toxins.query("~(basename.isin(@modex))")
toxins = toxins.query("~(basename == 'Tox-URI1.1' and source == 'aravind')")
#
# Setting up the evalue to avoid toxins models would identify cellular components
toxins = toxins.query("~(basename == 'REase-3.1' and evalue >= 1e-25)")
toxins = toxins.query("~(basename == 'STox_46.2' and evalue > 2.7e-28)")
toxins = toxins.query("~(basename == 'STox_37.1' and evalue > 1.2e-62)")          
toxins = toxins.query("~(basename == 'Peptidase_M23.1' and evalue > 7.2e-34)")
toxins = toxins.query("~(basename == 'Tox-Caspase.1' and evalue > 4.4e-19)")
toxins = toxins.query("~(basename == 'STox_46.1' and evalue > 1.3e-18)")

# Separate T6SS markers and toxins
rotifer.logger.warning("Isolating toxins and T6SS markers...")
toxins.function = np.where(toxins.basename.isin(t6ssmap.query('Evolved == 1').basename), "Fused_to", toxins.function)
toxins.function.replace("toxin","Toxin", inplace=True)
toxins = toxins[toxins.function.isin(['Toxin','Fused_to'])] # This important filter removes noise using domain2architecture's best models!
t6fused = toxins.query('function == "Fused_to"').copy()
toxins = toxins.query('function == "Toxin"').copy()
t6fused = t6fused.merge(ssg, on='c100i100', how='left')
t6fused = t6fused.filter(['c100i100','pid','Toxin','source'], axis=1)
t6fused.rename({'Toxin':'Fused_to','source':'EvolvedStrategy'}, axis=1, inplace=True)
t6fused.drop_duplicates(inplace=True)
t6fused.Fused_to = t6fused.Fused_to.fillna("") + ":" + t6fused.pid.fillna("")

# Finding identical proteins
rotifer.logger.warning("Load identical proteins and identify best toxin models...")
toxins = toxins.merge(ssg, on='c100i100', how='left')
toxins['score'] = toxins.source.map({'st':1, 'aravind':1, 'pfam':1, 'cdd':3, 'rocha':4})
toxins = riu.filter_nonoverlapping_regions(toxins, reference=['pid'], start='start', end='end', criteria={'score':True, 'evalue':True, 'region_length':False})

# Merging
rotifer.logger.warning("Merge genomic locations and best toxin models...")
c = ['c100i100','pid','source', 'basename', 'function', 'Toxin']
toxins = toxins.filter(c).drop_duplicates()
toxins = toxins.merge(toxnei.filter(['pid','locus','assembly','nucleotide','feature_order']).drop_duplicates(), on='pid', how='left')
toxins['nei_c'] = toxins.locus.map(dfj.set_index('locus').nei_c.to_dict())
nei_c = t6loci[t6loci.type.notna()]
for x,y in nei_c.set_index('nei_c').sort_values("type").type.to_dict().items():
    if y not in toxins.columns:
        toxins[y] = np.NaN
    toxins[y] = np.where(toxins['nei_c'] == x, toxins.locus, toxins[y])
nei_c = nei_c.type.drop_duplicates().sort_values().tolist()

# Filter genomes of the 10KSGs project
toxins = toxins[toxins.assembly.str.contains("^FD", regex=True)].copy()

# Analysing Cargo and Evolved status
rotifer.logger.warning("Identifying cargo and evolved toxins...")
toxins['Evolved'] = toxins.pid.isin(t6fused.pid)
toxins['Cargo'] = ~toxins['Evolved']
toxnei['Cargo'] = toxnei.pid.isin(toxins[toxins.Cargo].pid)
tmp = ~coord.ID.isin(toxins.pid) 
tmp = tmp & coord.basename.isin(t6ssmap.query('Cargo == 1').basename)
tmp = coord[tmp].ID.drop_duplicates().tolist()
tmp = sql.fetch_identical(gnc, tmp).pid.drop_duplicates()
toxnei['t6ss'] = toxnei.pid.isin(tmp)
tmp = toxnei.groupby('block_id').agg(Cargo=('Cargo','any'), t6ss=('t6ss','any'), pid=('pid','nunique'))
tmp = tmp.query('Cargo & t6ss').index.to_series()
toxnei['hasCargo'] = toxnei.block_id.isin(tmp)
toxins['Cargo'] = toxins.pid.isin(toxnei[toxnei.Cargo & toxnei.hasCargo].pid)

# Calculate overall distribution
def count_fusedto(x):
    df = x.str.split(":", expand=True)
    if len(df.columns) < 2:
        return ""
    df.rename({0:'c0',1:'c1'}, axis=1, inplace=True)
    df = df.query('c0.notna() and c0 != ""').groupby('c0').c1.nunique()
    df = df.astype(str).reset_index()
    df = df.apply(lambda x: x["c0"] + "(" + x["c1"] + ")", axis=1)
    df = ", ".join(df).replace("c0, c1","")
    return df
rotifer.logger.warning("Calculating overall toxin distribution...")
toxdist = toxins.merge(t6fused, on="pid", how='left')
toxdist['All']     = np.where(toxdist['Evolved'] | toxdist['Cargo'],   toxdist.locus, np.NaN)
toxdist['Evolved'] = np.where(toxdist['Evolved'], toxdist.locus, np.NaN)
toxdist['Cargo']   = np.where(toxdist['Cargo'],   toxdist.locus, np.NaN)
toxdist = toxdist.groupby(['basename','source']).agg(
        Toxin=('Toxin','first'),
        Models=('basename','nunique'),
        Genomes=('assembly','nunique'),
        Loci=('locus','nunique'),
        All=('All','nunique'),
        Evolved=('Evolved','nunique'),
        Cargo=('Cargo','nunique'),
        **{ x: (x,'nunique') for x in nei_c },
        EvolvedStrategy=('EvolvedStrategy',gf.count_series),
        Fused_to=('Fused_to', count_fusedto)
).rename({'source':'Strategy'}, axis=1)
#        Strategy=('source',gf.count_series),
toxdist.reset_index(inplace=True)
toxdist['T6SS'] = np.sum(toxdist.loc[:,nei_c], axis=1)

# Add Toxin annotation columns
tmp = toxann.groupby('Toxin').agg(
        FOLD=('FOLD',lambda x: ", ".join(x.dropna().drop_duplicates().sort_values())),
        Putative_Function=('Putative Function',lambda x: ", ".join(x.dropna().drop_duplicates().sort_values())),
        Activity=('Activity',lambda x: ", ".join(x.dropna().drop_duplicates().sort_values())),
        Target=('Target',lambda x: ", ".join(x.dropna().drop_duplicates().sort_values())),
).reset_index()
toxdist = toxdist.merge(tmp, on="Toxin", how="left")

# Save toxdist
toxdist.sort_values(['Evolved','Cargo'] + nei_c + ['source','Genomes','Toxin'], ascending=False, inplace=True)
if os.path.exists("toxdist.xlsx"):
    os.rename("toxdist.xlsx", "toxdist.old.xlsx") 
toxdist = toxdist[['basename','Toxin','source','FOLD','Genomes','Loci','All','Evolved','Cargo','T6SS','i1','i2','i3','i4b','Fused_to','Putative_Function','Activity','Target','EvolvedStrategy']]
toxdist.rename({'basename':'Model', 'source':'Source'}, axis=1, inplace=True)
toxdist.to_excel("toxdist.xlsx", "toxdist", index=False)

# Save intermediary tables
rotifer.logger.warning("Saving other tables...")
xlsx = pd.ExcelWriter("FDsEvolvedCargo.xlsx")
toxins.to_excel(xlsx, "EvolvedCargo", index=False)
t6fused.to_excel(xlsx,"T6SS domains", index=False)
xlsx.close()

rotifer.logger.warning("Done!")

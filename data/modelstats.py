#!/usr/bin/python3
# coding: utf-8

import pandas as pd
import sys
import os

def find_project_root(target_folder='10ksgt6ss'):
    current = os.path.abspath(os.getcwd())
    while True:
        if os.path.basename(current) == target_folder:
            return current
        parent = os.path.dirname(current)
        if parent == current:
            raise FileNotFoundError(f"'{target_folder}' not found in path hierarchy.")
        current = parent
project_root = find_project_root('10ksgt6ss')
sources_path = os.path.join(project_root, 'sources')
sys.path.insert(0, sources_path)
from working_dfs import j1, j2

phagemodels = pd.read_csv("phagemodels.tsv", header=None)[0].tolist()
vndf = pd.read_table("10k_vizinho_novo_df_jaccard.tsv", low_memory=False)
vndf = vndf.query('assembly.str.startswith("FD")').copy()

t6ss = pd.read_excel("t6ss.xlsx","all")
t6ssdict = t6ss.set_index('model')['replace'].to_dict()

# Processing all models
df = vndf.filter(['nei_c','assembly','block_id','cdd','rocha','pfam','aravind'])
df = df.melt(id_vars=['assembly','nei_c','block_id'], value_vars=['cdd','rocha','pfam','aravind'], value_name='arch')
df.arch = df.arch.fillna("").str.split("+")
df = df.explode('arch')
df['arch2'] = df.arch.replace(t6ssdict)

# Processing phage models
phagedf = df.query('arch in @phagemodels')
phagedf = phagedf.drop_duplicates(['block_id','arch'])
phagedist = phagedf.pivot_table(index='nei_c', columns='arch', values='block_id', aggfunc='nunique', fill_value=0, dropna=False)

# Processing T6SS models
t6ssmodels = t6ss['replace'].drop_duplicates().to_list()
t6ssdf = df.copy()
t6ssdf = t6ssdf.query('arch2 in @t6ssmodels')
t6ssdf = t6ssdf.drop_duplicates(['block_id','arch'])
t6ssdist = t6ssdf.pivot_table(index='nei_c', columns='arch2', values='block_id', aggfunc='nunique', fill_value=0, dropna=False)

# Join phage and t6ss counts
modeldist = t6ssdist.join(phagedist, how='left')
locistats = vndf.groupby('nei_c').agg(genomes=('assembly','nunique'), loci=('block_id','nunique'))
stats = locistats.join(modeldist, how='left').fillna(0).astype(int).reset_index()
stats.insert(1, 'revised', stats.nei_c.map(j1).tolist())
stats.insert(2, 'final', stats.revised.map(j2).tolist())
stats.loc[stats.nei_c.isin([24,30,36]),'final'] = 'Phage / Tailocin'
stats.loc[stats.nei_c.isin([1,6,11]),'final'] = 'Phage / Tailocin'
stats.final.replace("-","Orphan",inplace=True)
stats.sort_values(['final','genomes'], ascending=[False,False], inplace=True)

# Convert to percentual
percent = (stats.set_index('nei_c').drop(['genomes','loci','revised','final'], axis=1).div(stats.loci, axis=0) * 100).reset_index()
percent = stats.filter(['nei_c','revised','final','genomes','loci'], axis=1).merge(percent, on='nei_c', how='left')
percent.insert(5,'T6SS',(percent.set_index('nei_c').loc[:,t6ss.query('core == 1')['replace'].drop_duplicates().tolist()] >= 50).sum(axis=1).to_list())
percent.insert(6,'Phage Tail',(percent.set_index('nei_c').loc[:,phage.query('firstrow == "Phage tail components"').model.drop_duplicates().tolist()] >= 50).sum(axis=1).to_list())
percent.insert(7,'Phage Other',(percent.set_index('nei_c').loc[:,phage.query('firstrow == "Other phage components"').model.drop_duplicates().tolist()] >= 50).sum(axis=1).to_list())

xlsx = pd.ExcelWriter("model_gene_distribution_10ksg.xlsx")
percent.to_excel(xlsx,"distribution",index=False)
stats.to_excel(xlsx,"counts",index=False)
xlsx.close()

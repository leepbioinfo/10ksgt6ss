import yaml
import pandas as pd
import numpy as np
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



data_path = f"{project_root}/data/"

meta = pd.read_excel(f'{data_path}/meta_s3.xlsx')
serovar_dict = meta[['Barcode', 'Calculated Salmonella serovar']].set_index('Barcode')['Calculated Salmonella serovar'].fillna('unknown').to_dict()

with open(f'{data_path}/jaccard_type_classification.yaml', 'r') as f:
        j1 = yaml.load(f, Loader=yaml.SafeLoader)

with open(f'{data_path}/jaccard_group_classification.yaml', 'r') as f:
        j2 = yaml.load(f, Loader=yaml.SafeLoader)

#Loading dict to map toxin to effector type (nuclease, proforming etc...)
with open(f'{data_path}/model2function.yaml', 'r') as f:
        type_dict = yaml.load(f, Loader=yaml.SafeLoader)
        

with open(f'{data_path}/function2color.yaml', 'r') as f:
        type_to_color = yaml.load(f, Loader=yaml.SafeLoader)
        

#df10 = pd.read_pickle(f'{data_path}/10k_vizinho_novo_df_jaccard.pkl')
df10 = pd.read_csv(f'{data_path}/10k_vizinho_novo_df_jaccard.tsv', sep="\t", low_memory=False)
df10 = df10[df10['assembly'].str.startswith('FD')]
df10.nei_c = df10.nei_c.replace(j1)
df10.nei_c = df10.nei_c.replace(j2)


#Adding a  rule to classify the T6SS (It must have at least two identified component)
sf = df10.groupby('block_id').agg(t6 = ('t6ss', 'sum'), nei_c = ('nei_c', 'unique')).explode('nei_c')
sf['keep'] =  np.where(sf.nei_c =='Orphan','Orphan', np.where(sf.t6 >2, sf.nei_c, "Orphan"))
df10.nei_c = df10.block_id.map(sf.keep.to_dict())

# Making a list of genomes by T6SS types
assembly_to_t6ss_types = df10.drop_duplicates(subset=['block_id']).dropna(subset=['nei_c']).query('nei_c in ["i3", "i1", "i2", "i4b"]').groupby('assembly').agg(t6 = ('nei_c', lambda x : " ".join(pd.Series(x).drop_duplicates().sort_values().to_list())))
li3 = assembly_to_t6ss_types.where(lambda x : x =="i3").dropna().index.tolist()
li1 = assembly_to_t6ss_types.where(lambda x : x =="i1").dropna().index.tolist()
li2 = assembly_to_t6ss_types.where(lambda x : x =="i2").dropna().index.tolist()
li1i3 = assembly_to_t6ss_types.where(lambda x : x =="i1 i3").dropna().index.tolist()
li3i4b = assembly_to_t6ss_types.where(lambda x : x =="i3 i4b").dropna().index.tolist()
li1i3i4b = assembly_to_t6ss_types.where(lambda x : x =="i1 i3 i4b").dropna().index.tolist()

# List of genomes with one or two types
l1type = li1 + li2 + li3
l2type = li1i3 + li3i4b
#Creation of the dataframe to analyze effectors distribution
#The FDEvolvedCargo4 table contains the new toxins nomenclature
final = pd.read_excel(f'{data_path}/FDEvolvedCargo5.xlsx')


final.i1 = np.where(final.i1.isin(df10.query('nei_c == "i1"').pid.dropna().drop_duplicates()), final.i1, np.nan)
final.i2 = np.where(final.i2.isin(df10.query('nei_c == "i2"').pid.dropna().drop_duplicates()), final.i2, np.nan)
final.i3 = np.where(final.i3.isin(df10.query('nei_c == "i3"').pid.dropna().drop_duplicates()), final.i3, np.nan)
final.i4b = np.where(final.i4b.isin(df10.query('nei_c == "i4b"').pid.dropna().drop_duplicates()), final.i4b, np.nan)

#Adding info from evolved domains
te = pd.read_excel(f'{data_path}/FDEvolvedCargo5.xlsx', sheet_name=1)
te['cargo'] = te.Fused_to.str.split(':', expand=True)[0]
final['evolved_domain'] = final.pid.map(te.set_index('pid').cargo.to_dict())
final.evolved_domain = final.evolved_domain.replace({"DUF2345" :"VgrG"})


t2 = final[['pid', 'source', 'basename', 'assembly']]
t2 = t2.rename(columns={'assembly': 'genome'})

#Genes coded inside T6SS type specific regions
insidei1 = final.query('i1.notna()')
insidei2 = final.query('i2.notna()')
insidei3 = final.query('i3.notna()')


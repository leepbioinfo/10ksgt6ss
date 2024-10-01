import yaml
import pandas as pd
import numpy as np
from matplotlib.gridspec import GridSpec
# %load fig_2C_velha_e_boa.py
# %load fig_2C.py
from working_dfs import t2,li1, li3,li1i3, li2, type_dict, type_to_color, insidei1, insidei2, insidei3, final
import matplotlib.gridspec as gridspec
from matplotlib.colors import  ListedColormap
import seaborn as sns
import matplotlib.pyplot as plt


# *******Renaming some models names ********
data_path = "/panfs/pan1/proteinworld/People/gian/projects/10k/data/final_files"
with open(f'{data_path}/rename_model_name.yaml', 'r') as f:
    rename_model = yaml.load(f, Loader=yaml.SafeLoader)

t2.basename = t2.basename.replace(rename_model)
#********************************************

df3 = t2
pivot_table = df3.pivot_table(index="genome", columns="basename", aggfunc=lambda x:1, fill_value=0).astype(int)
te = pivot_table['pid'].reset_index()
te2 = te.set_index('genome').apply(lambda x :''.join( x.replace({0: 'a', 1:'b'}).tolist()), axis=1)
te3 = te2.value_counts()
te4 = (te3.reset_index().rename({'index': 'tox'}, axis=1).apply(lambda x : pd.Series(list(x.tox)).where(lambda x: x=='b').dropna().index.tolist(), axis=1).explode().drop_duplicates() +1).sort_values().tolist()
te4.insert(0, 0)
te5 = te.iloc[:, te4]
te6 = te2.where(lambda x : x.isin(te3.index.tolist())).dropna().index.tolist()
te7 = te5.query('genome in @te6')
te8 = te7.set_index('genome').apply(lambda x :''.join( x.replace({0: 'a', 1:'b'}).tolist()), axis=1).value_counts()
t9 = te8.reset_index().rename({'index': 'tox', 0: 'number_of_genomes'}, axis=1).set_index('number_of_genomes').tox.reset_index()
t10 = t9.tox.apply(lambda x: pd.Series(list(x)))
t10.columns = te.iloc[:, te4].columns[1:].tolist()
t10.index = t9.number_of_genomes
t10df = pd.Series(t10.columns, name='model').to_frame().eval('tox_type = model.map(@type_dict)').merge(pd.DataFrame(type_to_color).T, left_on='tox_type', right_index=True).sort_values(['position', 'model'])
counter = 0
for x in t10df.tox_type.unique():
    counter += 100
    x = t10df.query('tox_type == @x').model.to_list()
    t10.loc[:,x] = (t10.loc[:,x ] == 'b').replace({True: counter, False :0})

cm = ListedColormap(["#FFFFFF"] + t10df.color.drop_duplicates().tolist())

t10 = t10.astype(int)
t10 = t10[t10df.model.tolist()]
a4_dims = (25, 50)
o = t10.reset_index().reset_index().iloc[:,0:2].rename({'index':'a'}, axis=1)
bar_values  = o.number_of_genomes.tolist()
fig = plt.figure(figsize=a4_dims)

# Configuração do GridSpec

gs = GridSpec(40, 40, figure=fig)
# Heatmap
ax1 = fig.add_subplot(gs[0:39, :-1])
sns.heatmap(t10, linewidths=0, cmap=cm, cbar=False, ax=ax1)
v = ax1.get_xticks()
xt = np.linspace(v.min(), v.max(), t10.shape[1])
ax1.set_xticks(xt, labels=t10.columns.tolist())
ax1.set_xticklabels(t10.columns.tolist(), rotation=45, horizontalalignment='left', fontsize=6)
ax1.tick_params(labelbottom=False, labeltop=True, labelright=False,labelleft=False, left=False, right=False, top=False, bottom=False)

ax1.set(ylabel=None) 

# Gráfico de barras
ax2 = fig.add_subplot(gs[0:39, -1], sharey=ax1)
bar_positions = np.arange(len(bar_values)) +1
ax2.barh(bar_positions, bar_values, color='gray', height=0.8, align='center')
ax2.set_yticks(bar_positions, labels=bar_values)
ax2.plot(bar_values, bar_positions, color='red', linestyle='--', linewidth=1)


# Remove os frames do gráfico de barras
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.spines['bottom'].set_visible(False)

# Esconde os eixos e ticks do gráfico de barras
ax2.tick_params(labelbottom=False, labeltop=False, labelright=True, labelleft=False, left=False, right=False, top=False, bottom=False, labelsize=6)
ax2.set_xlabel('      # of Genomes', fontsize=6)
ax2.xaxis.set_label_position("top")

# Gráfico de barras na parte inferior
ax3 = fig.add_subplot(gs[-1, :-1], sharex=ax1)
## getting count per genome of each toxin:
tox_genomes = t2.groupby('basename').agg(g = ('genome','nunique')).reset_index().set_index('basename').loc[t10.columns.to_list()]

ax3.bar(np.arange(len(v)),tox_genomes.g.tolist(), color='gray', width=0.8, align='edge')
# Remove os frames do gráfico de barras
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['left'].set_visible(True)
ax3.spines['bottom'].set_visible(True)
ax3.tick_params(labelbottom=False, labeltop=False, labelright=False, labelleft=True, left=True, right=False, top=False, bottom=False, labelsize=6)
ax3.set_ylabel('# of Genomes', fontsize=6)
fig.tight_layout()
plt.savefig("./Figure_S1.pdf")


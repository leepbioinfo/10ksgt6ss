# coding: utf-8
# %load fig_2C.py
from working_dfs import type_to_color, li1, li3, li1i3, t2, type_dict, type_to_color
from matplotlib.colors import  ListedColormap
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def heat_map(df, rows, columns,top_rows=False,group=False, output = "heatmap.pdf", orderby=False, top_columns=False):
    df = df[[rows, columns]]
    pivot_table = df.pivot_table(index=rows, columns=columns, aggfunc=lambda x:1, fill_value=0).astype(int)
    pt = pivot_table.apply(lambda x :''.join( x.replace({0: 'a', 1:'b'}).tolist()), axis=1)
    if group:
        if top_rows:
            te3 = pt.value_counts().head(top_rows)
        else:
            te3 = pt.value_counts()
        te4 = (te3.reset_index().rename({'index': 'tox'}, axis=1).apply(lambda x : pd.Series(list(x.tox)).where(lambda x: x=='b').dropna().index.tolist(), axis=1).explode().drop_duplicates()).sort_values().tolist()
        te5 = pivot_table.iloc[:, te4]
        te6 = pt.where(lambda x : x.isin(te3.index.tolist())).dropna().index.tolist()
        te7 = te5.query('genome in @te6')
        te8 = te7.apply(lambda x :''.join( x.replace({0: 'a', 1:'b'}).tolist()), axis=1).value_counts()
        t9 = te8.reset_index().rename({'index': 'tox', 0: 'number_of_genomes'}, axis=1).set_index('number_of_genomes').tox.reset_index()
        t10 = t9.tox.apply(lambda x: pd.Series(list(x)))
        if top_rows:
            t10.columns = pivot_table.iloc[:, te4].columns.tolist()
        else:
            t10.columns = pivot_table.columns.tolist()

        t10.index = t9.number_of_genomes
    else:
        t10 = pivot_table.replace({0:'a',1:'b'})
    
    if orderby:
        t10 = t10.loc[orderby]

    t10df = pd.Series(t10.columns, name='model').to_frame().eval('tox_type = model.map(@type_dict)').merge(pd.DataFrame(type_to_color).T, left_on='tox_type', right_index=True).sort_values(['position', 'model'])
    counter = 0
    for x in t10df.tox_type.unique():
        counter += 100
        x = t10df.query('tox_type == @x').model.to_list()
        t10.loc[:,x] = (t10.loc[:,x ] == 'b').replace({True: counter, False :0})

    cm = ListedColormap(["#FFFFFF"] + t10df.color.drop_duplicates().tolist())

    t10 = t10.astype(int)
    if top_columns:
        filtered_columns = t10df[t10df.model.isin(pivot_table[t10df.model.tolist()].sum().sort_values(ascending=False).head(top_columns).index.tolist())].model.tolist()
    else :
        filtered_columns =  t10df.model.tolist()
    t10 = t10[filtered_columns]
    a4_dims = (10, 10)
    fig, ax = plt.subplots(figsize=a4_dims)
    ax = sns.heatmap(t10, linewidths=0,cmap=cm, cbar=False)
    v = ax.get_xticks()
    xt = np.linspace(v.min(), v.max(), t10.shape[1])
    ax.set_xticks(xt, labels=t10.columns.tolist())
    ax.set_xticklabels(t10.columns.tolist(), rotation=45, horizontalalignment='left')
    ax.tick_params(labelbottom=False,labeltop=True, labelright=False , left=False, right=False, top=True)
    ax.yaxis.set_label_position("right")
    ax.get_yaxis().set_ticks([])
    plt.savefig(f"./{output}")

from working_dfs import t2,li1, li3,li1i3, li2, type_dict, type_to_color, insidei1, insidei2, insidei3, meta, assembly_to_t6ss_types
from matplotlib.colors import  ListedColormap
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import pandas as pd
import math

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

    
def padding_df(df, how='right'):
    cdf = df.copy()
    c = df.columns
    pad_col_name=[]
    for x in c:
        cdf[x] = cdf[x].fillna('').astype(str)
        w = cdf[x].str.len().max()
        cdf[x] = cdf[x].str.pad(width =w, side=how)
        if how == 'left':
            pad_col_name.append(x.rjust(w))
        else:
            pad_col_name.append(x.ljust(w))

    cdf.columns = pad_col_name
    return cdf


def serovar_fig(ser_list, output='ser_fig.pdf'):

    count = 0
    lines = int(math.ceil(len(ser_list)/7) *2)
    fig = plt.figure( figsize=(10,lines),tight_layout=True)
    gs = gridspec.GridSpec(int(lines/2), 7)
    for x,y in enumerate(ser_list):
        print(y)
        count += 1
        ax = fig.add_subplot(gs[x])
        to_plot = t2.query('serovar in @y').drop_duplicates(subset=["basename", "genome"]).basename.value_counts().reset_index().head()
        to_plot['func'] = to_plot['index'].map(type_dict)
        to_plot = to_plot.merge(pd.DataFrame(type_to_color).T, left_on='func', right_index=True).sort_values(['basename','position'], ascending=[False,True])
        to_count  = t2.query('serovar in @y').eval('t6_type = genome.map(@assembly_to_t6ss_types.t6.to_dict())')
        to_count.t6_type  = to_count.t6_type.fillna('ND')
        title = padding_df(to_count.drop_duplicates('genome').t6_type.value_counts().rename('#genomes').reset_index().rename({'index':'T6SS types'}, axis=1))
        title = title.to_string(index=None, justify='left')
        title = f'{y}\n\n{title}\n'
        tex_position_ref = (to_plot.basename.max() * 2/100)
        plt.bar(to_plot['index'], to_plot.basename, color=to_plot['index'].map(type_dict).map(pd.DataFrame(type_to_color).T.color.to_dict()))
        for i, value in to_plot.iterrows():
            plt.text(i, value.basename + tex_position_ref , str(value.basename), ha='center')
        plt.title(title,loc='left')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xticklabels(to_plot['index'].tolist(),rotation=45, horizontalalignment='right')
        ax.set_xlim(-1,5)
        print(count)

    plt.tight_layout()
    plt.savefig(output)    

def count_series(
        series,
        normalize=False,
        cut_off=False,
        ):
    '''
    Function to flatten a pd.Series.value_counts
    The normalize options is to write the results as frequency.
    The cut_off options is to print results only above a given threashold.
    If the normalize function is True, you should give the cutoff value as pct.
    '''
    import pandas as pd

    flattened_list = []
    s = series.dropna().value_counts(normalize=normalize)

    if cut_off:
        cut_off = cut_off/100
        s = s.where(lambda x: x >= cut_off).dropna()
    for y, z in s.items():
        if normalize:
            flattened_list.append(f'{y}({100 * z:.2f}%)')
        else:
            flattened_list.append(f'{y}({z})')
    return ', '.join(flattened_list)


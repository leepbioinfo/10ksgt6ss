# coding: utf-8
# %load ../../final/fig_2G.py
from working_dfs import t2,li1, li3,li1i3, li2, type_dict, type_to_color, insidei1, insidei2, insidei3, meta, assembly_to_t6ss_types
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import math
from rotifer.devel.alpha import gian_func as gf

#set font size and family for all the plot
font = {'family' : 'Arial',
        'weight' : 'medium',
        'style'  : 'normal',
        'size'   : 5}

matplotlib.rc('font', **font)

ser_map = meta[['Barcode', 'Calculated Salmonella serovar']].set_index('Barcode')['Calculated Salmonella serovar'].fillna('unknown').to_dict()
t2['serovar'] = t2.genome.map(ser_map)

to_fig_2G = ["S. Typhimurium", "S. Enteritidis", "S. Panama","S. Infantis","S. Agona","S. Dublin","S. Typhi"]
#Adding the S. to the serovar names
t2.loc[~t2.serovar.fillna("unknown").str.startswith(("I ", "II","IV", "unk")),'serovar'] = 'S. ' + t2.serovar

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
        to_plot['func'] = to_plot['basename'].map(type_dict)
        to_plot = to_plot.merge(pd.DataFrame(type_to_color).T, left_on='func', right_index=True).sort_values(['count','position'], ascending=[False,True])
        to_count  = t2.query('serovar in @y').eval('t6_type = genome.map(@assembly_to_t6ss_types.t6.to_dict())')
        to_count.t6_type  = to_count.t6_type.fillna('ND')
        title = gf.padding_df(to_count.drop_duplicates('genome').t6_type.value_counts().rename('#genomes').reset_index().rename({'#genomes':'T6SS types'}, axis=1))
        title = title.to_string(index=None, justify='left')
        title = f'{y}\n\n{title}\n'
        tex_position_ref = (to_plot['count'].max() * 2/100)
        plt.bar(to_plot['basename'], to_plot['count'], color=to_plot['basename'].map(type_dict).map(pd.DataFrame(type_to_color).T.color.to_dict()))
        for i, value in to_plot.iterrows():
            plt.text(i, value['count'] + tex_position_ref , str(value['count']), ha='center')
        plt.title(title,loc='left')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xticklabels(to_plot['basename'].tolist(),rotation=45, horizontalalignment='right')
        ax.set_xlim(-1,5)
        print(count)

    plt.tight_layout()
    plt.savefig(output)    


serovar_fig(t2.serovar.value_counts().index.tolist(), 'Supp_2G.pdf')
serovar_fig(to_fig_2G, 'Fig2G.pdf')

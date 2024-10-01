from working_dfs import t2, type_dict, type_to_color, final
#import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib
import pandas as pd
import matplotlib.pyplot as plt

#set font size and family for all the plot
font = {'family' : 'Arial',
        'weight' : 'medium',
        'style'  : 'normal',
        'size'   : 4}

matplotlib.rc('font', **font)

cm = 1/2.54

total = t2.groupby('genome').basename.value_counts().rename('xxx').reset_index().basename.value_counts().reset_index().rename({'index':"Effectors", "basename": "number_genomes"}, axis=1)
inside = final.query('i1.notna() or i2.notna() or i3.notna() or i4b.notna()').groupby('assembly').basename.value_counts().rename('xxx').reset_index().basename.value_counts().reset_index().rename({'index':"Effectors", "basename": "inside"}, axis=1)
in_out_count = total.merge(inside, how='left').fillna(0)
in_out_count['outside'] = in_out_count.number_genomes - in_out_count.inside

to_plot_2b = in_out_count.head(20)
to_plot_2b['toxin_type'] = to_plot_2b.Effectors.map(type_dict)
to_plot_2b['color_map'] = to_plot_2b.toxin_type.map(pd.DataFrame(type_to_color).T.color.to_dict())


fig, ax = plt.subplots()
fig.set_size_inches(7.4 * cm, 4 *cm)

fruits = to_plot_2b.Effectors.tolist()
counts1 = to_plot_2b.inside.tolist()
counts2 = to_plot_2b.outside.tolist()
bar_labels = to_plot_2b.number_genomes.tolist()
bar_colors = to_plot_2b.color_map.tolist()
for i, value in enumerate(bar_labels):
    plt.text(i, value + 100, str(value), ha='center')
ax.bar(fruits, counts1, color=bar_colors,label=bar_labels, bottom=counts2)
ax.bar(fruits, counts2,alpha=0.7,label ="orphan", color=bar_colors)
ax.set_xticklabels(fruits, rotation=45, horizontalalignment='right')
ax.set_ylabel('Total genomes')
plt.minorticks_on()
plt.gca().xaxis.set_minor_locator(plt.NullLocator())
ax.spines[['right', 'top']].set_visible(False)

#creating a custom legend:
legend_dict = pd.DataFrame(type_to_color).T.color.to_dict()
patchList = []
for key in legend_dict:
        data_key = mpatches.Patch(color=legend_dict[key], label=key)
        patchList.append(data_key)

plt.legend(handles=patchList, frameon=False)

plt.tight_layout()
plt.savefig("./Figure_2B.svg")

#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

data_path = "."
output = "figure_1E"
input_table = f'''{os.path.join(data_path,output)}.tsv'''

# Load, concatenate and filter
overlap = pd.read_csv(input_table, sep="\t")
overlap.sort_values(['pshamax','model','pshared','overlap'], ascending=False, inplace=True)

# Start drawinf panels and histogramsi
sns.set(rc={'text.usetex' : True})
sns.set_style(style='white')
joint_kws=dict(color=None)
marginal_kws=dict(bins=20, color='gainsboro', edgecolor='black')
g = sns.jointplot(data=overlap, x="pshared", y='levdiff', ylim=(-47,66), marginal_kws=marginal_kws, joint_kws=joint_kws)
g.fig.set_figwidth(16)
g.fig.set_figheight(8)
g.ax_joint.cla()

# Isolate pshared > 0 and add regression line for this subset
missing = overlap.query('pshared == 0').copy()
found = overlap.query('pshared > 0').copy()
#sns.regplot(x=found.pshared, y=found.levdiff, scatter=False, ci=90, line_kws=dict(color='black', linewidth=0.5, linestyle='dashed'))

# Add scatter plot
sc = g.ax_joint.scatter(found.pshared, found.levdiff, c='black', edgecolor='black', marker="o", alpha=0.7)
g.ax_joint.scatter(missing.pshared, missing.levdiff, edgecolor='black', color='yellow', marker="o", alpha=0.7)
g.ax_joint.scatter(missing.query('model == "STox_15.1"').pshared, missing.query('model == "STox_15.1"').levdiff, edgecolor='black', color='red', marker="o", alpha=0.7)

# Vertical lines at pshared == 25% and pshared == 75%
g.ax_joint.axvline(x=0.25, c="gray", ls=':', linewidth=0.8, zorder=0, clip_on=False)
g.ax_joint.axvline(x=0.75, c="gray", ls=':', linewidth=0.8, zorder=0, clip_on=False)
g.ax_joint.axhline(y=0, c="gray", ls='dotted', linewidth=0.6, zorder=0, clip_on=False)

# Saving
g.ax_joint.set_xlabel('Frequency of shared hits between new and published HMMs', fontsize=12, fontweight='bold')
g.ax_joint.set_ylabel('HMM divergence score', fontsize=12)
g.savefig(f"{output}.svg")
g.savefig(f"{output}.pdf")
plt.clf()

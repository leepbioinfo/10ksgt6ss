#!/usr/bin/env python3

from scipy.stats import norm
from working_dfs import t2,l1type, l2type 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from functions import count_series

#set font size and family for all the plot
font = {'family' : 'Arial',
        'weight' : 'medium',
        'style'  : 'normal',
        'size'   : 5}

matplotlib.rc('font', **font)


t3 = t2.groupby('genome').agg(basename = ('basename', lambda x : count_series(x, count='architecture')), effector_count = ('basename', 'nunique'))

fig, ax1 = plt.subplots()
fig.set_size_inches(4, 3)
t1 = t3.query('genome in @l1type')
t2 = t3.query('genome in @l2type')
linspaced1 = np.linspace(t1.effector_count.min()-1, t1.effector_count.max(), 40)
fit1 = norm.pdf(linspaced1, t1.effector_count.mean(), t1.effector_count.std()) * len(t1)
plt.hist(t1.effector_count, bins=np.arange(20)-0.5, alpha=0.5, label="One T6SS", color="steelblue")
plt.plot(linspaced1, fit1 , '-o', color="steelblue", label="Fit", markersize=3)
plt.legend(loc=2,frameon=False)
ax1.set_ylabel('Number of genomes')
ax1.set_xlabel('Number of effectors per genome')
ax1.set(xticks=range(16), xlim=[0.5, 15.5], yticks=range(250, 1800, 250))
plt.text(0.7, 1600, f'count = {int(t1.effector_count.describe()["count"])}\nmean = {t1.effector_count.describe()["mean"]:.2f} \nstd = {t1.effector_count.describe()["std"]:.2f}', ha='left', color='steelblue')
plt.text(11.7, 1600, f'count = {int(t2.effector_count.describe()["count"])}\nmean = {t2.effector_count.describe()["mean"]:.2f} \nstd = {t2.effector_count.describe()["std"]:.2f}', ha='left', color='darkorange')

ax2 = ax1.twinx()
linspaced2 = np.linspace(t2.effector_count.min()-1, t1.effector_count.max(), 40)
fit2 = norm.pdf(linspaced2, t2.effector_count.mean(), t2.effector_count.std()) * len(t2)
plt.hist(t2.effector_count, bins=np.arange(20)-0.5,  alpha=0.5, label="2 or more T6SS", color="darkorange")
ax2.set_ylabel('Number of genomes')
plt.plot(linspaced2, fit2, '-o', label="fit", color="darkorange", markersize=3)
plt.legend(loc=1,frameon=False)

ax2.set(yticks=range(20, 140, 20))


#plt.fill_between(linspaced1, fit1/13, fit2 , where=((fit1/12 <= fit2) & (linspaced1 >2) ), color='red', alpha=0.1, interpolate=False, label='Above')


ax1.spines[['right', 'top']].set_visible(False)
ax2.spines[['left', 'top']].set_visible(False)
ax2.spines['right'].set_color('darkorange')
ax1.spines['left'].set_color('steelblue')


plt.tight_layout()
plt.savefig('./Old_Figure_2A.svg')



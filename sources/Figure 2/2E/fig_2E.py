#!/usr/bin/env python3

import matplotlib.gridspec as gridspec
import matplotlib.gridspec as gridspec
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

import os
import sys

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
sys.path.insert(0, sources_path)



from working_dfs import t2,li1, li3,li1i3, li2, type_dict, type_to_color, insidei1, insidei2, insidei3

#set font size and family for all the plot
font = {'family' : 'Arial',
        'weight' : 'medium',
        'style'  : 'normal',
        'size'   : 5}

matplotlib.rc('font', **font)

pci1 = insidei1.eval(' func_type = basename.map(@type_dict)').func_type.value_counts(normalize=True).reset_index()
pci1['colors'] = pci1['func_type'].map(pd.DataFrame(type_to_color).T.color.to_dict())
pci2 = insidei2.eval(' func_type = basename.map(@type_dict)').func_type.value_counts(normalize=True).reset_index()
pci2['colors'] = pci2['func_type'].map(pd.DataFrame(type_to_color).T.color.to_dict())
pci3 = insidei3.eval(' func_type = basename.map(@type_dict)').func_type.value_counts(normalize=True).reset_index()
pci3['colors'] = pci3['func_type'].map(pd.DataFrame(type_to_color).T.color.to_dict())

fig = plt.figure(tight_layout=True)
gs = gridspec.GridSpec(1, 3)
ax = fig.add_subplot(gs[0, 0])
ax.set_title('i3')
ax.pie(pci3.proportion,
       colors=pci3.colors,autopct='%1.1f%%', pctdistance=1.25)
ax = fig.add_subplot(gs[0, 1])
ax.set_title('i1')
ax.pie(pci1.proportion,
       colors=pci1.colors,autopct='%1.1f%%', pctdistance=1.25)
ax = fig.add_subplot(gs[0, 2])
ax.set_title('i2')
ax.pie(pci2.proportion,
       colors=pci2.colors,autopct='%1.1f%%', pctdistance=1.25)
plt.savefig('./Figure_2E.pdf')

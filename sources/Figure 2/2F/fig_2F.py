#!/usr/bin/env python3
import os
import sys
import matplotlib.gridspec as gridspec
from venn import venn
from matplotlib.colors import  ListedColormap
import matplotlib.pyplot as plt

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


from working_dfs import t2,li1, li3,li1i3, li2, type_dict, type_to_color, insidei1, insidei2, insidei3, final
cmv = ListedColormap(['#74C3E3','#D184E1', '#FF958E', '#00C469' ])

si1 = set(final[['basename', 'i1', 'i2', 'i3', 'i4b']].set_index('basename').notna().eval('orphan  = (i1 + i2 + i3 + i4b) ==False').astype(int).reset_index().query('i1 ==1').basename)
si2 = set(final[['basename', 'i1', 'i2', 'i3', 'i4b']].set_index('basename').notna().eval('orphan  = (i1 + i2 + i3 + i4b) ==False').astype(int).reset_index().query('i2 ==1').basename)
si3 = set(final[['basename', 'i1', 'i2', 'i3', 'i4b']].set_index('basename').notna().eval('orphan  = (i1 + i2 + i3 + i4b) ==False').astype(int).reset_index().query('i3 ==1').basename)
si4b = set(final[['basename', 'i1', 'i2', 'i3', 'i4b']].set_index('basename').notna().eval('orphan  = (i1 + i2 + i3 + i4b) ==False').astype(int).reset_index().query('i4b ==1').basename)
sorphan = set(final[['basename', 'i1', 'i2', 'i3', 'i4b']].set_index('basename').notna().eval('orphan  = (i1 + i2 + i3 + i4b) ==False').astype(int).reset_index().query('orphan ==1').basename)
vend = {'i1':si1,'i2': si2, 'i3':si3,'orphan': sorphan}
venn(vend, cmap=cmv)
plt.savefig("./Figure_2F.png")

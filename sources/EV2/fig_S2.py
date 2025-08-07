#!/usr/bin/env python3
import os
import sys
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import math

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

from working_dfs import t2,li1, li3,li1i3, li2, type_dict, type_to_color, insidei1, insidei2, insidei3, meta, assembly_to_t6ss_types
from functions import serovar_fig

#set font size and family for all the plot
font = {'family' : 'Arial',
        'weight' : 'medium',
        'style'  : 'normal',
        'size'   : 5}

matplotlib.rc('font', **font)

ser_map = meta[['Barcode', 'Calculated Salmonella serovar']].set_index('Barcode')['Calculated Salmonella serovar'].fillna('unknown').to_dict()
t2['serovar'] = t2.genome.map(ser_map)

#Adding the S. to the serovar names
t2.loc[~t2.serovar.fillna("unknown").str.startswith(("I ", "II","IV", "unk")),'serovar'] = 'S. ' + t2.serovar

serovar_fig(t2.serovar.value_counts().index.tolist(), 'Figure_S2.pdf')

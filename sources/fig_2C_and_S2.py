#!/usr/bin/env python3

from working_dfs import t2,li1, li3,li1i3, li2, type_dict, type_to_color, insidei1, insidei2, insidei3, meta, assembly_to_t6ss_types
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import math
from functions import serovar_fig

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

serovar_fig(t2.serovar.value_counts().index.tolist(), 'Supp_2G.pdf')
serovar_fig(to_fig_2G, 'Fig2G.pdf')

#!/usr/bin/env python3
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
from working_dfs import li1, li3, li1i3, t2, type_dict, type_to_color
from fig_2B_function import heat_map
l = li3+li1i3+li1
l.remove('FD01875449')
l.remove('FD01872508')
heat_map(t2,'genome','basename', orderby=l,top_columns=15, output='Figure_2D.pdf')

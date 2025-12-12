#!/usr/bin/env python3

import os
import sys
import pickle
import numpy as np
import pandas as pd
from glob import glob
import rotifer
from rotifer.genome import data as rgd
from rotifer.interval import utils as riu
from rotifer.pandas import functions as rpf
from rotifer.devel.beta import sequence as rdbs
from rotifer.devel.alpha import gian_func as gf
from rotifer.devel.alpha import collection as rdac
from other_functions import *

# Setup
today = "20240211"
reuse = False

# Load annotated neighborhoods
rotifer.logger.warning("Loading Gian results...")
dfj = pickle.load(open("../10k_vizinho_novo_df_jaccard.pk","rb"))

# Load model annotations
rotifer.logger.warning("Loading model annotations...")
models = pd.read_excel("models.xlsx")

# List of T6SS markers for evolved effectors
rotifer.logger.warning("Loading evolved models...")
xlsx = pd.ExcelFile("control.xlsx")
t6ssmap = pd.read_excel(xlsx, "t6ssmap")
toxann = pd.read_excel(xlsx, "toxann")

# Load seqscan results
rotifer.logger.warning("Loading architectures...")
if reuse and os.path.exists("info.pkl"):
    info = pd.read_pickle("info.pkl")
else:
    arch = pd.concat([ pd.read_csv(x, sep="\t").eval(f'source = "{x.replace(".scan.arch","").replace(today + ".","")}"') for x in glob("*.scan.arch") ])
    arch = arch.pivot(index='ID', columns='source', values='architecture')
    arch = arch.reset_index().rename({'ID':'pid'}, axis=1)
    info = gf.load_seq_scan(today, folder=".", haldane=True)
    info.rename({'c80e3':'c80i0'}, axis=1, inplace=True)
    info = info.drop(['pfam','aravind'], axis=1).merge(arch, on='pid', how='left')
    info.sort_values(['c80i0','c80i70','c100i100'], inplace=True)
    info["function"] = np.NaN
    info.to_pickle("info.pkl")
    del(arch)

# Load HMMscan coordinates
if reuse and os.path.exists("coord.pkl"):
    coord = pd.read_pickle("coord.pkl")
else:
    rotifer.logger.warning("Loading libraries...")
    coord = "st aravind cdd rocha pfam".split(" ")
    coord = [ today + "." + x + ".scan.arch.tsv" for x in coord ]
    coord = load_library(coord, prefix=f'{today}.', suffix=".scan.arch.tsv")
    map_column(coord, models, origin='domain', destination='basename', subsets='source')
    coord.loc[coord.domain == "SP","basename"] = "SP"
    coord.loc[coord.domain == "TM","basename"] = "TM"
    coord.loc[coord.domain == "FD01875780_04639_sg1048","basename"] = "FD01875780_04639"
    coord.loc[coord.basename == "SP","function"] = "localization"
    coord.loc[coord.basename == "TM","function"] = "localization"
    map_column(coord, models, origin='basename', destination='function', through=['basename'])
    coord.to_pickle("coord.pkl")

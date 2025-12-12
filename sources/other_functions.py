import os
import sys
import pickle
import sqlite3
import numpy as np
import pandas as pd
from glob import glob
import rotifer
from rotifer.genome import data as rgd
from rotifer.interval import utils as riu
from rotifer.pandas import functions as rpf
from rotifer.devel.beta import sequence as rdbs
from rotifer.devel.alpha import gian_func as gf
from rotifer.devel.alpha import rodolfo as rdar
from rotifer.devel.alpha import collection as rdac
import models as mparser
import sql

# Set basename
def map_column(df, mapdf, origin='model', destination='basename', through=['basename','name','id','acc'], subsets='source', inplace=True):
    if destination not in df.columns:
        df[destination] = np.NaN
    if subsets in df.columns:
        subsetNames = df[subsets].drop_duplicates()
    else:
        subsetNames = [None]
    for setName in subsetNames:
        missing = True
        mapsubset = mapdf
        if setName:
            missing = (df[subsets] == setName)
            if setName in set(mapdf[subsets].drop_duplicates()):
                mapsubset = mapdf[mapdf[subsets] == setName]
        for column in through:
            missing = missing & df[destination].isna()
            if missing.any():
                if column == destination:
                    dmap = mapsubset.eval(f'replace = {column}').set_index(column)['replace'].to_dict()
                else:
                    dmap = mapsubset.set_index(column)[destination].to_dict()
                update = df.loc[missing, origin].map(dmap).tolist()
                df.loc[missing,destination] = update
    del(column)
    del(mapsubset)
    del(missing)
    del(setName)
    if not inplace:
        return df

def load_library(files, prefix=None, suffix=None, source_method=None):
    coord = []
    for table in files:
        if source_method:
            source = source_method(table)
        else:
            source = os.path.basename(table)
            if prefix:
                source = source.replace(f'{prefix}',"")
            if suffix:
                source = source.replace(f'{suffix}',"")
        if not os.path.exists(table):
            rotifer.logger.warning(f'No file {table} for libray {source}')
            continue
        library = pd.read_csv(table, sep="\t")
        rotifer.logger.warning(f'Library {table} loaded')
        library['source'] = source
        coord.append(library)
    coord = pd.concat(coord, ignore_index=True)
    return coord

def compare_hmmsearch_hmmscan(hmmsearch, hmmscan, models):
    result = hmmsearch.groupby(['basename','source']).sequence.nunique().reset_index()
    df2 = hmmscan.groupby(['basename','source']).ID.nunique().reset_index()
    result = result.merge(df2, on=['basename','source'], how='outer')
    result = result.fillna(0).astype({'sequence':int,'ID':int})
    result['diff'] = result.ID - result.sequence
    df2 = models.filter(['basename','source','function','description']).drop_duplicates()
    result = result.merge(df2, on=['basename','source'], how='left')
    result.sort_values(['diff','ID','sequence'], ascending=False, inplace=True)
    result.reset_index(drop=True, inplace=True)
    result['description'] = rpf.pad(result.description, side="right")
    del(df2)
    return result


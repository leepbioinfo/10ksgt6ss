#!/usr/bin/env python3
# File: /projects/salmonella/work/20210212/source1.py

import os
import datetime
import pandas
import pickle
from rotifer.genome import io as rgio
from rotifer.genome import utils as rgutils
from rotifer.genome.data import NeighborhoodDF

# Globals
_o = dir()
TODAY = datetime.datetime.now().strftime("%Y%m%d")
MIN_CTERMINAL_LENGTH = 100 # Minimal for putative new C-terminal domains
selected = ['assembly', 'nucleotide', 'c80e3', 'c100i100', 'pid', 'locus', 'plen', 't6ss', 'pfam', 'aravind', 'rocha', 'cdd']

# Aggregator for list with NaN
def scat(x):
    return ", ".join(sorted(np.unique(", ".join(sorted([ x for x in { str(y):1 for y in x if isinstance(y,str) }.keys() if x != "" ])).split(", "))))

### Data loading and preparation

# Load manually curated toxin annotation (prepared from ../20200803/domdist.tsv)
if os.path.exists("/projects/salmonella/work/20210212/toxann.tsv"):
    toxann = pd.read_csv("/projects/salmonella/work/20210212/toxann.tsv", sep="\t")
else:
    toxann = pd.read_csv("/home/leep/rfsouza/projects/salmonella/work/20200905/salmonella_t6ss_toxins.txt", sep="\t")
    toxann = toxann.drop('comment', axis=1)  # Toxic effector domain's annotations
    toxann = pd.concat([toxann,pd.DataFrame([['Peptidoglycan hydrolase','--','DUF1131','DUF1131','pfam','--',1]], columns=list(toxann.columns))])
    toxann.reset_index(drop=True, inplace=True)
    toxann.to_csv("toxann.tsv", sep="\t", index=False)

# Load list of toxin models
if os.path.exists("/projects/salmonella/work/20210212/toxsig.lst"):
    toxsig = pd.read_csv("/projects/salmonella/work/20210212/toxsig.lst", sep="\t")
    toxsig = set(toxsig.name)
else:
    toxsig = pd.read_csv("/home/leep/rfsouza/projects/toxins/data/Pfam-A.20210208.txt", sep="\t")
    toxsig = toxsig.query('Toxin == 1 and name != "RHS_repeat" and name != "AAA_21"').sort_values(['name'])
    toxsig = pd.concat([toxsig.filter(["name"]), toxann.filter(["model"]).rename({"model":"name"}, axis=1)]).drop_duplicates()
    toxsig.to_csv("toxsig.lst", sep="\t", index=False)
    toxsig = set(toxsig.name) # List of toxins

# Loading genome annotations
if os.path.exists("genome.pkl"):
    genome = pickle.load(open("genome.pkl","rb"))
else:
    # Importing MMseqs clusters: (100% id, 100% cov, clusthash) and (80% cov, 0.001 e-value, easy-clust)
    # Marking all T6SS components
    # ssg: sequence similarity groups
    ssg = pd.read_csv("/projects/salmonella/data/merged2/merged2.c80e3_cluster.tsv", sep="\t", names=['c80e3','c100i100']).drop_duplicates()
    c100i100 = pd.read_csv("/projects/salmonella/data/merged2/merged2.c100i100.tsv", sep="\t", names=['c100i100','pid']).drop_duplicates()
    ssg = ssg.merge(c100i100, on='c100i100', how='left')
    ssg = ssg[['c80e3','c100i100','pid']]
    del(c100i100)

    # Load domain architectures
    f = ['/projects/salmonella/work/20210212/merged2.c100i100.new.aravind.scan.arch',
         '/projects/salmonella/work/20210212/merged2.c100i100.new.cdd.hmmsearch.arch',
         '/projects/salmonella/work/20210212/merged2.c100i100.new.pfam.hmmscan.arch',
         '/projects/salmonella/work/20210212/merged2.c100i100.new.rocha.hmmsearch.arch',
         '/projects/salmonella/work/20200813/mergednrdb_rep.aravind.arch',
         '/projects/salmonella/work/20200803/mergednrdb_rep.cdd.arch',
         '/projects/salmonella/work/20200813/mergednrdb_rep.pfam.arch',
         '/projects/salmonella/work/20200803/mergednrdb_rep.rocha.arch']
    d = pd.DataFrame()
    for s in ['pfam','aravind','rocha','cdd']:
        tmp = pd.DataFrame()
        for p in [ x for x in f if s in x ]:
            p = pd.read_csv(p, sep="\t")
            p.rename({'ID':'pid', 'architecture':s}, axis=1, inplace=True)
            tmp = pd.concat([tmp,p])
        if s == "cdd":
            tmp['cdd'] = tmp.cdd.replace({'MIX_III':'MIX'})
        if d.empty:
            d = tmp.filter(['pid',s]).drop_duplicates() 
        else:
            d = d.merge(tmp.filter(['pid',s]).drop_duplicates(), on='pid', how='outer')
    d = ssg.filter(['c100i100','pid']).drop_duplicates().merge(d, on='pid', how='right')
    d = d.drop(['pid'], axis=1).drop_duplicates()
    ssg = ssg.merge(d, on="c100i100", how="left")
    del(tmp)
    del(d, f, p, s)

    # Add protein domain architecture to the genome dataframe
    genome = NeighborhoodDF(pd.read_csv("/projects/salmonella/data/merged2/merged2.tsv", sep="\t"))
    genome = genome.merge(ssg, on='pid', how='left')
    del(ssg)

    # Import T6SS components detected using both Pfam and Rocha's profiles
    t6ss = pd.read_csv("/projects/salmonella/work/20200803/mergednrdb_rep.t6ss.acc", sep="\t", names=['pid'])
    t6ss = pd.concat([t6ss,pd.read_csv("/projects/salmonella/work/20210212/merged2.c100i100.new.t6ss.acc", sep="\t", names=['pid'])])
    genome['t6ss'] = genome.c100i100.isin(genome[genome.pid.isin(t6ss.pid)].c100i100)

    # Save genome annotation
    genome = genome.query('type != "gene"')
    pickle.dump(genome, open("genome.pkl",'wb'))

# Table of protein domain coordinates (architecture2tabke + source column)
if os.path.exists("domain.pkl"):
    domain = pickle.load(open("domain.pkl","rb"))
else:
    # Series for selecting CDS rows: genome DF wont change its number of rows
    cds = (genome.type == "CDS")
    domain = pd.read_csv("/projects/salmonella/work/20210212/merged2.scan.arch.tsv", sep="\t").drop_duplicates()
    domain.rename({domain.columns[0]:'pid'}, axis=1, inplace=True)
    tmp = genome[cds].filter(['c100i100','pid']).drop_duplicates()
    domain = tmp.merge(domain, on="pid", how="right")        # Replace column pid with c100i100
    domain = domain.drop(['pid'], axis=1).drop_duplicates()  #
    tmp = genome[cds].filter(selected).drop_duplicates()
    domain = tmp.merge(domain, on="c100i100", how="left")
    pickle.dump(domain, open("domain.pkl",'wb'))
    del(tmp, cds)

### Analysis: find target loci and neighbors
if os.path.exists("cgrf.pkl"):
    cgrf = pickle.load(open("cgrf.pkl","rb"))
else:
    cgrf = genome.neighbors(genome.t6ss, after=5, before=5, min_block_distance=2)
    pickle.dump(cgrf, open("cgrf.pkl","wb"))

### Analysis: identify toxins

# Load table that maps domains to T6SS markers
if os.path.exists("/projects/salmonella/work/20210212/t6ssmap.tsv"):
    t6ssmap = pd.read_csv("/projects/salmonella/work/20210212/t6ssmap.tsv", sep="\t")
else:
    t6ssmap = pd.read_csv("/home/leep/rfsouza/projects/salmonella/work/20200804/t6ss.map.txt", sep="\t")
    for rhs in domain[domain.domain.fillna("-").str.contains("rhs", case=False)].domain.unique():
            t6ssmap = t6ssmap.append(pd.DataFrame([["pfam",rhs,"RHS",0,"robson"]], columns=t6ssmap.columns))
    t6ssmap = t6ssmap.query('model != "RHSP-X2"')
    t6ssmap.to_csv("t6ssmap.tsv", sep="\t", index=False)
    del(rhs)

# Selecting putative effectors: fusions to T6SS-related domains and BastionX results
if os.path.exists("toxins.pkl"):
    toxins = pickle.load(open("toxins.pkl","rb"))
else:
    toxarch = []
    for m in ["PAAR","VgrG","RHS","Hcp","MIX"]:
        l = [m] + t6ssmap[t6ssmap['replace']==m].model.tolist()
        v = domain[domain.domain.isin(l)].c80e3.unique()
        v = domain[domain.c80e3.isin(v)]
        v['source'] = v.source.replace({'aravind':'A','pfam':'P','cdd':'C','rocha':'R'})
        v['fused'] = m
        v["final_score"] = 0
        toxarch.append(v)
    del(m,l,v)

    # Load BastionX T6SS results and make sure all proteins within the same
    # cluster are assigned the same BastionX prediction. Keep in mind that
    # only representatives of each c80e3 cluster have been analyzed and we,
    # therefore, keep only the best BastionX score in each c80e3 cluster.
    bx = pd.read_csv("/projects/salmonella/work/20200820/bastionx.tsv", sep="\t").drop_duplicates()
    bx.rename({bx.columns[0]:'pid'}, axis=1, inplace=True)
    # Replace pid by c80e3
    bx = genome[genome.pid.isin(bx.pid)].filter(['c80e3','pid']).drop_duplicates().merge(bx, on="pid", how="right")
    bx = bx.drop(['pid'], axis=1).sort_values(['c80e3','final_score'], ascending=[True,False]).drop_duplicates("c80e3")
    # Merge with genome annotation and filter columns to the match toxins dataframe
    bx = genome[genome.c80e3.isin(bx.c80e3)].filter(selected).merge(bx, on="c80e3", how="right") # Prapagate BastionX prediction
    bx6 = bx.query('final_label == "VI"').filter([*selected, "final_score"])
    bx6['domain'] = np.NaN
    bx6['start'] = 0
    bx6['end'] = 0
    bx6['evalue'] = 10000
    bx6['source'] = 'B'
    bx6['fused'] = np.NaN
    bx6 = bx6[list(toxarch[0].columns)]
    del(bx)

    # Aggregate toxin annotations
    toxarch.append(bx6)
    toxarch = pd.concat(toxarch)
    toxarch.loc[toxarch['start'].isna(),'start'] = 0
    toxarch.loc[toxarch['end'].isna(),'end'] = 0
    del(bx6)

    # Find toxins
    toxins = toxarch.groupby(['c80e3', 'c100i100'])
    toxins = toxins.agg({'assembly':'nunique','locus':'nunique','pfam':'max','aravind':'max','rocha':'max','cdd':'max','fused':scat,'source':scat,'final_score':'max'})
    toxins = toxins.reset_index().sort_values(['fused','c80e3','locus','pfam'], ascending=[True,True,False,True])
    toxins['bastionx'] = toxins.source.str.contains('B')

    # Find proteins with uncharacterized C-terminal regions and tagging known toxins
    #
    # Note that, in our previous analysis (../20200803/source2.py), proteins 
    # containing previously characterized toxin domains were removed,
    # regardless of whether they contained C-terminal extensions or not
    #
    # This time, I'll annotate all uncharacterized C-terminal extensions
    # We also use this oportunity to mark sequences with known toxin domains
    #
    # See ../20200815 and ~rfsouza/projects/salmonalla/work/20200815 for first attempt
    #tmp = toxins.c100i100.unique() # Focus on the fusions to T6SS markers
    #ctcoord = domain.query('c100i100 in @tmp')
    #ctcoord = ctcoord.filter(['c100i100', 'domain', 'start', 'end', 'source', 'plen'])
    ctcoord = toxarch.filter(['c100i100', 'domain', 'start', 'end', 'source', 'plen'])
    ctcoord.drop_duplicates(inplace=True)
    ctcoord['toxin'] = ctcoord.domain.isin(toxsig)
    ctcoord.start.fillna(0, inplace=True)
    ctcoord.end.fillna(0, inplace=True)
    ctcoord.rename({'plen':'ctend'}, axis=1, inplace=True)
    ctcoord.eval('ctstart = end + 1', inplace=True)
    ctcoord = ctcoord.groupby(['c100i100','source']).agg({'toxin':'max','ctstart':'max','ctend':'first'})
    ctcoord.eval('ctlen = ctend - ctstart + 1', inplace=True)
    ctcoord.reset_index(inplace=True)
    ctcoord = ctcoord.groupby(['c100i100']).agg({'toxin':'max','ctstart':'max','ctend':'max','ctlen':'min'}).reset_index()
    toxins = toxins.merge(ctcoord, on='c100i100', how='left').reset_index(drop=True)
    del(ctcoord)

    # Merge Andre's annotation
    andre = pd.read_csv("/home/leep/rfsouza/projects/salmonella/work/20201013/Analysis_Unknowns_Domains_HHpred_Searchs.txt", sep="\t")
    toxins = toxins.merge(andre.drop(['c8e3'], axis=1), on="c100i100", how="left")
    toxins.Hits = toxins.Hits.fillna(0)
    toxins['Probability(%)'] = toxins['Probability(%)'].fillna(0)
    del(andre)

    # Save
    pickle.dump(toxins, open("toxins.pkl","wb"))

# Calculating global toxin prevalence
if os.path.exists("toxdist.tsv"):
    toxdist = pd.read_csv("toxdist.tsv", sep="\t")
else:
    istoxin = set(domain[domain.domain.isin(toxsig)].c80e3)
    istoxin = set(toxins.c80e3.unique()).union(istoxin)
    istoxin = domain.c80e3.isin(istoxin)
    toxdist = domain[istoxin].filter(['c80e3','c100i100','assembly','locus','domain','source']).drop_duplicates()
    toxdist.rename({'domain':'model'}, axis=1, inplace=True)
    toxdist = toxdist.merge(toxann, left_on=['model','source'], right_on=['model','source'], how="left")
    toxdist.drop(['source'], axis=1, inplace=True)
    toxdist.rename({'domain':'toxin'}, axis=1, inplace=True)
    toxdist = toxdist.merge(toxins[['c100i100','fused','source']].drop_duplicates(), on="c100i100", how="left")
    toxdist = toxdist.groupby(['Activity','toxin'])
    toxdist = toxdist.agg({'assembly':'nunique','c80e3':'nunique','c100i100':'nunique','locus':'nunique','fused':scat,'source':scat}).reset_index()
    toxdist = toxdist.sort_values(['Activity','assembly','locus'], ascending=[True,False,False])

    # Calculating toxin prevalence around T6SS loci or related genes
    t6 = cgrf.pid.dropna().unique()
    t6 = istoxin & domain.pid.isin(t6)
    t6 = domain[t6].filter(['c80e3','c100i100','assembly','pid','domain','source'])
    t6.rename({'domain':'model'}, axis=1, inplace=True)
    t6 = t6.merge(toxann, left_on=['model','source'], right_on=['model','source'], how="left")
    t6.drop(['source'], axis=1, inplace=True)
    t6.rename({'domain':'toxin'}, axis=1, inplace=True)
    t6 = t6.groupby(['toxin'])
    t6 = t6.agg({'assembly':'nunique','c80e3':'nunique','c100i100':'nunique','pid':'nunique'}).reset_index()
    # 
    toxdist = toxdist.merge(t6, left_on=['toxin'], right_on=['toxin'], how='outer').set_index('toxin').reset_index()
    toxdist[toxdist.select_dtypes(float).columns] = toxdist.select_dtypes(float).fillna(0).astype(int)
    toxdist['assembly_yx'] = 100 * (toxdist.assembly_y / toxdist.assembly_x)
    toxdist['assembly_yx'] = toxdist['assembly_yx'].fillna(0).round(0).astype(int)
    toxdist.to_csv("toxdist.tsv", sep="\t", index=False)
    del(t6,istoxin)

# Load and identify selected toxins and genomes (curated dataframe)
#toxsel = pd.read_csv("/home/aluizsilva/projects/salmonella/data/workcopy_toxins.20201026.2_ethel_28.tsv", sep="\t")
if os.path.exists("curated.tsv"):
    curated = pd.read_csv("curated.tsv", sep="\t")
else:
    toxsel = pd.read_csv("/home/ebsantos/aluizsilva/projects/salmonella/data/20210118/workcopy_toxins.20201026.2_ethel_28.tsv", sep="\t")
    ethel = [2, 3, 5, 10, 11, 12, 15, 17, 21, 27, 28]
    choices = toxsel[toxsel["Group #"].isin(ethel)]
    curated = choices.filter(["Group #","c80e3","c100i100"]).drop_duplicates()
    curated.rename({'Group #':'group','c80e3':'oldc80e3','c100i100':'pid'}, axis=1, inplace=True)
    tmp = curated.filter(["group","oldc80e3","oldc80e3"]).drop_duplicates()
    tmp.columns = curated.columns
    curated = pd.concat([curated,tmp], ignore_index=True).drop_duplicates()
    curated.to_csv("curated.tsv", sep="\t", index=False)
    del(choices, toxsel, ethel, tmp)

# Import genome metadata
if os.path.exists("md10k.tsv"):
    md10k = pd.read_csv("md10k.tsv", sep="\t")
else:
    md10k = pd.read_csv("/projects/salmonella/data/10k_genomes/data/metadata.tsv", sep="\t")
    md10k.rename({'Barcode':'assembly','Species':'organism'}, axis=1, inplace=True)
    newstrains = (genome.type == "source")
    newstrains = newstrains & genome.assembly.str.contains("GC[AF]_", regex=True)
    newstrains = genome[newstrains].filter(['assembly','nucleotide','organism','end'])
    newstrains = newstrains.groupby(['assembly','organism','nucleotide']).agg({'end':'max'}).reset_index()
    newstrains = newstrains.groupby('assembly').agg({'nucleotide':'nunique','organism':'first','end':'sum'}).reset_index()
    newstrains.rename({'nucleotide':'Num_contigs','end':'Assembly_len'}, axis=1, inplace=True)
    newstrains['Host'] = 'Human'
    newstrains['Country'] = 'Brazil'
    newstrains['Contact'] = 'Ethel'
    md10k = pd.concat([md10k, newstrains], ignore_index=True)
    del(newstrains)

    # Analysis: number of toxins per genome

    # Count T6SS components
    t6loci = cgrf.groupby(['assembly','block_id']).agg({'locus':'nunique','query':'sum'}).reset_index()
    t6loci = t6loci.groupby(['assembly']).agg({'block_id':'nunique','locus':'sum','query':['sum','max']})
    t6loci.reset_index(inplace=True)
    t6loci.columns = ['assembly','number_of_t6ss_loci','number_of_genes_per_t6ss_loci','t6ss_components_total','t6ss_components_max_per_loci']
    md10k = md10k.merge(t6loci, on="assembly", how="left")
    md10k[t6loci.columns[1:]] = md10k[t6loci.columns[1:]].fillna(0).astype(int)
    del(t6loci)

    # Count any toxin
    istoxin = set(domain[domain.domain.isin(toxsig)].c80e3)
    istoxin = set(toxins.c80e3.unique()).union(istoxin)
    istoxin = genome.c80e3.isin(istoxin)
    toxperg = genome[istoxin].groupby(['assembly']).agg({'locus':'nunique'}).reset_index()
    toxperg.rename({'locus':'all_toxins'}, axis=1, inplace=True)
    md10k = md10k.merge(toxperg, on='assembly', how='left')
    md10k.all_toxins.fillna(0, inplace=True)
    md10k.all_toxins = md10k.all_toxins.astype(int)

    # Count toxins near genes with T6SS signatures
    t6 = cgrf.pid.dropna().unique()
    toxperg = genome[istoxin & genome.pid.isin(t6)].groupby(['assembly']).agg({'locus':'nunique'}).reset_index()
    toxperg.rename({'locus':'toxins_near_t6ss'}, axis=1, inplace=True)
    md10k = md10k.merge(toxperg, on='assembly', how='left')
    md10k.toxins_near_t6ss.fillna(0, inplace=True)
    md10k.toxins_near_t6ss = md10k.toxins_near_t6ss.astype(int)
    del(t6, istoxin, toxperg)

    # Count selected toxins
    tmp = genome[genome.pid.isin(curated.pid)].c80e3.unique()
    toxperg = genome[genome.c80e3.isin(tmp)].groupby(['assembly']).agg({'locus':'nunique'}).reset_index()
    toxperg.rename({'locus':'curated_anywhere'}, axis=1, inplace=True)
    md10k = md10k.merge(toxperg, on='assembly', how='left')
    md10k.curated_anywhere.fillna(0, inplace=True)
    md10k.curated_anywhere = md10k.curated_anywhere.astype(int)
    toxperg = cgrf[cgrf.c80e3.isin(tmp)].groupby(['assembly']).agg({'locus':'nunique'}).reset_index()
    toxperg.rename({'locus':'curated_near_t6ss'}, axis=1, inplace=True)
    md10k = md10k.merge(toxperg, on='assembly', how='left')
    md10k.curated_near_t6ss.fillna(0, inplace=True)
    md10k.curated_near_t6ss = md10k.curated_near_t6ss.astype(int)
    del(tmp)

    # Save md10k
    import openpyxl
    md10k.to_csv("md10k.tsv", sep="\t", index=False)
    md10k.to_excel("md10k.xlsx", sheet_name="curated", index=False)

# Manually selected Moreira's genomes
if os.path.exists("moreira.tsv"):
    moreira = pd.read_csv("moreira.tsv", sep="\t")
else:
    moreira = """FD01848829
    FD01848851
    FD01848865
    FD01848866
    FD01848874
    FD01848875
    FD01848827
    FD01848834
    FD01848876
    FD01848877
    FD01848883
    FD01848817
    FD01848819
    FD01848820
    FD01848821
    FD01848825
    FD01848826
    FD01848835
    FD01848836
    FD01848837
    FD01848841
    FD01848842
    FD01848845
    FD01848850
    FD01848857
    FD01848860
    FD01848861
    FD01848867
    FD01848873
    FD01848882""".split()
    tmp = genome[genome.pid.isin(curated.pid)].c80e3.unique()
    moreira = genome[(genome.assembly.isin(moreira)) & (genome.c80e3.isin(tmp))]
    moreira.to_csv("moreira.tsv", index=False, sep="\t")
    del(tmp)

# Report
_o = [ x for x in dir() if x not in _o and x[0] != "_" ]
print(f'The following variables were loaded: {_o}')

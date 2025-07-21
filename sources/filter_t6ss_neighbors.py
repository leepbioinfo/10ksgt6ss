#!/usr/bin/env python3

df = pd.read_pickle('../data/genome.pkl')
df.loc[df.rocha == 'T6SSi_tssH', 't6ss'] = False
df.loc[df.rocha == 'T6SSiii_tssH+T6SSi_tssH', 't6ss'] = False
df.loc[df.rocha == 'T6SSi_tssH+T6SSiii_tssH', 't6ss'] = False
t6 = df[df.t6ss ==True].t6ss
df2 = df.neighbors(t6, after=10, before=10, min_block_distance=4)
df2.to_pickle('vizinhos_t6_novo.pkl')

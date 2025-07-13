# coding: utf-8
import community
import networkx as nx
df = pd.read_pickle('./vizinhos_t6_novo.pkl')
viz = df.groupby('block_id').agg(vizi = ('c80e3', lambda x: set(pd.Series(x).dropna().drop_duplicates().sort_values().to_list())), vc = ('c80e3', 'nunique')).reset_index()
viz2 = viz.query('vc >= 3')
viz2 = viz2[~viz2.vizi.astype(str).duplicated(keep='first')]
viz2['t'] = 1
t2 = viz2.merge(viz2, on='t')
t2.block_id_x, t2.block_id_y = np.where(t2.block_id_x > t2.block_id_y , [t2.block_id_x, t2.block_id_y], [t2.block_id_y, t2.block_id_x])
t2 = t2.drop_duplicates(['block_id_x', 'block_id_y'])
print(f'tamanho t2 depois da ordenacao:  {len(t2)}')
t2 = t2[t2.block_id_x != t2.block_id_y]
t2['uni'] = t2.apply(lambda x: x.vizi_x.union(x.vizi_y), axis=1)
t2['inter'] = t2.apply(lambda x: x.vizi_x.intersection(x.vizi_y), axis=1)
t2['Jaccard_index'] = t2.apply(lambda x: len(x.inter)/len(x.uni), axis=1)
t2 = t2.sort_values('Jaccard_index')
t3 = t2[['block_id_x', 'block_id_y', 'Jaccard_index']].rename({'block_id_x': 'source', 'block_id_y':'target', 'Jaccard_index':'ji'}, axis=1)
G2 = nx.from_pandas_edgelist(t3[['source', 'target', 'ji']].query('ji >= 0.3'), edge_attr='ji')
partition = community.best_partition(G2,weight='ji')
c = pd.DataFrame.from_dict(partition,orient='index').reset_index().rename(
    {'index': 'block_id', 0: 'nei_c'}, axis=1)
viz.vizi = viz.vizi.astype('str')
v3 = viz.merge(c)
v3.vizi = v3.vizi.astype('str')
cc = [c  for c in sorted(nx.connected_components(G2), key=len, reverse=True)]
ncc = pd.Series(cc).reset_index().explode([0]).rename({'index':'ncc', 0: 'block_id'}, axis=1)
v4 = v3.merge(ncc, how='left')
v5 = viz.merge(v4[['vizi', 'nei_c', 'ncc']], how='left')
df.to_pickle('./10k_vizinho_novo_df_jaccard.pk')
t3.to_pickle('./jaccard_df_vizinho_novo.pk')

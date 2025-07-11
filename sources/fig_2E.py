from working_dfs import t2,li1, li3,li1i3, li2, type_dict, type_to_color, insidei1, insidei2, insidei3, final
import matplotlib.gridspec as gridspec
from venn import venn
from matplotlib.colors import  ListedColormap
import matplotlib.pyplot as plt
cmv = ListedColormap(['#74C3E3','#D184E1', '#FF958E', '#00C469' ])

si1 = set(final[['basename', 'i1', 'i2', 'i3', 'i4b']].set_index('basename').notna().eval('orphan  = (i1 + i2 + i3 + i4b) ==False').astype(int).reset_index().query('i1 ==1').basename)
si2 = set(final[['basename', 'i1', 'i2', 'i3', 'i4b']].set_index('basename').notna().eval('orphan  = (i1 + i2 + i3 + i4b) ==False').astype(int).reset_index().query('i2 ==1').basename)
si3 = set(final[['basename', 'i1', 'i2', 'i3', 'i4b']].set_index('basename').notna().eval('orphan  = (i1 + i2 + i3 + i4b) ==False').astype(int).reset_index().query('i3 ==1').basename)
si4b = set(final[['basename', 'i1', 'i2', 'i3', 'i4b']].set_index('basename').notna().eval('orphan  = (i1 + i2 + i3 + i4b) ==False').astype(int).reset_index().query('i4b ==1').basename)
sorphan = set(final[['basename', 'i1', 'i2', 'i3', 'i4b']].set_index('basename').notna().eval('orphan  = (i1 + i2 + i3 + i4b) ==False').astype(int).reset_index().query('orphan ==1').basename)
vend = {'i1':si1,'i2': si2, 'i3':si3,'orphan': sorphan}
venn(vend, cmap=cmv)
plt.savefig("./Figure_2E.png")

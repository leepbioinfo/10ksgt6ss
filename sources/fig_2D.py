from working_dfs import li1, li3, li1i3, t2
from functions import heat_map
l = li3+li1i3+li1
l.remove('FD01875449')
l.remove('FD01872508')
heat_map(t2,'genome','basename', orderby=l,top_columns=15, output='Figure_2D.pdf')

from functions import serovar_fig
from working_dfs import t2
import matplotlib
#set font size and family for all the plot
font = {'family' : 'Arial',
        'weight' : 'medium',
        'style'  : 'normal',
        'size'   : 5}

matplotlib.rc('font', **font)


serovar_fig(t2.serovar.value_counts().index.tolist(), 'Figure_S2.pdf')


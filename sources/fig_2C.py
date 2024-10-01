import matplotlib
from functions import serovar_fig

#set font size and family for all the plot
font = {'family' : 'Arial',
        'weight' : 'medium',
        'style'  : 'normal',
        'size'   : 5}

matplotlib.rc('font', **font)


to_fig_2C = ["S. Typhimurium", "S. Enteritidis", "S. Panama","S. Infantis","S. Agona","S. Dublin","S. Typhi"]



serovar_fig(to_fig_2C, 'Figure_2C.pdf')


import os

files = os.listdir('.')
if not 'data' in files: od.mkdir('./data')
if not 'plots' in files: od.mkdir('./plots')

import analysis.plot_afforestation

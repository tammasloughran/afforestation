# Main entry point of ananlysis package
import os
import sys

from analysis.plot_afforestation import make_plots

files = os.listdir('.')
if not 'data' in files: od.mkdir('./data')
if not 'plots' in files: od.mkdir('./plots')

sys.exit(make_plots())


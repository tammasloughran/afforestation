# Main entry point of ananlysis package
import os
import sys

from analysis.plot_afforestation import make_plots
from analysis.plot_land_cover_fractions import make_land_cover_plot

files = os.listdir('.')
if not 'data' in files: od.mkdir('./data')
if not 'plots' in files: od.mkdir('./plots')

sys.exit(make_plots())


# Main entry point of ananlysis package
import os
import sys

from analysis.plot_afforestation import make_veg_plots, make_clim_plots
from analysis.plot_land_cover_fractions import make_land_cover_plot

files = os.listdir('.')
if not 'data' in files: od.mkdir('./data')
if not 'plots' in files: od.mkdir('./plots')


def main():
    make_veg_plots()
    make_clim_plots()
    make_land_cover_plot()


sys.exit(main())


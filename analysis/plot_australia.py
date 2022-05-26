#!/usr/bin/env python3
import matplotlib.pyplot as plt

if __name__ != 'analysis.plot_australia':
    # plot_afforestation.py is main program or imported as a module from another script.
    from constants import ENSEMBLES, TABLES, VARIABLES, CLIM_VARIABLES
else:
    # plot_afforestation.py imported as a module of the analysis package.
    from analysis.constants import ENSEMBLES, TABLES, VARIABLES, CLIM_VARIABLES



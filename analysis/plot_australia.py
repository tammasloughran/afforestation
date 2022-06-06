#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

if __name__ != 'analysis.plot_australia':
    # plot_afforestation.py is main program or imported as a module from another script.
    from constants import ENSEMBLES, TABLES, VARIABLES, CLIM_VARIABLES
    from cmip_files import get_filename
    from cdo_calc_load import load_aus_pool
else:
    # plot_afforestation.py imported as a module of the analysis package.
    from analysis.constants import ENSEMBLES, TABLES, VARIABLES, CLIM_VARIABLES
    from analysis.cmip_files import get_filename
    from analysis.cdo_calc_load import load_aus_pool

NTIMES = 86
NENS = 10

var = np.ones((NENS,NTIMES))*np.nan
for i,ens in enumerate(ENSEMBLES):
    filename = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, 'Lmon', 'cVeg')[0]
    var[i,...] = load_aus_pool(input=filename, var='cVeg')
    plt.plot(var[i,...], color='grey')

var_mean = var.mean(axis=0)
plt.plot(var_mean, color='black')

plt.show()


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

years = list(range(2015, 2101))

data = np.ones((NENS,NTIMES))*np.nan
for i,ens in enumerate(ENSEMBLES):
    filename = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, 'Lmon', 'cVeg')[0]
    data[i,...] = load_aus_pool(input=filename, var='cVeg')
    plt.plot(years, data[i,...], color='grey')

data_mean = data.mean(axis=0)
data_std = data.std(axis=0, ddof=1)
plt.fill_between(years, data_mean+data_std, data_mean-data_std, color='grey', alpha=0.4)
plt.plot(years, data_mean, color='black')

plt.show()


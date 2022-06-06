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


def plot_region(years, data, data_mean, data_std):
    for i,ens in enumerate(ENSEMBLES):
        plt.plot(years, data[i,...], color='grey')
    plt.fill_between(years, data_mean+data_std, data_mean-data_std, color='grey', alpha=0.4)
    plt.plot(years, data_mean, color='black')
    plt.hlines(0, years[0], years[-1], linestyle='dotted', color='black')
    plt.ylabel('$\Delta$ Pg(C)')
    plt.xlabel('Time (Year)')
    plt.title("ACCES-ESM1.5 Australia cVeg")
    plt.show()


ref_file = '[ data/cVeg_ACCESS-ESM1-5_esm-hist-aff_ensmean_200501-202412mean.nc ]'
reference = load_aus_pool(input=ref_file, var='cVeg')

aff_data = np.ones((NENS,NTIMES))*np.nan
for i,ens in enumerate(ENSEMBLES):
    filenames = ' '.join(get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, 'Lmon', 'cVeg'))
    filenames = '[ '+filenames+' ]'
    aff_data[i,...] = load_aus_pool(input=filenames, var='cVeg')

anom_data = aff_data - reference

ssp585_data = np.ones((NENS,NTIMES))*np.nan
for i,ens in enumerate(ENSEMBLES):
    filenames = ' '.join(get_filename('C4MIP', 'esm-ssp585', ens, 'Lmon', 'cVeg'))
    filenames = '[ '+filenames+' ]'
    ssp585_data[i,...] = load_aus_pool(input=filenames, var='cVeg')

diff_data = ssp585_data - aff_data

data_mean = anom_data.mean(axis=0)
data_std = anom_data.std(axis=0, ddof=1)
plot_region(years, anom_data, data_mean, data_std)

data_mean = diff_data.mean(axis=0)
data_std = diff_data.std(axis=0, ddof=1)
plot_region(years, diff_data, data_mean, data_std)


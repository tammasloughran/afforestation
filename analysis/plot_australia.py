#!/usr/bin/env python3
import glob
import os

import matplotlib.pyplot as plt
import numpy as np

if __name__ != 'analysis.plot_australia':
    # plot_afforestation.py is main program or imported as a module from another script.
    from cdo_calc_load import load_aus_base_flux, load_aus_flux, load_aus_pool, load_aus_clim
    from cmip_files import get_filename
    from constants import CLIM_VARIABLES, ENSEMBLES, TABLES, VARIABLES, SEC_IN_DAY
else:
    # plot_afforestation.py imported as a module of the analysis package.
    from analysis.cdo_calc_load import (load_aus_base_flux, load_aus_flux,
                                        load_aus_pool, load_aus_clim)
    from analysis.cmip_files import get_filename
    from analysis.constants import CLIM_VARIABLES, ENSEMBLES, TABLES, VARIABLES, SEC_IN_DAY

NTIMES = 86
NENS = 10
COLORS = {'gpp':'green',
        'npp':'olive',
        'ra':'orange',
        'rh':'saddlebrown',
        'nbp':'purple',
        'cVeg':'darkgreen',
        'cLitter':'chocolate',
        'cSoil':'black',
        'tas':'black',
        'pr':'blue'}
PLOTS_DIR = './plots'
DATA_DIR = './data'

files = glob.glob('./*')
if PLOTS_DIR not in files:
    os.mkdir(PLOTS_DIR)
if DATA_DIR not in files:
    os.mkdir(DATA_DIR)

# Control flag
files = glob.glob('./data/*')
if any(['.npy' in f for f in files]):
    load_npy_files = True
else:
    load_npy_files = False
#load_npy_files = False # Uncomment to override previous check.

years = list(range(2015, 2101))


def plot_veg_region(years, data, data_mean, data_std, var, label=''):
    plt.figure()
    for i,ens in enumerate(ENSEMBLES):
        plt.plot(years, data[i,...], color='grey', alpha=0.4)
    plt.fill_between(years, data_mean+data_std, data_mean-data_std, color=COLORS[var], alpha=0.4)
    plt.plot(years, data_mean, color=COLORS[var])
    plt.hlines(0, years[0], years[-1], linestyle='dotted', color='black')
    plt.xlim(left=years[0], right=years[-1])
    plt.ylabel('$\Delta$ Pg(C)')
    plt.xlabel('Time (Year)')
    plt.title(f"ACCES-ESM1.5 Australia {var}")
    plt.savefig('plots/'+var+'_aus_'+label+'.png')


def plot_clim_region(years, data, data_mean, data_std, var, label=''):
    plt.figure()
    for i,ens in enumerate(ENSEMBLES):
        plt.plot(years, data[i,...], color='grey', alpha=0.4)
    plt.fill_between(years, data_mean+data_std, data_mean-data_std, color=COLORS[var], alpha=0.4)
    plt.plot(years, data_mean, color=COLORS[var])
    plt.xlim(left=years[0], right=years[-1])
    if var=='pr':
        plt.ylabel('mm/day')
    elif var=='tas':
        plt.ylabel('$^\circ$C')
    plt.xlabel('Time (Year)')
    plt.title(f"ACCES-ESM1.5 Australia {var}")
    plt.savefig('plots/'+var+'_aus_'+label+'.png')


def make_australia_plots():
    for table in TABLES:
        for var in VARIABLES[table]:
            print("Plotting region australia for ", var)
            # Get the reference period value.
            ref_file = '[ data/'+var+'_ACCESS-ESM1-5_esm-hist-aff_ensmean_200501-202412mean.nc ]'
            if var in ['cVeg', 'cLitter', 'cSoil']:
                reference = load_aus_pool(input=ref_file, var=var)
            elif var in ['gpp', 'npp', 'ra', 'rh', 'nbp']:
                reference = load_aus_base_flux(input=ref_file, var=var)

            # Load the afforestation data.
            if load_npy_files:
                aff_data = np.load(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585-ssp126Lu_aus.npy')
            else:
                aff_data = np.ones((NENS,NTIMES))*np.nan
                for i, ens in enumerate(ENSEMBLES):
                    filenames = ' '.join(
                            get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, table, var))
                    filenames = '[ '+filenames+' ]'
                    if var in ['cVeg', 'cLitter', 'cSoil']:
                        aff_data[i,:] = load_aus_pool(input=filenames, var=var)
                    elif var in ['gpp', 'npp', 'ra', 'rh', 'nbp']:
                        aff_data[i,:] = load_aus_flux(input=filenames, var=var)
                np.save(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585-ssp126Lu_aus.npy', aff_data)

            # Calculate anomaly.
            anom_data = aff_data - reference

            # Load the ssp585 data.
            if load_npy_files:
                ssp585_data = np.load(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585_aus.npy')
            else:
                ssp585_data = np.ones((NENS,NTIMES))*np.nan
                for i, ens in enumerate(ENSEMBLES):
                    filenames = ' '.join(get_filename('C4MIP', 'esm-ssp585', ens, table, var))
                    filenames = '[ '+filenames+' ]'
                    if var in ['cVeg', 'cLitter', 'cSoil']:
                        ssp585_data[i,:] = load_aus_pool(input=filenames, var=var)
                    elif var in ['gpp', 'npp', 'ra', 'rh', 'nbp']:
                        ssp585_data[i,:] = load_aus_flux(input=filenames, var=var)
                np.save(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585_aus.npy', ssp585_data)

            # Clalculate the difference.
            diff_data = aff_data - ssp585_data

            # Plot
            anom_data_mean = anom_data.mean(axis=0)
            anom_data_std = anom_data.std(axis=0, ddof=1)
            plot_veg_region(years, anom_data, anom_data_mean, anom_data_std, var, label='anom')

            diff_data_mean = diff_data.mean(axis=0)
            diff_data_std = diff_data.std(axis=0, ddof=1)
            plot_veg_region(years, diff_data, diff_data_mean, diff_data_std, var, label='diff')

    for var in ['tas', 'pr']:
        print("Plotting region australia for ", var)
        # Load the afforestation data
        if load_npy_files:
            aff_data = np.load(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585-ssp126Lu_aus.npy')
            ssp585_data = np.load(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585_aus.npy')
        else:
            aff_data = np.ones((NENS,NTIMES))*np.nan
            ssp585_data = np.ones((NENS,NTIMES))*np.nan
            for i, ens in enumerate(ENSEMBLES):
                filenames = ' '.join(get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, 'Amon', var))
                filenames = '[ '+filenames+' ]'
                aff_data[i,:] = load_aus_clim(input=filenames, var=var)
                filenames = ' '.join(get_filename('C4MIP', 'esm-ssp585', ens, 'Amon', var))
                filenames = '[ '+filenames+' ]'
                ssp585_data[i,:] = load_aus_clim(input=filenames, var=var)
            np.save(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585-ssp126Lu_aus.npy', aff_data)
            np.save(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585_aus.npy', ssp585_data)

        if var=='tas':
            aff_data -= 273.15
            ssp585_data -= 273.15
        elif var=='pr':
            aff_data *= SEC_IN_DAY
            ssp585_data *= SEC_IN_DAY

        # Calculate the difference
        diff_data = aff_data - ssp585_data

        # Plot
        plot_clim_region(years, aff_data, aff_data.mean(axis=0), aff_data.std(axis=0, ddof=1),
                var, label='aff_'+var)
        plot_clim_region(years, ssp585_data, ssp585_data.mean(axis=0),
                ssp585_data.std(axis=0, ddof=1), var, label='ssp585_'+var)
        plot_clim_region(years, diff_data, diff_data.mean(axis=0), diff_data.std(axis=0, ddof=1),
                var, label='diff_'+var)

    plt.show()

    # Clean up
    temp_files = glob.glob('./cdoPy*')
    for f in temp_files:
        os.remove(f)

if __name__ != 'analysis.plot_australia':
    make_australia_plots()


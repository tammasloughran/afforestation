#!/usr/bin/env python3
# plot_afforestation.py generates plots for the afforestation simulation.
# In general "aff" is used to refer to the afforestation simulation esm-ssp585-ssp126Lu (high
# fossil fuels emissions and low land use change scenario), and "ssp585" means the emissions driven
# esm-ssp585 simulation (high fossil fuel emisssions scenario).
import glob
import os
import pdb
import sys

import cdo as cdo_module
import matplotlib.pyplot as plt
import numpy as np

if __name__ == 'analysis.plot_afforestation':
    # plot_afforestation.py imported as a module of the analysis package.
    from analysis.baseline import global_sum_baselines
    from analysis.cdo_calc_load import cdo_fetch_ensembles
    from analysis.constants import (
            CLIM_VARIABLES,
            DATA_DIR,
            ENSEMBLES,
            PLOTS_DIR,
            SEC_IN_DAY,
            TABLES,
            VARIABLES,
            )
else:
    # plot_afforestation.py is main program or imported as a module from another script.
    from baseline import global_sum_baselines
    from cdo_calc_load import cdo_fetch_ensembles
    from constants import (
            CLIM_VARIABLES,
            DATA_DIR,
            ENSEMBLES,
            PLOTS_DIR,
            SEC_IN_DAY,
            TABLES,
            VARIABLES,
            )

# Local constants
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

files = glob.glob('./*')
if PLOTS_DIR not in files:
    os.mkdir(PLOTS_DIR)
if DATA_DIR not in files:
    os.mkdir(DATA_DIR)

# Control flag
files = glob.glob(f'{DATA_DIR}/cVeg_*_global.npy')
if any(['.npy' in f for f in files]):
    load_npy_files = True
else:
    load_npy_files = False
#load_npy_files = True # Uncomment to override previous check.


def plot_ensembles(years, data, data_mean, data_std, var):
    """Plot all ensemble members with ensemble mean and standard deviation.
    """
    model = 'ACCESS-ESM1.5'
    plt.figure()
    for ens in ENSEMBLES:
        plt.plot(years, data[int(ens)-1], color='lightgray', linewidth=0.6, alpha=0.4)
    plt.plot(years, data_mean, color=COLORS[var], label="Ensemble mean")
    plt.fill_between(years, data_mean+data_std, data_mean-data_std, color=COLORS[var],
            label="+-1$\sigma$", alpha=0.4)
    plt.hlines(0, years[0], years[-1], colors='black', linestyles='dotted')
    plt.xlim(left=years[0], right=years[-1])
    plt.xlabel('Year')
    plt.ylabel(f'{var.upper()} anomaly (PgC/year)')
    plt.title(f"{model} {var.upper()}")


def plot_ensembles_clim(years, data, data_mean, data_std, var):
    """Plot all ensemble members with ensemble mean and standard deviation.
    """
    model = 'ACCESS-ESM1.5'
    plt.figure()
    for ens in ENSEMBLES:
        plt.plot(years, data[int(ens)-1], color='lightgray', linewidth=0.6, alpha=0.4)
    plt.plot(years, data_mean, color=COLORS[var], label="Ensemble mean")
    plt.fill_between(years, data_mean+data_std, data_mean-data_std, color=COLORS[var],
            label="+-1$\sigma$", alpha=0.4)
    plt.xlim(left=years[0], right=years[-1])
    plt.xlabel('Year')
    if var=='tas':
        plt.ylabel(var.upper()+' ($^{\circ}$C)')
    elif var=='pr':
        plt.ylabel(f'{var.upper()} (mm/day)')
    plt.title(f"{model} {var.upper()}")


def make_veg_plots():
    model = 'ACCESS-ESM1.5'
    for table in TABLES:
        for var in VARIABLES[table]:
            print(f"Processing {var}")
            # Load data
            if load_npy_files:
                aff_data = np.load(f'{DATA_DIR}/{var}_{model}_esm-ssp585-ssp126Lu_global.npy')
                ssp585_data = np.load(f'{DATA_DIR}/{var}_{model}_esm-ssp585_global.npy')
            else:
                aff_data = cdo_fetch_ensembles('LUMIP', 'esm-ssp585-ssp126Lu', table, var)
                ssp585_data = cdo_fetch_ensembles('C4MIP', 'esm-ssp585', table, var)
                np.save(f'{DATA_DIR}/{var}_{model}_esm-ssp585-ssp126Lu_global.npy', aff_data)
                np.save(f'{DATA_DIR}/{var}_{}model_esm-ssp585_global.npy', ssp585_data)


            # Anomaly, mean and standard deviation relative to 2015 baseline. Demonstrates overall
            # impact of climate, CO2 forcing and afforestation on pool/flux.
            data_anomaly = aff_data - global_sum_baselines[var]
            data_ens_mean = np.mean(data_anomaly, axis=0)
            data_ens_std = np.std(data_anomaly, axis=0, ddof=1)

            # Anomaly relative to the esm-ssp585 scenario from C4MIP. Demonstrates only the changes
            # due to afforestation.
            data_aff_diff = aff_data - ssp585_data
            data_aff_diff_mean = np.mean(data_aff_diff, axis=0)
            data_aff_diff_std = np.std(data_aff_diff, axis=0, ddof=1)

            # Plot the graphs for anomalies relative to 2015.
            years = list(range(2015, 2101))
            plot_ensembles(years, data_anomaly, data_ens_mean, data_ens_std, var)
            plt.savefig(f'{PLOTS_DIR}/'+ \
                    f'{var}_{model}_esm-ssp585-ssp126Lu_ensembles_anomalies.png')
            plt.close()

            # Plot the graphs for anomalies relative to the esm-ssp585
            plot_ensembles(years, data_aff_diff, data_aff_diff_mean, data_aff_diff_std, var)
            plt.savefig(f'{PLOTS_DIR}/'+ \
                    f'{var}_{model}_esm-ssp585-ssp126Lu_ensembles_diff.png')
            plt.close()


def make_clim_plots():
    model = 'ACCESS-ESM1.5'
    table = 'Amon'
    for var in CLIM_VARIABLES[table]:
        print(f"Processing {var}")
        if load_npy_files:
            aff_data = np.load(f'{DATA_DIR}/{var}_{model}_esm-ssp585-ssp126Lu_global.npy')
            ssp585_data = np.load(f'{DATA_DIR}/{var}_{model}_esm-ssp585_global.npy')
        else:
            aff_data = cdo_fetch_ensembles('LUMIP', 'esm-ssp585-ssp126Lu', table, var)
            ssp585_data = cdo_fetch_ensembles('C4MIP', 'esm-ssp585', table, var)
            np.save(f'{DATA_DIR}/{var}_{model}_esm-ssp585-ssp126Lu_global.npy', aff_data)
            np.save(f'{DATA_DIR}/{var}_{model}_esm-ssp585_global.npy', ssp585_data)

        # Calculate mean and standar deviation
        if var == 'tas':
            aff_data -= 273.15
            ssp585_data -= 273.15
        elif var=='pr':
            aff_data *= SEC_IN_DAY
            ssp585_data *= SEC_IN_DAY
        aff_mean = np.mean(aff_data, axis=0)
        ssp585_mean = np.mean(ssp585_data, axis=0)
        aff_std = np.std(aff_data, axis=0, ddof=1)
        ssp585_std = np.std(ssp585_data, axis=0, ddof=1)

        # Plot
        years = list(range(2015, 2101))
        plot_ensembles_clim(years, aff_data, aff_mean, aff_std, var)
        plt.savefig(f'{PLOTS_DIR}/'+\
                f'{var}_{model}_esm-ssp585-ssp126Lu_ensembles.png')
        plt.close()
        plot_ensembles_clim(years, ssp585_data, ssp585_mean, ssp585_std, var)
        plt.savefig(f'{PLOTS_DIR}/'+\
                f'{var}_{model}_esm-ssp585_ensembles.png')
        plt.close()
        plot_ensembles_clim(years, ssp585_data-aff_data, ssp585_mean-aff_mean,
                aff_std, var)
        plt.savefig(f'{PLOTS_DIR}/'+\
                f'{var}_{model}_esm-ssp585_ensembles_diff.png')
        plt.close()


if __name__ != 'analysis.plot_afforestation':
    make_veg_plots()
    make_clim_plots()

    # Clean up
    temp_files = glob.glob('./cdoPy*')
    for f in temp_files:
        os.remove(f)

    plt.show()


#!/usr/bin/env python3
# plot_afforestation.py generates plots for the afforestation simulation.
# In general "aff" is used to refer to the afforestation simulation esm-ssp585-ssp126Lu (high
# fossil fuels emissions and low land use change scenario), and "ssp585" means the emissions driven
# esm-ssp585 simulation (high fossil fuel emisssions scenario).
import glob
import os
import sys

import cdo as cdo_module
import matplotlib.pyplot as plt
import numpy as np

if __name__ != 'analysis.plot_afforestation':
    # plot_afforestation.py is main program or imported as a module from another script.
    from baseline import global_sum_baselines
    from cdo_calc_load import cdo_fetch_ensembles
    from constants import ENSEMBLES, TABLES, VARIABLES, CLIM_VARIABLES
else:
    # plot_afforestation.py imported as a module of the analysis package.
    from analysis.baseline import global_sum_baselines
    from analysis.cdo_calc_load import cdo_fetch_ensembles
    from analysis.constants import ENSEMBLES, TABLES, VARIABLES, CLIM_VARIABLES

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
PLOTS_DIR = 'plots'


def plot_ensembles(years, data, data_mean, data_std, var):
    """Plot all ensemble members with ensemble mean and standard deviation.
    """
    plt.figure()
    for ens in ENSEMBLES:
        plt.plot(years, data[int(ens)-1], color='lightgray', linewidth=0.6)
    plt.plot(years, data_mean, color=COLORS[var], label="Ensemble mean")
    plt.plot(years, data_mean+data_std, color=COLORS[var], linewidth=0.8,
            label="+-1$\sigma$")
    plt.plot(years, data_mean-data_std, color=COLORS[var], linewidth=0.8)
    plt.hlines(0, years[0], years[-1], colors='black', linestyles='dotted')
    plt.xlabel('Year')
    plt.ylabel(f'{var.upper()} anomaly (PgC/year)')
    plt.title(f"ACCESS-ESM1-5 {var.upper()}")


def plot_ensembles_clim(years, data, data_mean, data_std, var):
    """Plot all ensemble members with ensemble mean and standard deviation.
    """
    plt.figure()
    for ens in ENSEMBLES:
        plt.plot(years, data[int(ens)-1], color='lightgray', linewidth=0.6)
    plt.plot(years, data_mean, color=COLORS[var], label="Ensemble mean")
    plt.plot(years, data_mean+data_std, color=COLORS[var], linewidth=0.8,
            label="+-1$\sigma$")
    plt.plot(years, data_mean-data_std, color=COLORS[var], linewidth=0.8)
    plt.xlabel('Year')
    plt.ylabel(f'{var.upper()}')
    plt.title(f"ACCESS-ESM1-5 {var.upper()}")


def make_veg_plots():
    for table in TABLES:
        for var in VARIABLES[table]:
            print(f"Processing {var}")
            # Load data
            aff_data = cdo_fetch_ensembles('LUMIP', 'esm-ssp585-ssp126Lu', table, var)
            ssp585_data = cdo_fetch_ensembles('C4MIP', 'esm-ssp585', table, var)

            # Anomaly, mean and standard deviation relative to 2015 baseline. Demonstrates overall
            # impact of climate, CO2 forcing and afforestation on pool/flux.
            data_anomaly = aff_data - global_sum_baselines[var]
            data_ens_mean = np.mean(data_anomaly, axis=0)
            data_ens_std = np.std(data_anomaly, axis=0, ddof=1)

            # Anomaly relative to the esm-ssp585 scenario from C4MIP. Domenstrates only the changes
            # due to afforestation.
            data_aff_diff = aff_data - ssp585_data
            data_aff_diff_mean = np.mean(data_aff_diff, axis=0)
            data_aff_diff_std = np.std(data_aff_diff, axis=0, ddof=1)

            # Plot the graphs for anomalies relative to 2015.
            years = [y for y in range(2015, 2101)]
            plot_ensembles(years, data_anomaly, data_ens_mean, data_ens_std, var)
            plt.savefig(f'{PLOTS_DIR}/'+ \
                    f'{var}_ACCESS-ESM1-5_esm-ssp585-ssp126Lu_ensembles_anomalies.svg')

            # Plot the graphs for anomalies relative to the esm-ssp585
            plot_ensembles(years, data_aff_diff, data_aff_diff_mean, data_aff_diff_std, var)
            plt.savefig(f'{PLOTS_DIR}/'+ \
                    f'{var}_ACCESS-ESM1-5_esm-ssp585-ssp126Lu_ensembles_diff.svg')
            plt.close()


def make_clim_plots():
    table = 'Amon'
    for var in CLIM_VARIABLES[table]:
        print(f"Processing {var}")
        aff_data = cdo_fetch_ensembles('LUMIP', 'esm-ssp585-ssp126Lu', table, var)
        ssp585_data = cdo_fetch_ensembles('C4MIP', 'esm-ssp585', table, var)

        # Calculate mean and standar deviation
        if var == 'tas':
            aff_data -= 273.15
            ssp585_data -= 273.15
        aff_mean = np.mean(aff_data, axis=0)
        ssp585_mean = np.mean(ssp585_data, axis=0)
        aff_std = np.std(aff_data, axis=0, ddof=1)
        ssp585_std = np.std(ssp585_data, axis=0, ddof=1)

        # Plot
        years = [y for y in range(2015, 2101)]
        plot_ensembles_clim(years, aff_data, aff_mean, aff_std, var)
        plt.savefig(f'{PLOTS_DIR}/'+\
                f'{var}_ACCESS-ESM1-5_esm-ssp585-ssp126Lu_ensembles.svg')
        plot_ensembles_clim(years, ssp585_data, ssp585_mean, ssp585_std, var)
        plt.savefig(f'{PLOTS_DIR}/'+\
                f'{var}_ACCESS-ESM1-5_esm-ssp585_ensembles.svg')
        plot_ensembles_clim(years, ssp585_data-aff_data, ssp585_mean-aff_mean,
                aff_std, var)
        plt.savefig(f'{PLOTS_DIR}/'+\
                f'{var}_ACCESS-ESM1-5_esm-ssp585_ensembles_diff.svg')
        plt.close()


if __name__ != 'analysis.plot_afforestation':
    make_veg_plots()
    make_clim_plots()

    # Clean up
    temp_files = glob.glob('./cdoPy*')
    for f in temp_files:
        os.remove(f)

    #plt.show()


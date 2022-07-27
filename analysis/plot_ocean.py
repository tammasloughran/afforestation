#!/usr/bin/env python3
"""Plot the ocean carbon in the afforestation scenario and the reference ssp585 scenario.
"""
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from cdo import Cdo
import cdo_decorators as cdod
if __name__ == 'analysis.plot_afforestation':
    # plot_afforestation.py imported as a module of the analysis package.
    from analysis.cmip_files import get_filename
    from analysis.constants import (
            DATA_DIR,
            ENSEMBLES,
            PLOTS_DIR,
            )
else:
    # plot_afforestation.py is main program or imported as a module from another script.
    from cmip_files import get_filename
    from constants import (
            DATA_DIR,
            ENSEMBLES,
            PLOTS_DIR,
            )
import ipdb

cdo = Cdo()
cdo.debug = True

OCEAN_VARIABLES = {
        'Omon':[
                'fgco2', # Surface Downward Flux of Total CO2. Units kg m-2 s-1
                #'fgco2nat', # Surface Downward Flux of Natural CO2. Units kg m-2 s-1
                #'intpp', # Primary Organic Carbon Production by All Types of Phytoplankton.
                          # Units mol m-2 s-1
                #'spco2', # Surface Aqueous Partial Pressure of CO2. Units Pa
                #'spco2nat', # Surface Aqueous Partial Pressure of Natural CO2. Units Pa
                ],
        }
NENS = len(ENSEMBLES)
NTIMES = 2101 - 2015

load_cdo = False
load_npy = not load_cdo


@cdod.cdo_cat(input2='') # Concatenate all files in input1
@cdod.cdo_fldsum
@cdod.cdo_yearmonmean
def cdo_load_ocean(var:str, input:str)->np.ma.MaskedArray:
    """Use CDO to load an ocean carbon variable.
    """
    return cdo.copy(input=input, returnCdf=True, options='-L').variables[var][:].squeeze()


def plot_ocean_carbon(aff_data, ssp585_data, label:str='')->None:
    years = list(range(2015, 2015+NTIMES))
    plt.figure()
    diff = aff_data - ssp585_data
    for e,ens in enumerate(ENSEMBLES):
        plt.plot(years, diff[e,:], color='lightgray')
    plt.plot(years, diff.mean(axis=0), color='black')
    plt.xlabel('Year')
    plt.ylabel('C')
    plt.savefig(f'{PLOTS_DIR}/ocean_carbon_aff_ssp585_{label}.png')


def make_ocean_carbon_plot():
    # Load Data
    table = 'Omon'
    for var in OCEAN_VARIABLES[table]:
        if load_cdo:
            aff_data = np.ones((NENS,NTIMES))*np.nan
            ssp585_data = np.ones((NENS,NTIMES))*np.nan
            for e,ens in enumerate(ENSEMBLES):
                aff_files = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, table, var)
                aff_files = '[ '+' '.join(aff_files)+' ]'
                aff_data[e,:] = cdo_load_ocean(var=var, input=aff_files)
                ssp585_files = get_filename('C4MIP', 'esm-ssp585', ens, table, var)
                ssp585_files = '[ '+' '.join(ssp585_files)+' ]'
                ssp585_data[e,:] = cdo_load_ocean(var=var, input=ssp585_files)
            # Save.
            np.save(f'{DATA_DIR}/{var}_aff_ocean_carbon_data.npy', aff_data)
            np.save(f'{DATA_DIR}/{var}_ssp585_ocean_carbon_data.npy', ssp585_data)
        elif load_npy:
            aff_data = np.load(f'{DATA_DIR}/{var}_aff_ocean_carbon_data.npy')
            ssp585_data = np.load(f'{DATA_DIR}/{var}_ssp585_ocean_carbon_data.npy')

        # Plot
        plot_ocean_carbon(aff_data, ssp585_data, label='flux')
        plot_ocean_carbon(
                np.cumsum(aff_data, axis=1),
                np.cumsum(ssp585_data, axis=1),
                label='cumulative',
                )


if __name__ != 'analysis.plot_ocean':
    make_ocean_carbon_plot()

    # Clean up temp files.
    temp_files = glob.glob('./cdoPy*')
    for f in temp_files:
        os.remove(f)

    plt.show()


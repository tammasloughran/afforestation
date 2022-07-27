#!/usr/bin/env python3
"""Plot the atmosphere co2 concentrations.
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
cdo.debug = False

CO2_VARIABLES = {
        'Emon':'co23D', # CO2 concentration (mixing ratio) in kg/kg.
        #'Amon':'co2', # Mole fraction of CO2. esm-ssp585 only
        }
MIL = 1000000
MOLMASS_CO2 = 0.0440095 # kg/mol
MOLMASS_AIR = 0.0289652 # kg/mol
KGKG_TO_MOLMOL = MOLMASS_AIR/MOLMASS_CO2 # Converion of kg/kg to mol/mol
NENS = len(ENSEMBLES)
NTIMES = 2101 - 2015

load_cdo = False
load_npy = not load_cdo


@cdod.cdo_vertmean
@cdod.cdo_fldmean
@cdod.cdo_yearmonmean
def global_average(input:str, output:str)->None:
    """Create intermediate files for field and vertical aggregated yearly data using CDO.
    """
    cdo.copy(input=input, output=output, options='-L')


@cdod.cdo_cat(input2='')
def load_co2(var:str, input:str):
    """Concatenate and load the data using CDO.
    """
    return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


def make_co2_plot():
    """Load the CO2 data either from netcdf files with CDO or numpy files and plot.
    """
    fig = plt.figure()
    aff_co2_data = np.ones((NENS,NTIMES))*np.nan
    ssp585_co2_data = np.ones((NENS,NTIMES))*np.nan

    # Load data.
    if load_cdo:
        table = 'Emon'
        var = CO2_VARIABLES[table]
        for e,ens in enumerate(ENSEMBLES):
            aff_co2_files = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, table, var)
            for f in aff_co2_files:
                global_average(input=f, output=f.split('/')[-1]+'_global.nc')
            aff_co2_files = '[ '+' '.join(sorted(glob.glob('*_global.nc')))+' ]'
            aff_co2_data[e,:] = load_co2(var=var, input=aff_co2_files)
            junk_files = glob.glob('*_global.nc')
            for f in junk_files: os.remove(f)

            ssp585_co2_files = get_filename('C4MIP', 'esm-ssp585', ens, table, var)
            for f in ssp585_co2_files:
                global_average(input=f, output=f.split('/')[-1]+'_global.nc')
            ssp585_co2_files = '[ '+' '.join(sorted(glob.glob('*_global.nc')))+' ]'
            ssp585_co2_data[e,:] = load_co2(var=var, input=ssp585_co2_files)
            junk_files = glob.glob('*_global.nc')
            for f in junk_files: os.remove(f)
        np.save(f'{DATA_DIR}/aff_co2_data.npy', aff_co2_data)
        np.save(f'{DATA_DIR}/ssp585_co2_data.npy', ssp585_co2_data)
    elif load_npy:
        aff_co2_data = np.load(f'{DATA_DIR}/aff_co2_data.npy')
        ssp585_co2_data = np.load(f'{DATA_DIR}/ssp585_co2_data.npy')

    # Plot each ensemble member.
    for e,ens in enumerate(ENSEMBLES):
        years = list(range(2015, 2015+NTIMES))
        plt.plot(
                years,
                aff_co2_data[e,:]*KGKG_TO_MOLMOL*MIL,
                color='blue',
                linewidth=0.6,
                alpha=0.4,
                )
        plt.plot(
                years,
                ssp585_co2_data[e,:]*KGKG_TO_MOLMOL*MIL,
                color='red',
                linewidth=0.6,
                alpha=0.4,
                )
    # Plot ensemble mean.
    plt.plot(
            years,
            aff_co2_data.mean(axis=0)*KGKG_TO_MOLMOL*MIL,
            color='blue',
            label='esm-ssp585-ssp126Lu',
            )
    plt.plot(
            years,
            ssp585_co2_data.mean(axis=0)*KGKG_TO_MOLMOL*MIL,
            color='red',
            label='esm-ssp585',
            )
    plt.legend()
    plt.xlabel('Year')
    plt.ylabel('CO$_2$ concentration (ppm)')
    plt.savefig(f'{PLOTS_DIR}/CO2_aff_and_ssp585.png')


if __name__ != 'analysis.plot_atmosphere':
    make_co2_plot()

    # Clean up temp files.
    temp_files = glob.glob('./cdoPy*')
    for f in temp_files:
        os.remove(f)

    plt.show()


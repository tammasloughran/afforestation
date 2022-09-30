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
            C_IN_CO2_RATIO,
            DATA_DIR,
            ENSEMBLES,
            KG_IN_PG,
            MOLMASS_O2, # kg/mol
            NENS,
            NTIMES,
            PLOTS_DIR,
            SEC_IN_DAY,
            )
else:
    # plot_afforestation.py is main program or imported as a module from another script.
    from cmip_files import get_filename
    from constants import (
            C_IN_CO2_RATIO,
            DATA_DIR,
            ENSEMBLES,
            KG_IN_PG,
            MOLMASS_O2, # kg/mol
            NENS,
            NTIMES,
            PLOTS_DIR,
            SEC_IN_DAY,
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
EXAMPLE_FILE = '/g/data/fs38/publications/CMIP6/C4MIP/CSIRO/ACCESS-ESM1-5/esm-ssp585/r1i1p1f1/'+\
        'Omon/fgco2/gn/latest/fgco2_Omon_ACCESS-ESM1-5_esm-ssp585_r1i1p1f1_gn_201501-210012.nc'

load_cdo = False
load_npy = not load_cdo

# If ocean grid area file does not exist, create it.
if 'ocean_gridarea.nc' not in os.listdir(DATA_DIR):
    cdo.gridarea(
            input=EXAMPLE_FILE,
            output=f'{DATA_DIR}/ocean_gridarea.nc',
            options='-L', # Lock I/O (sequential access).
            )


@cdod.cdo_cat(input2='') # Concatenate all files in input1
@cdod.cdo_mul(input2=f'{DATA_DIR}/ocean_gridarea.nc') # Convert from /m to per gridcell.
@cdod.cdo_fldsum # Spatial aggregation
@cdod.cdo_mulc(str(SEC_IN_DAY)) # Convertfrom /s into /day.
@cdod.cdo_muldpm # Convert from /day into /month.
@cdod.cdo_yearsum # Convert from /month into /year
@cdod.cdo_divc(str(KG_IN_PG)) # Convert from kg into Pg.
def cdo_load_ocean_flux(var:str, input:str)->np.ma.MaskedArray:
    """Use CDO to load an ocean carbon variable.
    """
    return cdo.copy(input=input, returnCdf=True, options='-L').variables[var][:].squeeze()


def plot_ocean_carbon(aff_data:np.ndarray, ssp585_data:np.ndarray, var, label:str='flux')->None:
    """Create a plot of ocean carbon for all ensembles.
    Input data is in CO2 but plot is converted to C 
    """
    years = list(range(2015, 2015+NTIMES))
    plt.figure()
    diff = (aff_data - ssp585_data)*C_IN_CO2_RATIO
    for e,ens in enumerate(ENSEMBLES):
        plt.plot(years, diff[e,:], color='lightblue', alpha=0.5)
    plt.plot(years, diff.mean(axis=0), color='darkblue', label=f'{var} Aff. - SSP585')
    plt.hlines(0, years[0], years[-1], color='black', linestyle='dashed')
    plt.xlabel('Year')
    if label=='flux':
        plt.title('ACCESS-ESM1.5 Difference Ocean Downward Carbon flux')
        plt.ylabel('$\Delta$C [Pg(C)/year]')
    else:
        plt.title('ACCESS-ESM1.5 Difference Cumulative Ocean Carbon')
        plt.ylabel('$\Delta$C [Pg(C)]')
    plt.legend()
    plt.savefig(f'{PLOTS_DIR}/{var}_ocean_carbon_aff_ssp585_{label}.png')


def make_ocean_carbon_plot()->None:
    """Function to load ocean carbon data and execute plotting routine.
    """
    # Load Data
    table = 'Omon'
    for var in OCEAN_VARIABLES[table]:
        if load_cdo:
            aff_data = np.ones((NENS,NTIMES))*np.nan
            ssp585_data = np.ones((NENS,NTIMES))*np.nan
            for e,ens in enumerate(ENSEMBLES):
                aff_files = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, table, var)
                aff_files = '[ '+' '.join(aff_files)+' ]'
                aff_data[e,:] = cdo_load_ocean_flux(var=var, input=aff_files)
                ssp585_files = get_filename('C4MIP', 'esm-ssp585', ens, table, var)
                ssp585_files = '[ '+' '.join(ssp585_files)+' ]'
                ssp585_data[e,:] = cdo_load_ocean_flux(var=var, input=ssp585_files)
            # Save.
            np.save(f'{DATA_DIR}/{var}_aff_ocean_carbon_data.npy', aff_data)
            np.save(f'{DATA_DIR}/{var}_ssp585_ocean_carbon_data.npy', ssp585_data)
        elif load_npy:
            aff_data = np.load(f'{DATA_DIR}/{var}_aff_ocean_carbon_data.npy')
            ssp585_data = np.load(f'{DATA_DIR}/{var}_ssp585_ocean_carbon_data.npy')

        # Plot
        plot_ocean_carbon(aff_data, ssp585_data, var, label='flux')
        plot_ocean_carbon(
                np.cumsum(aff_data, axis=1),
                np.cumsum(ssp585_data, axis=1),
                var,
                label='cumulative',
                )


if __name__ != 'analysis.plot_ocean':
    make_ocean_carbon_plot()

    # Clean up temp files.
    temp_files = glob.glob('./cdoPy*')
    for f in temp_files:
        os.remove(f)

    plt.show()


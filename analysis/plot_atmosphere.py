#!/usr/bin/env python3
"""Plot the atmosphere co2 concentrations.
"""
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
from cdo import Cdo
import cdo_decorators as cdod
if __name__ == 'analysis.plot_afforestation':
    # plot_afforestation.py imported as a module of the analysis package.
    from analysis.cmip_files import get_filename
    from analysis.constants import (
            DATA_DIR,
            ENSEMBLES,
            KG_IN_PG,
            MASS_ATM, # kg
            PLOTS_DIR,
            MIL,
            MOLMASS_CO2, # kg/mol
            MOLMASS_O2, # kg/mol
            MOLMASS_AIR, # kg/mol
            KGKG_TO_MOLMOL, # Converion of kg/kg to mol/mol
            NENS,
            NTIMES,
            )
else:
    # plot_afforestation.py is main program or imported as a module from another script.
    from cmip_files import get_filename
    from constants import (
            DATA_DIR,
            ENSEMBLES,
            KG_IN_PG,
            MASS_ATM, # kg
            PLOTS_DIR,
            MIL,
            MOLMASS_CO2, # kg/mol
            MOLMASS_O2, # kg/mol
            MOLMASS_AIR, # kg/mol
            KGKG_TO_MOLMOL, # Converion of kg/kg to mol/mol
            NENS,
            NTIMES,
            )
import ipdb

cdo = Cdo()
cdo.debug = False

CO2_VARIABLES = {
        'Emon':'co23D', # CO2 concentration (mixing ratio) in kg/kg.
        #'Amon':'co2', # Mole fraction of CO2. esm-ssp585 only
        }

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
def load_co2(var:str, input:str)->np.ma.MaskedArray:
    """Concatenate and load the data using CDO.
    """
    return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


@cdod.cdo_mulc()
def load_mass(var:str, input:str)->None:
    """Load atmospheric CO2 as total mass in Pg(CO$_2$).
    """
    return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


def make_co2_plot()->None:
    """Load the CO2 data either from netcdf files with CDO or numpy files and plot.
    """
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
    fig = plt.figure()
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

    # Plot mass difference.
    #fig = plt.figure()
    #for e,ens in enumerate(ENSEMBLES):
    #    years = list(range(2015, 2015+NTIMES))
    #    plt.plot(
    #            years,
    #            (aff_co2_data[e,:] - ssp585_co2_data[e,:])*(MASS_ATM/KG_IN_PG)/MOLMASS_O2,
    #            color='royalblue',
    #            linewidth=0.6,
    #            alpha=0.4,
    #            )
    #plt.plot(
    #    years,
    #    (aff_co2_data - ssp585_co2_data).mean(axis=0)*(MASS_ATM/KG_IN_PG)/MOLMASS_O2,
    #    color='royalblue',
    #    label='aff - ssp585',
    #    )
    #plt.hlines(0, years[0], years[-1], linestyle='dashed', color='black')
    #plt.legend()
    #plt.xlabel('Year')
    #plt.ylabel('Carbon mass [Pg(C)]')
    #plt.savefig(f'{PLOTS_DIR}/CO2_atm_mass_diff.png')

    # Carbon budget
    # Ocean
    for_ocean = np.load(f'{DATA_DIR}/fgco2_aff_ocean_carbon_data.npy')*0.27291
    ssp585_ocean = np.load(f'{DATA_DIR}/fgco2_ssp585_ocean_carbon_data.npy')*0.27291
    # Land
    for_cland = np.load(f'{DATA_DIR}/cLand_ACCESS-ESM1.5_esm-ssp585-ssp126Lu_global.npy')
    for_cland_flux = np.concatenate([np.repeat(np.nan, 10)[:,None], np.diff(for_cland, 1)], axis=1)
    ssp585_cland = np.load(f'{DATA_DIR}/cLand_ACCESS-ESM1.5_esm-ssp585_global.npy')
    ssp585_cland_flux = np.concatenate([np.repeat(np.nan, 10)[:,None], np.diff(ssp585_cland, 1)], axis=1)
    # Fossil fuel emissions
    emissions_file = nc.Dataset(f'{DATA_DIR}/ssp585_emissions_global.nc', 'r')
    ffemissions = emissions_file.variables['field1561'][1:].squeeze()*0.27291
            # 0.27291 %mass of C in co2 to convert from co2 to C

    for_cland_flux_mean = for_cland_flux.mean(axis=0)
    for_cland_flux_mean[0] = 0
    ssp585_cland_flux_mean = ssp585_cland_flux.mean(axis=0)
    for_ocean_mean = for_ocean.mean(axis=0)
    ssp585_ocean_mean = ssp585_ocean.mean(axis=0)
    for_atmosphere = ffemissions - for_cland_flux_mean - for_ocean_mean
    ssp585_atmosphere = ffemissions - ssp585_cland_flux_mean - ssp585_ocean_mean

    fig = plt.figure()
    years = list(range(2015, 2015+NTIMES))
    plt.fill_between(years, ffemissions, color='grey', label='Fossil fuel emissions')
    plt.fill_between(years, -1*for_atmosphere-for_ocean_mean-for_cland_flux_mean, color='lightblue', label='Atmospheric accumulation')
    plt.fill_between(years, -1*for_ocean_mean-for_cland_flux_mean, color='blue', label='Ocean sink')
    plt.fill_between(years, -1*for_cland_flux_mean, color='green', label='Land sink')
    plt.hlines(0, years[0], years[-1], linestyle='dashed', color='black')
    plt.legend(loc='lower left')
    plt.xlabel('Year')
    plt.ylabel('Carbon mass [Pg(C)]')
    plt.savefig(f'{PLOTS_DIR}/carbon_budget_esm-ssp585-ssp126Lu.png')


if __name__ != 'analysis.plot_atmosphere':
    make_co2_plot()

    # Clean up temp files.
    temp_files = glob.glob('./cdoPy*')
    for f in temp_files:
        os.remove(f)

    plt.show()


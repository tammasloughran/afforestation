#!/usr/bin/env python3
"""This was originally supposed to plot the atmosphere co2 concentrations. But then it turned into
an analysis of the global carbon budget in ACCESS forrestation simulation.
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
    from analysis.cmip_files import get_filename, get_archive_filename, SSP585126LU_ENS
    from analysis.constants import (
            C_IN_CO2_RATIO,
            DATA_DIR,
            DPI,
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
            SEC_IN_DAY,
            )
else:
    # plot_afforestation.py is main program or imported as a module from another script.
    from cmip_files import get_filename, get_archive_filename, SSP585126LU_ENS
    from constants import (
            C_IN_CO2_RATIO,
            DATA_DIR,
            DPI,
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
            SEC_IN_DAY,
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
@cdod.cdo_fldmean(weights='TRUE')
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


#@cdod.cdo_mulc()
#def load_mass(var:str, input:str)->np.ma.MaskedArray:
#    """Load atmospheric CO2 as total mass in Pg(CO$_2$).
#    """
#    return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


@cdod.cdo_mul(input2=f'{DATA_DIR}/gridarea.nc') # Converts from m-2 to per gridcell
@cdod.cdo_divc(str(KG_IN_PG)) # Converts from kg to Pg
@cdod.cdo_mulc(str(SEC_IN_DAY)) # Converts from s-1 to day-1
@cdod.cdo_muldpm # Converts from day-1 to month-1
@cdod.cdo_yearsum # Converts from mon-1 to year-1
@cdod.cdo_fldsum # Converts from per gridcell to global total
@cdod.cdo_mulc(str(C_IN_CO2_RATIO)) # Convert from CO2 to C
def load_archive_carbon(var:str, input:str)->np.ma.MaskedArray:
    """Load archive carbon variable
    """
    return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


def make_co2_plot()->None:
    """Load the CO2 mixing ratio data and plot the concentration obsolute and difference.
    Also plot the values converted to atmospheric mass of C, which assumes a fixed total mass of
    the atmosphere.
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
    plt.title('Atmospheric concentration of CO$_2$')
    plt.xlabel('Year')
    plt.ylabel('CO$_2$ concentration (ppm)')
    plt.savefig(f'{PLOTS_DIR}/CO2_aff_and_ssp585.png')

    # Plot the concentration difference.
    fig = plt.figure()
    years = list(range(2015, 2015+NTIMES))
    diff = aff_co2_data - ssp585_co2_data
    diff *= KGKG_TO_MOLMOL*MIL
    plt.fill_between(years, diff.min(axis=0), diff.max(axis=0), color='royalblue', alpha=0.4)
    plt.plot(years, diff.mean(axis=0), color='royalblue', label='aff - ssp585')
    plt.hlines(0, years[0], years[-1], linewidth=0.5, color='black')
    plt.xlim(left=years[0], right=years[-1])
    plt.xlabel('Year')
    plt.ylabel('$\Delta$ CO$_2$ (ppm)')
    plt.title('Difference in Atmospheric concentration of CO$_2$')
    plt.savefig(f'{PLOTS_DIR}/CO2_atm_ppm_diff.png')

    # Plot mass difference.
    fig = plt.figure()
    years = list(range(2015, 2015+NTIMES))
    #for e,ens in enumerate(ENSEMBLES):
    #    plt.plot(
    #            years,
    #            (aff_co2_data[e,:] - ssp585_co2_data[e,:])*(MASS_ATM/KG_IN_PG)*C_IN_CO2_RATIO,
    #            color='royalblue',
    #            linewidth=0.6,
    #            alpha=0.4,
    #            )
    diff = aff_co2_data - ssp585_co2_data
    diff *= (MASS_ATM/KG_IN_PG)*C_IN_CO2_RATIO
    plt.fill_between(years, diff.min(axis=0), diff.max(axis=0), color='royalblue', alpha=0.4)
    plt.plot(
        years,
        diff.mean(axis=0),
        color='royalblue',
        label='aff - ssp585',
        )
    plt.hlines(0, years[0], years[-1], linewidth=0.5, color='black')
    plt.xlim(left=years[0], right=years[-1])
    plt.xlabel('Year')
    plt.ylabel('Carbon mass [Pg(C)]')
    plt.title('Difference in Atmospheric Mass of C')
    plt.savefig(f'{PLOTS_DIR}/CO2_atm_mass_diff.png')

    ### Carbon budget
    # This is commented out because I wrote a new carbon budget using the correct variables
    # from the ACCESS archive. See function make_carbon_budget_plot.
    ## Ocean
    #for_ocean = np.load(f'{DATA_DIR}/fgco2_aff_ocean_carbon_data.npy')*C_IN_CO2_RATIO
    #ssp585_ocean = np.load(f'{DATA_DIR}/fgco2_ssp585_ocean_carbon_data.npy')*C_IN_CO2_RATIO
    ## Land
    #for_cland = np.load(f'{DATA_DIR}/cLand_ACCESS-ESM1.5_esm-ssp585-ssp126Lu_global.npy')
    #for_cland_flux = np.concatenate([np.repeat(np.nan, 10)[:,None], np.diff(for_cland, 1)], axis=1)
    #ssp585_cland = np.load(f'{DATA_DIR}/cLand_ACCESS-ESM1.5_esm-ssp585_global.npy')
    #ssp585_cland_flux = np.concatenate(
    #        [np.repeat(np.nan, 10)[:,None], np.diff(ssp585_cland, 1)], axis=1)
    ## Fossil fuel emissions
    #emissions_file = nc.Dataset(f'{DATA_DIR}/ssp585_emissions_global.nc', 'r')
    #ffemissions = emissions_file.variables['field1561'][1:].squeeze()*C_IN_CO2_RATIO

    ## This is kind of cheating. I'm inferring the atmosphere component of the budget as the
    ## resiual of the other terms.
    #for_cland_flux_mean = for_cland_flux.mean(axis=0)
    #for_cland_flux_mean[0] = 0
    #ssp585_cland_flux_mean = ssp585_cland_flux.mean(axis=0)
    #for_ocean_mean = for_ocean.mean(axis=0)
    #ssp585_ocean_mean = ssp585_ocean.mean(axis=0)
    #for_atmosphere = ffemissions - for_cland_flux_mean - for_ocean_mean
    #ssp585_atmosphere = ffemissions - ssp585_cland_flux_mean - ssp585_ocean_mean

    #fig = plt.figure()
    #years = list(range(2015, 2015+NTIMES))
    #plt.fill_between(years, ffemissions, color='grey', label='Fossil fuel emissions')
    #plt.fill_between(
    #        years, -1*for_atmosphere-for_ocean_mean-for_cland_flux_mean,
    #        color='lightblue', label='Atmospheric accumulation')
    #plt.fill_between(years, -1*for_ocean_mean-for_cland_flux_mean, color='blue', label='Ocean sink')
    #plt.fill_between(years, -1*for_cland_flux_mean, color='green', label='Land sink')
    #plt.hlines(0, years[0], years[-1], linestyle='dashed', color='black')
    #plt.legend(loc='lower left')
    #plt.xlabel('Year')
    #plt.ylabel('Carbon flux [Pg(C)/year]')
    #plt.savefig(f'{PLOTS_DIR}/carbon_budget_esm-ssp585-ssp126Lu.png')

    # Plot the difference between my estimated GATM and the inferred residual to see how wrong
    # I am. Spoiler: I'm quite wrong.
    #plt.figure()
    #plt.plot(
    #        years[1:],
    #        np.diff(aff_co2_data.mean(axis=0), axis=0)*(MASS_ATM/KG_IN_PG)*C_IN_CO2_RATIO,
    #        color='blue',
    #        label='My calculated C to atm. converted from kg/kg to Pg(C)',
    #        )
    #ssp585_atmosphere[0] = 0
    #plt.plot(
    #        years,
    #        for_atmosphere,
    #        color='red',
    #        label='Inferred residual (EFF - SOCEAN - SLAND)',
    #        )
    #plt.title("GATM esm-ssp585")
    #plt.legend()
    #plt.show()


def make_carbon_budget_plot()->None:
    """Load emissions, and natural sinks data from the ACCESS archive data and plot the cabon
    budget for the esm-ssp595-ssp126Lu experiment.
    """
    variables = [
            'socean',
            'sland',
            'gatm',
            'eff',
            ]
    ARCHIVE_DIR = '/g/data/p66/tfl561/archive_data'
    # Load data.
    var_data = {}
    for var in variables:
        if load_cdo==True:
            var_data[var] = np.ones((NENS,NTIMES))*np.nan
            for e,ens in enumerate(SSP585126LU_ENS):
                filename = f'{ARCHIVE_DIR}/{var}_{ens}.nc'
                var_data[var][e,:] = load_archive_carbon(var=var, input=filename)
            np.save(f'{DATA_DIR}/{var}_PgCyear.npy', var_data[var])
        else:
            var_data[var] = np.load(f'{DATA_DIR}/{var}_PgCyear.npy')

    # Create stacked ensemble means for budget
    eff = var_data['eff'].mean(axis=0)
    gatm = -var_data['gatm'].mean(axis=0) - var_data['socean'].mean(axis=0) \
            + var_data['sland'].mean(axis=0)
    socean = var_data['sland'].mean(axis=0) - var_data['socean'].mean(axis=0)
    sland = var_data['sland'].mean(axis=0)

    # Plot data.
    plt.figure()
    years = list(range(2015, 2015+len(sland)))
    plt.fill_between(years, eff, color='grey', label='Fossil fuel emissions')
    plt.fill_between(years, gatm, color='lightblue', label='Atmospheric accumulation')
    plt.fill_between(years, socean, color='blue', label='Ocean sink')
    plt.fill_between(years, sland, color='green', label='Land sink')
    plt.xlim(left=years[0], right=years[-1])
    plt.xlabel('Years')
    plt.ylabel('Pg(C)/year')
    plt.legend(loc='upper left', frameon=False)
    plt.title('esm-ssp585-ssp126Lu Global Carbon Budget')
    plt.savefig(f'{PLOTS_DIR}/esm-ssp585-ssp126Lu_budget.png', dpi=DPI)


if __name__ != 'analysis.plot_atmosphere':
    #make_co2_plot()
    make_carbon_budget_plot()

    # Clean up temp files.
    temp_files = glob.glob('./cdoPy*')
    for f in temp_files:
        os.remove(f)

    plt.show()


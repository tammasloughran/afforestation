#!/usr/bin/env python3
"""Plot the atmosphere co2 concentrations for all models.
"""
import os
import sys
import glob
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
from cdo import Cdo
import cdo_decorators as cdod
if __name__ == 'analysis.plot_afforestation':
    # plot_afforestation.py imported as a module of the analysis package.
    from analysis.cmip_files import get_filenames
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
    from cmip_files import get_filenames
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
color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

# Constants
CO2_VARIABLES = {
        'Emon':'co23D', # CO2 concentration (mixing ratio) in kg/kg.
        'Amon':'co2', # Mole fraction of CO2. esm-ssp585 only
        }
MODELS = [ # Some models have missing co2 data for one of the experiments.
        'ACCESS-ESM1-5',
        'BCC-CSM2-MR',
        #'CanESM5',
        #'GFDL-ESM4',
        'MIROC-ES2L',
        'MPI-ESM1-2-LR',
        #'NorESM2-LM', # NorESM has been removed because it was run in concentration driven mode.
        #'UKESM1-0-LL',
        'CESM2',
        ]
INSTIT = {
        'ACCESS-ESM1-5':'CSIRO',
        'BCC-CSM2-MR':'BCC',
        'CanESM5':'CCma',
        'GFDL-ESM4':'NOAA-GFDL',
        'MIROC-ES2L':'MIROC',
        'MPI-ESM1-2-LR':'MPI-M',
        #'NorESM2-LM':'NCC',
        'UKESM1-0-LL':'MOHC',
        'CESM2':'NCAR',
        }
ENSEMBLES = {
        'ACCESS-ESM1-5':'r1i1p1f1',
        'BCC-CSM2-MR':'r1i1p1f1',
        'CanESM5':'r1i1p2f1',
        'MIROC-ES2L':'r1i1p1f2',
        'UKESM1-0-LL':'r1i1p1f2',
        'MPI-ESM1-2-LR':'r1i1p1f1',
        #'NorESM2-LM':'r1i1p1f1',
        'GFDL-ESM4':'r1i1p1f1',
        'CESM2':'r1i1p1f1',
        }
A_OR_AER = {
        'ACCESS-ESM1-5':'Amon',
        'BCC-CSM2-MR':'Amon',
        'CanESM5':'Amon',
        'GFDL-ESM4':'Amon',
        'MIROC-ES2L':'AERmon',
        'MPI-ESM1-2-LR':'Amon',
        #'NorESM2-LM':'AERmon',
        'UKESM1-0-LL':'Amon',
        'CESM2':'Amon',
        }
COLORS = {
        'CSIRO':color_cycle[0],
        'BCC':color_cycle[1],
        'CCma':color_cycle[2],
        'NOAA-GFDL':color_cycle[3],
        'MIROC':color_cycle[4],
        'MPI-M':color_cycle[5],
        #'NCC':color_cycle[6],
        'NCAR':color_cycle[6],
        'MOHC':color_cycle[8],
        }
TO_PPM = 1000000

# CMIP6 CO2 concentration data taken from
# https://tntcat.iiasa.ac.at/SspDb/dsd?Action=htmlpage&page=10
remind_magpie_co2 = [
        379.850,
        390.505,
        417.249,
        452.767,
        499.678,
        559.692,
        635.793,
        730.025,
        841.520,
        963.842,
        1088.970,
        ]
remind_magpie_years = [
        2005,
        2010,
        2020,
        2030,
        2040,
        2050,
        2060,
        2070,
        2080,
        2090,
        2100,
        ]


load_cdo = False
load_npy = not load_cdo


@cdod.cdo_cat(input2='')
@cdod.cdo_vertmean
@cdod.cdo_fldmean(weights='TRUE')
@cdod.cdo_yearmonmean
def global_average(var:str, input:str)->np.ma.MaskedArray:
    """Load files field and vertical aggregated yearly data using CDO.
    """
    return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


def make_co2_models_plot()->None:
    """Make the model intercomparison plot of co2 concentration.
    """
    global for_co2_data
    global ssp585_co2_data
    for_co2_data = {}
    ssp585_co2_data = {}

    # Load data
    for model in MODELS:
        print(f"Loading model {model}")
        if load_cdo:
            if model=='ACCESS-ESM1-5':
                print(model, 'will not be loaded from cdo. Use .npy files.')
                continue
            table = A_OR_AER[model]
            for_co2_files = sorted(get_filenames(
                    'LUMIP',
                    INSTIT[model],
                    model,
                    'esm-ssp585-ssp126Lu',
                    ENSEMBLES[model],
                    table,
                    'co2',
                    ))
            for_co2_files = '[ ' + ' '.join(for_co2_files) + ' ]'
            for_co2_data[model] = global_average(var='co2', input=for_co2_files)

            ssp585_co2_files = sorted(get_filenames(
                    'C4MIP',
                    INSTIT[model],
                    model,
                    'esm-ssp585',
                    ENSEMBLES[model],
                    table,
                    'co2',
                    ))
            ssp585_co2_files = '[ ' + ' '.join(ssp585_co2_files) + ' ]'
            ssp585_co2_data[model] = global_average(var='co2', input=ssp585_co2_files)

            np.save(f'{DATA_DIR}/{model}_co2_for.npy', for_co2_data[model].data)
            np.save(f'{DATA_DIR}/{model}_co2_ssp585.npy', ssp585_co2_data[model].data)
        else:
            if model=='ACCESS-ESM1-5':
                try:
                    access_for_co2 = np.load(f'{DATA_DIR}/aff_co2_data.npy')
                except:
                    print('Need to create ACCESS-ESM1-5 co2 .npy file first.')
                    sys.exit(1)
                for_co2_data[model] = access_for_co2.mean(axis=0)*KGKG_TO_MOLMOL
                access_ssp585_co2 = np.load(f'{DATA_DIR}/ssp585_co2_data.npy')
                ssp585_co2_data[model] = access_ssp585_co2.mean(axis=0)*KGKG_TO_MOLMOL
                access_diff_co2 = access_for_co2 - access_ssp585_co2
            else:
                for_co2_data[model] = np.load(f'{DATA_DIR}/{model}_co2_for.npy')
                ssp585_co2_data[model] = np.load(f'{DATA_DIR}/{model}_co2_ssp585.npy')

    # Plot absolute values.
    plt.figure
    for model in MODELS:
        years = np.arange(2015, 2015+len(for_co2_data[model]))
        plt.plot(
                years,
                for_co2_data[model]*TO_PPM,
                color=COLORS[INSTIT[model]],
                label=model+' forestation',
                )
        years = np.arange(2015, 2015+len(ssp585_co2_data[model]))
        plt.plot(
                years,
                ssp585_co2_data[model]*TO_PPM,
                color=COLORS[INSTIT[model]],
                label=model+' esm-ssp585',
                linestyle='dashed',
                )
    plt.plot(
            remind_magpie_years,
            remind_magpie_co2,
            color='black',
            label='REMIND-MAGPIE SSP5-8.5',
            linestyle='solid',
            )
    plt.xlim(2010, 2100)
    plt.xlabel('Years')
    plt.ylabel('CO$_2$ mixing ratio (ppm)')
    plt.legend()
    plt.savefig(f'{PLOTS_DIR}/models/models_co2.png', dpi=DPI)

    # Plot relative to ssp585
    plt.figure()
    for model in MODELS:
        # MPI is missing a year in ssp585 and NorESM is missing a year in forestation scenario.
        if model=='MPI-ESM1-2-LR':
            for_co2_data[model] = for_co2_data[model][:-1]
        if model=='NorESM2-LM':
            ssp585_co2_data[model] = ssp585_co2_data[model][:-1]
            continue # Actually don't even plot NorESM here. It's the wrong data anyway.
        years = np.arange(2015, 2015+len(for_co2_data[model]))
        plt.plot(
                years,
                for_co2_data[model]*TO_PPM - ssp585_co2_data[model]*TO_PPM,
                color=COLORS[INSTIT[model]],
                label=model,
                )
    plt.fill_between(
            np.arange(2015, 2015+len(access_for_co2[0,:])),
            access_diff_co2.min(axis=0)*TO_PPM*KGKG_TO_MOLMOL,
            access_diff_co2.max(axis=0)*TO_PPM*KGKG_TO_MOLMOL,
            color=COLORS['CSIRO'],
            alpha=0.5,
            )
    plt.hlines(0, 2015, 2100, color='black', linewidth=0.5)
    plt.xlim(2015, 2100)
    plt.xlabel('Years')
    plt.ylabel('$\Delta$ CO$_2$ mixing ratio (ppm)')
    plt.legend()
    plt.title('Difference between forestation and esm-ssp585')
    plt.savefig(f'{PLOTS_DIR}/models/models_co2_diff.png', dpi=DPI)


if __name__ != 'analysis.plot_models_co2':
    make_co2_models_plot()

    plt.show()


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
        #'BCC-CSM2-MR', # BCC-CSM2 has beed excluded at the request of BCC due to a bug.
        'CanESM5',
        #'GFDL-ESM4',
        'MIROC-ES2L',
        'MPI-ESM1-2-LR',
        #'NorESM2-LM', # NorESM has been removed because it was run in concentration driven mode.
        #'UKESM1-0-LL',
        'CESM2',
        ]
INSTIT = {
        'ACCESS-ESM1-5':'CSIRO',
        #'BCC-CSM2-MR':'BCC',
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
        'CanESM5':'r1i1p1f1',
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
        'NorESM2-LM':'AERmon',
        'UKESM1-0-LL':'Amon',
        'CESM2':'Amon',
        }
COLORS = {
        'CSIRO':color_cycle[0],
        #'BCC':color_cycle[1],
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

# UKESM data provided my Spencer Liddicoat
UKESM_co2 = [
        0.08184814, 0.03216553, 0.13171387, 0.16906738, 0.15847778, -0.02719116,
        0.28515625, 1.0713501, 1.5263062, 0.8100281, -0.2937622, -1.1992493, -0.73062134,
        -1.2385254, -1.44281, -1.4172974, -1.751648, -1.5945435, -2.4179688, -1.9859009,
        -1.3006592, -2.712555, -3.5486145, -4.0009155, -5.8648376, -6.3403015, -5.470093,
        -5.6828613, -5.6918945, -6.2141724, -7.425293, -7.461731, -7.150696, -8.360596,
        -9.320862, -9.034668, -8.496155, -9.218567, -10.516418, -10.931885, -9.960632,
        -9.759094, -9.646362, -9.090576, -7.7731323, -8.1466675, -9.414246, -9.557129,
        -9.194275, -9.434387, -8.627197, -8.28833, -9.507874, -8.981812, -8.380676, -9.296509,
        -10.569092, -12.085999, -12.627991, -13.179993, -14.441223, -15.57135, -15.703369,
        -15.689697, -14.785095, -15.866455, -14.553223, -13.411682, -12.780579, -13.000549,
        -13.348633, -12.764221, -10.759094, -10.843323, -11.480164, -10.855286, -11.055664,
        -13.439331, -14.272949, -12.908691, -11.814209, -12.122925, -11.6427, -11.89624,
        -11.100464, -10.079956,
        ]
UKESM_for=[
        409.86435, 412.86044, 415.64072, 418.2894, 421.2275, 424.06015, 427.32812, 431.074,
        434.88058, 437.86197, 440.70035, 443.98523, 448.55127, 452.0446, 455.5525, 459.756,
        463.45178, 467.59286, 471.4693, 476.62518, 481.93414, 485.70074, 490.51297, 495.59457,
        500.6811, 506.1226, 511.8401, 517.7143, 524.1846, 530.72314, 535.94226, 541.8053,
        548.64575, 555.2149, 562.07184, 569.54956, 577.0409, 584.3829, 591.7246, 599.5155,
        608.2901, 617.00494, 625.7843, 635.3482, 644.8851, 653.5727, 662.6138, 672.1841,
        682.93726, 693.18463, 703.55115, 714.43835, 725.6508, 737.3981, 748.36426, 758.944,
        770.6142, 783.0665, 794.81775, 806.43, 818.297, 830.9239, 843.0899, 855.3816, 868.58734,
        882.06555, 896.7623, 910.51666, 924.47833, 938.6327, 951.9604, 965.9483, 981.2608,
        995.74603, 1008.3736, 1022.0523, 1036.5997, 1049.4207, 1062.49, 1076.915, 1092.3136,
        1106.9558, 1121.4246, 1134.8689, 1148.6837, 1163.917,
        ]
UKESM_ssp585 = [
        409.7825, 412.82828, 415.509, 418.12033, 421.06903, 424.08734,
        427.04297, 430.00266, 433.35428, 437.05194, 440.9941, 445.18448,
        449.2819, 453.2831, 456.9953, 461.1733, 465.20343, 469.1874,
        473.88727, 478.61108, 483.2348, 488.4133, 494.06158, 499.5955,
        506.54593, 512.4629, 517.3102, 523.39716, 529.87646, 536.9373,
        543.36755, 549.267, 555.79645, 563.5755, 571.3927, 578.5842,
        585.53705, 593.60144, 602.241, 610.4474, 618.25073, 626.76404,
        635.43066, 644.4388, 652.6582, 661.71936, 672.028, 681.7412,
        692.13153, 702.619, 712.17834, 722.7267, 735.1587, 746.3799,
        756.74493, 768.2405, 781.1833, 795.1525, 807.44574, 819.61,
        832.7382, 846.49524, 858.7933, 871.0713, 883.37244, 897.932,
        911.31555, 923.92834, 937.2589, 951.63324, 965.309, 978.7125,
        992.0199, 1006.58936, 1019.85376, 1032.9076, 1047.6554, 1062.86,
        1076.763, 1089.8237, 1104.1278, 1119.0787, 1133.0673, 1146.7651,
        1159.7842, 1173.997,
        ]
CanESM5_for = [
        398.86, 401.59, 403.82, 406.59, 409.48, 412.67, 415.76, 417.95, 421.11, 424.80, 428.67,
        431.86, 433.82, 437.76, 442.69, 446.11, 449.16, 453.16, 458.66, 463.42, 466.84, 469.81,
        473.84, 478.93, 484.52, 487.80, 491.31, 496.70, 501.96, 506.75, 511.32, 516.78, 522.98,
        530.46, 535.77, 540.30, 545.33, 550.92, 556.86, 562.94, 571.41, 578.85, 585.32, 593.04,
        599.48, 605.32, 612.74, 619.51, 627.41, 635.36, 642.76, 650.40, 659.51, 668.86, 677.78,
        687.04, 697.39, 706.45, 714.96, 724.45, 734.61, 743.10, 752.27, 762.36, 773.11, 783.87,
        792.66, 801.54, 813.25, 824.13, 834.34, 844.55, 855.04, 865.47, 875.17, 885.38, 893.93,
        904.13, 914.35, 925.03, 935.78, 947.10, 958.02, 967.03, 976.36, 986.65,
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
            if model!='CanESM5':
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
            else:
                for_co2_data['CanESM5'] = np.array(CanESM5_for)/TO_PPM

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
            if model=='CanESM5': ssp585_co2_data[model] /= TO_PPM

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
    plt.figure()
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
    # Manually add the UKESM b/c since it was missing from the ESGF. Provided by Spencer Liddicoat.
    plt.plot(
            years,
            UKESM_for,
            color=COLORS['MOHC'],
            label='UKESM1-0-LL forestation',
            )
    plt.plot(
            years,
            UKESM_ssp585,
            color=COLORS['MOHC'],
            label='UKESM1-0-LL esm-ssp585',
            linestyle='dashed',
            )
    # Manually add the CanESM5 data. Provided by Vivek Aurora.
    #plt.plot(
    #        years,
    #        CanESM5_for,
    #        color=COLORS['CCma'],
    #        label='CanESM5 forestation',
    #        )
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
    plt.legend(frameon=False)
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
    plt.plot(
            years,
            UKESM_co2,
            color=COLORS['MOHC'],
            label='UKESM1-0-LL',
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
    plt.legend(frameon=False)
    plt.title('Difference between forestation and esm-ssp585')
    plt.savefig(f'{PLOTS_DIR}/models/models_co2_diff.png', dpi=DPI)


if __name__ != 'analysis.plot_models_co2':
    make_co2_models_plot()

    plt.show()


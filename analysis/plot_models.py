#!/usr/bin/env python3
# plot_afforestation.py generates plots for the afforestation simulation.
# In general "aff" is used to refer to the afforestation simulation esm-ssp585-ssp126Lu (high
# fossil fuels emissions and low land use change scenario), and "ssp585" means the emissions driven
# esm-ssp585 simulation (high fossil fuel emisssions scenario).
import glob
import os
import ipdb
import sys

import cdo as cdo_module
import matplotlib.pyplot as plt
import numpy as np
import cdo_decorators as cdod
import pymannkendall as pmk
import ipdb

if __name__ == 'analysis.plot_afforestation':
    # plot_afforestation.py imported as a module of the analysis package.
    from analysis.cmip_files import get_filenames
    from analysis.constants import (
            DATA_DIR,
            DPI,
            ENSEMBLES,
            KG_IN_PG,
            NENS,
            PLOTS_DIR,
            SEC_IN_DAY,
            TABLES,
            )
else:
    # plot_afforestation.py is main program or imported as a module from another script.
    from cmip_files import get_filenames
    from constants import (
            DATA_DIR,
            DPI,
            ENSEMBLES,
            KG_IN_PG,
            NENS,
            PLOTS_DIR,
            SEC_IN_DAY,
            TABLES,
            )

cdo = cdo_module.Cdo()
cdo.debug = False

# Local constants
COLORS = {
        'gpp':'green',
        'npp':'olive',
        'ra':'orange',
        'rh':'saddlebrown',
        'nbp':'purple',
        'cVeg':'darkgreen',
        'cLitter':'chocolate',
        'cSoil':'black',
        'tas':'black',
        'pr':'blue',
        }

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
load_npy_files = True # Uncomment to override previous check.

color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

MODELS = { # ACCESS is excluded here. Needs separate plotting for ensembles.
        #'CSIRO':'ACCESS-ESM1-5',
        #'BCC':'BCC-CSM2-MR', # BCC has been excluded.
        #'NCC':'NorESM2-LM', # NorESM has been excluded.
        'CCma':'CanESM5',
        'MIROC':'MIROC-ES2L',
        'MOHC':'UKESM1-0-LL',
        'MPI-M':'MPI-ESM1-2-LR',
        'NCAR':'CESM2',
        'NOAA-GFDL':'GFDL-ESM4',
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
ENSEMBLES = {
        #'BCC':'r1i1p1f1',
        #'NCC':'r1i1p1f1',
        'CCma':'r1i1p2f1',
        'MIROC':'r1i1p1f2',
        'MOHC':'r1i1p1f2',
        'MPI-M':'r1i1p1f1',
        'NCAR':'r1i1p1f1',
        'NOAA-GFDL':'r1i1p1f1',
        }
SSP585_ENSEMBLES = {
        #'BCC':'r1i1p1f1',
        #'NCC':'r1i1p1f1',
        'CCma':'r1i1p1f1',
        'MIROC':'r1i1p1f2',
        'MOHC':'r1i1p1f2',
        'MPI-M':'r1i1p1f1',
        'NCAR':'r1i1p1f1',
        'NOAA-GFDL':'r1i1p1f1',
        }
VARIABLES = {
        'Lmon':[
            'gpp',
            'nbp',
            'npp',
            'ra',
            'rh',
            'cLitter',
            'cVeg',
            ],
        'Emon':[
            'cLand',
            'cSoil',
            ],
        'Amon':[
            'pr',
            'tas',
            ],
        }


def make_model_plots()->None:
    """Create plots of carbon pools and climate trends for all models.
    """
    print('Starting model intercomparison analysis.')
    if not os.path.exists(f'{PLOTS_DIR}/models'): os.mkdir(f'{PLOTS_DIR}/models')
    # Create container dicts and variables for aggregating soil and litter pools.
    litter_soil_for = {}
    litter_soil_ssp585 = {}
    for m in MODELS.keys():
        litter_soil_for[m] = 0
        litter_soil_ssp585[m] = 0
    access_litter_soil_for = 0
    access_litter_soil_ssp585 = 0
    for table in VARIABLES.keys():
        for var in VARIABLES[table]:
            print(f'\nProcessing {var}:')
            aff_data = {}
            ssp585_data = {}
            for instit in MODELS.keys():
                if instit=='MOHC' and var=='cLitter': continue # Skip missing data.
                if instit=='BCC' and var=='nbp': continue # Skip missing data.
                if instit=='NOAA-GFDL' and not table=='Amon': continue # GFDL only has pr and tas
                print('    -', instit, MODELS[instit])
                # Create the gridarea.nc file for this model.
                data_files = glob.glob(f'{DATA_DIR}/*')
                if f'{DATA_DIR}/gridarea_{MODELS[instit]}.nc' not in data_files:
                    print('       - Making grid area file')
                    afile = get_filenames(
                            'LUMIP',
                            instit,
                            MODELS[instit],
                            'esm-ssp585-ssp126Lu',
                            ENSEMBLES[instit],
                            table,
                            var,
                            )[0]
                    cdo.gridarea(input=afile, output=f'{DATA_DIR}/gridarea_{MODELS[instit]}.nc')

                model_land_frac = f'/g/data/p66/tfl561/CMIP6/C4MIP/{instit}/{MODELS[instit]}' \
                        f'/esm-ssp585/{SSP585_ENSEMBLES[instit]}/fx/sftlf/gn/latest' \
                        f'/sftlf_fx_{MODELS[instit]}_esm-ssp585_{SSP585_ENSEMBLES[instit]}_gn.nc'
                if instit=='BCC': # BCC land frac file is 0-1
                    frac_unit = 1
                else:
                    frac_unit = 100


                # The loader functions need to be defined in this loop to account for the model
                # resolution.
                @cdod.cdo_cat(input2='')
                @cdod.cdo_mul(input2=model_land_frac)
                @cdod.cdo_divc(str(frac_unit))
                @cdod.cdo_mul(input2=f'{DATA_DIR}/gridarea_{MODELS[instit]}.nc')
                @cdod.cdo_fldsum
                @cdod.cdo_yearmonmean
                @cdod.cdo_divc(str(KG_IN_PG))
                def cdo_pool_load_model(var:str, input:str)->np.ma.MaskedArray:
                    """Load global climate variable using CDO. Please refer to the following
                    decorators:
                    @cdod.cdo_cat(input2='')
                    @cdod.cdo_mul(input2=model_land_frac)
                    @cdod.cdo_divc(str(frac_unit))
                    @cdod.cdo_mul(input2=f'{DATA_DIR}/gridarea_{MODELS[instit]}.nc')
                    @cdod.cdo_fldsum
                    @cdod.cdo_yearmonmean
                    @cdod.cdo_divc(str(KG_IN_PG))
                    """
                    return cdo.copy(
                            input=input,
                            options='-L',
                            returnCdf=True,
                            ).variables[var][:].squeeze()


                @cdod.cdo_cat(input2='')
                @cdod.cdo_mul(input2=model_land_frac)
                @cdod.cdo_divc(str(frac_unit))
                @cdod.cdo_mul(input2=f'{DATA_DIR}/gridarea_{MODELS[instit]}.nc')
                @cdod.cdo_fldsum
                @cdod.cdo_mulc(str(SEC_IN_DAY))
                @cdod.cdo_muldpm
                @cdod.cdo_yearsum
                @cdod.cdo_divc(str(KG_IN_PG))
                def cdo_flux_load_model(var:str, input:str)->np.ma.MaskedArray:
                    """Load global climate variable using CDO. Please refer to the following
                    decorators:
                    @cdod.cdo_cat(input2='')
                    @cdod.cdo_mul(input2=model_land_frac)
                    @cdod.cdo_divc(str(frac_unit))
                    @cdod.cdo_mul(input2=f'{DATA_DIR}/gridarea_{MODELS[instit]}.nc')
                    @cdod.cdo_fldsum
                    @cdod.cdo_mulc(str(SEC_IN_DAY))
                    @cdod.cdo_muldpm
                    @cdod.cdo_yearsum
                    @cdod.cdo_divc(str(KG_IN_PG))
                    """
                    return cdo.copy(
                            input=input,
                            options='-L',
                            returnCdf=True,
                            ).variables[var][:].squeeze()


                @cdod.cdo_cat(input2='')
                #@cdod.cdo_masklonlatbox('-180','180','-60','90') # Exclude Antarctica
                #@cdod.cdo_ifthen(input1=model_land_frac) # Land only. Comment for land+ocean
                @cdod.cdo_fldmean(weights='TRUE')
                @cdod.cdo_yearmonmean
                def cdo_clim_load_model(var:str, input:str)->np.ma.MaskedArray:
                    """Load global climate variable using CDO. Please refer to the following
                    decorators:
                    @cdod.cdo_cat(input2='')
                    #@cdod.cdo_masklonlatbox('-180','180','-60','90') # Mask to remove Antarctica.
                    @cdod.cdo_ifthen(input1=LAND_FRAC_FILE) # Mask for climate over land only.
                    @cdod.cdo_fldmean(weights='TRUE')
                    @cdod.cdo_yearmonmean
                    """
                    return cdo.copy(
                            input=input,
                            options='-L',
                            returnCdf=True,
                            ).variables[var][:].squeeze()


                # Select pool or climate loader for this variable.
                if var in ['pr', 'tas']:
                    loader = cdo_clim_load_model
                elif var in ['gpp','npp','nbp','ra','rh']:
                    loader = cdo_flux_load_model
                else:
                    loader = cdo_pool_load_model

                # Load the CMIP6 models.
                if not load_npy_files:
                    file_list = sorted(get_filenames(
                            'LUMIP',
                            instit,
                            MODELS[instit],
                            'esm-ssp585-ssp126Lu',
                            ENSEMBLES[instit],
                            table,
                            var,
                            ))
                    filenames = '[ ' + ' '.join(file_list) + ' ]'
                    aff_data[instit] = loader(input=filenames, var=var)
                    np.save(
                            f'{DATA_DIR}/{var}_{MODELS[instit]}_aff_global.npy',
                            aff_data[instit].data,
                            )
                    file_list = sorted(get_filenames(
                            'C4MIP',
                            instit,
                            MODELS[instit],
                            'esm-ssp585',
                            SSP585_ENSEMBLES[instit],
                            table,
                            var,
                            ))
                    filenames = '[ ' + ' '.join(file_list) + ' ]'
                    ssp585_data[instit] = loader(input=filenames, var=var)
                    np.save(
                            f'{DATA_DIR}/{var}_{MODELS[instit]}_ssp585_global.npy',
                            ssp585_data[instit].data,
                            )
                else:
                    aff_data[instit] = np.load(
                            f'{DATA_DIR}/{var}_{MODELS[instit]}_aff_global.npy',
                            )
                    ssp585_data[instit] = np.load(
                            f'{DATA_DIR}/{var}_{MODELS[instit]}_ssp585_global.npy',
                            )
                # Correct units.
                if var=='tas' and aff_data[instit][0]>100: aff_data[instit] -= 273.15
                if var=='tas' and ssp585_data[instit][0]>100: ssp585_data[instit] -= 273.15

                # Aggregate soil and litter pools
                if var=='cSoil':
                    litter_soil_for[instit] += aff_data[instit]
                    litter_soil_ssp585[instit] += ssp585_data[instit]
                if var=='cLitter':
                    litter_soil_for[instit] += aff_data[instit]
                    litter_soil_ssp585[instit] += ssp585_data[instit]


            # Load the ACCESS-ESM1-5 data. I expect that the npy files already exist.
            access_aff = np.load(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585-ssp126Lu_global.npy')
            access_ssp585 = np.load(f'{DATA_DIR}/{var}_ACCESS-ESM1.5_esm-ssp585_global.npy')
            if var=='tas': access_aff -= 273.15
            if var=='tas': access_ssp585 -= 273.15
            diff = access_aff - access_ssp585
            ens_mean = np.mean(diff, axis=0)

            # Aggregate ACCESS soil and litter.
            if var=='cSoil':
                access_litter_soil_for += access_aff
                access_litter_soil_ssp585 += access_ssp585
            if var=='clitter':
                access_litter_soil_for += access_aff
                access_litter_soil_ssp585 += access_ssp585

            # Some models have a slightly different time period for each experiment.
            for m in MODELS.keys():
                if m=='MOHC' and var=='cLitter':
                    continue # Skip missing data.
                if m=='MPI-M':
                    aff_data[m] = aff_data[m][:-1]
                if m=='NCC':
                    ssp585_data[m] = ssp585_data[m][:-1]

            # Plot the difference between simulations for cpools.
            if not var=='tas' and not var=='pr':
                plt.figure()
                # Plot ACCESS ensemble mean and range.
                years = list(range(2015, 2015 + ens_mean.shape[0]))
                plt.plot(years, ens_mean, label='ACCESS-ESM1-5', color=COLORS['CSIRO'])
                plt.fill_between(years, diff.max(axis=0), diff.min(axis=0), color='lightblue')
                # The following plots individual ensemble members as gray lines.
                #for e in range(10):
                #    plt.plot(years, diff[e,:], color='gray', alpha=0.4)
                # Plot other models.
                for m in MODELS.keys():
                    if m=='MOHC' and var=='cLitter': continue
                    if m=='BCC' and var=='nbp': continue # Skip missing data.
                    if m=='NOAA-GFDL': continue # GFDL only has pr and tas
                    try:
                        diff_model = aff_data[m] - ssp585_data[m]
                    except:
                        ipdb.set_trace()
                    if m=='CCma' and (var=='cVeg' or var=='cSoil'):
                        # CanESM5 cVeg has a large initial bias in cVeg and cSoil. Remove this bias.
                        diff_model = diff_model - diff_model[0]
                    years = list(range(2015, 2015 + aff_data[m].shape[0]))
                    plt.plot(years, diff_model, color=COLORS[m], label=MODELS[m])
                # Plot features.
                plt.xlim(left=years[0], right=years[-1])
                plt.hlines(0, years[0], years[-1], color='black', linewidth=0.5)
                if var in ['gpp', 'npp', 'nbp', 'ra', 'rh']:
                    plt.ylabel('Pg(C)/year')
                else:
                    plt.ylabel('Pg(C)')
                plt.xlabel('Year')
                plt.title(f"{var} esm-ssp585-ssp126Lu - esm-ssp585")
                plt.legend(frameon=False)
                plt.savefig(f'{PLOTS_DIR}/models/{var}_model_intercomparison_diff.png', dpi=DPI)
                plt.close()

            # Plot the trends for tas and pr.
            if var=='tas' or var=='pr':
                trend = np.array(['']*NENS)
                h = np.ones(NENS).astype(bool)
                p = np.ones(NENS)*np.nan
                z = p.copy()
                tau = p.copy()
                s = p.copy()
                var_s = p.copy()
                slope = p.copy()
                intercept = p.copy()
                # Regression of difference first
                for e in range(10):
                    trend[e], h[e], p[e], z[e], tau[e], s[e], var_s[e], slope[e], intercept[e] = \
                                pmk.original_test(
                                        diff[e,:], # ACCESS-ESM1.5
                                        alpha=0.05,
                                        )
                # Regression of experiments first, then difference.
                # It's not a good idea to do this because the temperature change for each
                # experiment is not linear but the difference between them is.
                #for e in range(10):
                #    trenda, ha, pa, za, taua, sa, var_sa, slopea, intercepta = \
                #                pmk.original_test(
                #                        access_aff[e,:], # ACCESS-ESM1.5
                #                        alpha=0.05,
                #                        )
                #    years = np.arange(ens_mean.shape[0])
                #    aff_trend = slopea*years+intercepta
                #    trendb, hb, pb, zb, taub, sb, var_sb, slopeb, interceptb = \
                #                pmk.original_test(
                #                        access_ssp585[e,:], # ACCESS-ESM1.5
                #                        alpha=0.05,
                #                        )
                #    ssp585_trend = slopeb*years+interceptb
                #
                #    difference_of_trends = aff_trend - ssp585_trend
                #
                #    plt.figure()
                #    plt.plot(years, access_aff[e,:], label='aff data')
                #    plt.plot(years, access_ssp585[e,:], label='ssp585 data')
                #    plt.plot(years, aff_trend, label='aff_trend')
                #    plt.plot(years, ssp585_trend, label='ssp585_trend')
                #    plt.plot(years, slopec*years+interceptc, label='regress first')
                #    plt.plot(years, access_aff[e,:], label='aff')
                #    plt.plot(years, access_ssp585[e,:], label='ssp585')
                #    plt.plot(years, trend_line[e,:])
                #    plt.plot(years, diff[e,:], color='grey')
                #    plt.plot(years, slope[e]*years+intercept[e], label='diff first')
                #    plt.legend()
                #    plt.show()

                years = list(range(2015, 2015 + ens_mean.shape[0]))
                trend_line = slope.mean()*np.arange(len(years)) + intercept.mean()
                trend_lines = slope[:,None]*np.arange(len(years))[None,:] + intercept[:,None]
                plt.figure()
                for e in range(10):
                    if h[e]:
                        plt.plot(years, diff[e,:], color='orange')
                        plt.plot(years, trend_lines[e,:], color='red')
                    else:
                        plt.plot(years, diff[e,:], color='grey')
                        plt.plot(years, trend_lines[e,:], color='black')
                plt.show()

                fig, axes = plt.subplots(nrows=4, ncols=2, sharex=True, sharey=True)
                axes = axes.flatten()
                axes[0].plot(years, ens_mean, label='ACCESS-ESM1-5')
                axes[0].fill_between(years, diff.max(axis=0), diff.min(axis=0), color='lightblue')
                if h.all():
                    print("    - All ACCESS trends are significant")
                    axes[0].plot(years, trend_line, color=COLORS['CSIRO'])
                else:
                    print("    - Not all ACCESS trends are significant")
                    print("        - They are:", h)
                    print("        - deltas: ", trend_lines[:,-1] - trend_lines[:,0])
                    axes[0].plot(years, trend_line, color=COLORS['CSIRO'], linestyle='dotted')
                if slope.mean()>0: sign = '+'
                else: sign = '-'
                axes[0].annotate(
                        f'ACCESS-ESM1-5 ({sign})',
                        (0.3,1),
                        xycoords='axes fraction',
                        fontsize=8,
                        )
                axes[0].hlines(0, years[0], years[-1], color='black', linewidth=0.5)
                axes[1].set_axis_off() # Disable top right axes for odd number of models.
                for i,m in enumerate(MODELS.keys()):
                    if m=='MPI-M' and var=='pr': aff_data[m] = aff_data[m][0:-1]
                    diff_model = aff_data[m] - ssp585_data[m]
                    if m=='CCma' and var=='tas':
                        # CanESM5 has a large initial bias. Remove this bias.
                        diff_model = diff_model - diff_model[0]
                    years = list(range(2015, 2015 + aff_data[m].shape[0]))
                    axes[i+2].plot(years, diff_model, color=COLORS[m], label=MODELS[m])
                    #if m=='MIROC':
                    #    trenda, ha, pa, za, taua, sa, var_sa, slopea, intercepta = \
                    #            pmk.original_test(
                    #                    aff_data[m],
                    #                    alpha=0.05,
                    #                    )
                    #    trendb, hb, pb, zb, taub, sb, var_sb, slopeb, interceptb = \
                    #            pmk.original_test(
                    #                    ssp585_data[m],
                    #                    alpha=0.05,
                    #                    )
                    #    trend_linea = slopea*np.arange(len(years)) + intercepta
                    #    trend_lineb = slopeb*np.arange(len(years)) + interceptb
                    #    trend_line = trend_linea - trend_lineb
                    #    slope = slopea - slopeb
                    #else:
                    trend, h, p, z, tau, s, var_s, slope, intercept = pmk.original_test(
                            diff_model,
                            alpha=0.05,
                            )
                    trend_line = slope*np.arange(len(years)) + intercept
                    print(f'    - {MODELS[m]} trend={trend}, p={p}, h={h}, delta={trend_line[-1]-trend_line[0]}')
                    if h: # Hypothesis that there exists a trend is true.
                        axes[i+2].plot(years, trend_line, color=COLORS[m])
                    else:
                        axes[i+2].plot(years, trend_line, color=COLORS[m],
                                linestyle='dotted')
                    if slope>0: sign = '+'
                    else: sign = '-'
                    axes[i+2].annotate(
                            f'{MODELS[m]} ({sign})',
                            (0.3,1),
                            xycoords='axes fraction',
                            fontsize=8,
                            )
                    axes[i+2].set_xlim(left=years[0], right=years[-1])
                    axes[i+2].hlines(0, years[0], years[-1], color='black', linewidth=0.5)
                # This disables the lower right subfigure. Might need when adding/removing models.
                #axes[7].set_axis_off()
                # Add invisible subplot for common axes labels.
                fig.add_subplot(111, frameon=False)
                plt.tick_params(
                        labelcolor='none',
                        which='both',
                        top=False,
                        bottom=False,
                        left=False,
                        right=False,
                        )
                #plt.plot([0], [0], color=COLORS['CSIRO'], label='ACCESS-ESM1-5') # For legend.
                #for m in MODELS.keys():
                #    plt.plot([0], [0], color=COLORS[m], label=MODELS[m])
                plt.xlabel('Year')
                if var=='tas':
                    plt.ylabel('Temperature (°C)')
                    plt.title('Global mean 2 m surface air temperature difference')
                else: plt.ylabel('Precipitation (mm/day)')
                #fig.suptitle(f'{var} global mean trends')
                #plt.tight_layout()
                #plt.subplots_adjust(left=0.15,) # tight_layout leaves some extra space on the left.
                plt.savefig(f'{PLOTS_DIR}/models/{var}_trends.png', dpi=DPI)
                plt.show()
                plt.close()

            # Plot only the afforestation scenario for absolute values.
            plt.figure()
            # Plot ACCESS ensemble mean and range
            ens_mean = np.mean(access_aff, axis=0)
            years = list(range(2015, 2015 + ens_mean.shape[0]))
            plt.plot(years, ens_mean, color=COLORS['CSIRO'], label='ACCESS-ESM1-5')
            plt.fill_between(
                    years,
                    access_aff.max(axis=0),
                    access_aff.min(axis=0),
                    color='lightblue',
                    )
            # The following plots individual ensemble members as grey lines.
            #for e in range(10):
            #    plt.plot(years, access_aff[e,:], color='gray', alpha=0.4)
            for m in MODELS.keys():
                if m=='MOHC' and var=='cLitter': continue # Skip missing data.
                if m=='BCC' and var=='nbp': continue # Skip missing data.
                if m=='NOAA-GFDL' and var not in ['tas', 'pr']: continue # GFDL only has pr and tas
                years = list(range(2015, 2015 + aff_data[m].shape[0]))
                plt.plot(years, aff_data[m], color=COLORS[m], label=MODELS[m])
            # Plot features.
            plt.xlim(left=years[0], right=years[-1])
            if var=='tas':
                plt.ylabel('Temperature (°C)')
            elif var=='pr':
                plt.ylabel('Precipitation (mm/day)')
            elif var in ['gpp','npp','ra','rh','nbp']:
                plt.ylabel('Pg(C)/year')
            else:
                plt.ylabel('Pg(C)')
            plt.xlabel('Year')
            plt.title(f"{var} esm-ssp585-ssp126Lu")
            plt.legend(frameon=False)
            plt.savefig(
                    f'{PLOTS_DIR}/models/{var}_model_intercomparison_esm-ssp585-ssp126Lu.png',
                    dpi=DPI,
                    )
            # Plot only the reference scenario for absolute values.
            plt.figure()
            # Plot ACCESS ensemble mean and range
            ens_mean = np.mean(access_ssp585, axis=0)
            years = list(range(2015, 2015 + ens_mean.shape[0]))
            plt.plot(years, ens_mean, color=COLORS['CSIRO'], label='ACCESS-ESM1-5')
            plt.fill_between(
                    years,
                    access_ssp585.max(axis=0),
                    access_ssp585.min(axis=0),
                    color='lightblue',
                    )
            # The following plots individual ensemble members as grey lines.
            #for e in range(10):
            #    plt.plot(years, access_ssp585[e,:], color='gray', alpha=0.4)
            for m in MODELS.keys():
                if m=='MOHC' and var=='cLitter': continue # Skip missing data.
                if m=='BCC' and var=='nbp': continue # Skip missing data.
                if m=='NOAA-GFDL' and var not in ['tas', 'pr']: continue # GFDL only has pr and tas
                years = list(range(2015, 2015 + ssp585_data[m].shape[0]))
                plt.plot(years, ssp585_data[m], color=COLORS[m], label=MODELS[m])
            # Plot features.
            plt.xlim(left=years[0], right=years[-1])
            if var=='tas':
                plt.ylabel('Temperature (°C)')
            elif var=='pr':
                plt.ylabel('Precipitation (mm/day)')
            elif var in ['gpp','npp','ra','rh','nbp']:
                plt.ylabel('Pg(C)/year')
            else:
                plt.ylabel('Pg(C)')
            plt.xlabel('Year')
            plt.title(f"{var} esm-ssp585-ssp126Lu")
            plt.legend(frameon=False)
            plt.savefig(
                    f'{PLOTS_DIR}/models/{var}_model_intercomparison_esm-ssp585.png',
                    dpi=DPI,
                    )

    # Plot the soil+litter outside of the variable loops.
    plt.figure()
    # Bias correct CanESM.
    litter_soil_for['CCma'] = litter_soil_for['CCma'] - litter_soil_for['CCma'][0]
    litter_soil_ssp585['CCma'] = litter_soil_ssp585['CCma'] - litter_soil_ssp585['CCma'][0]
    litter_soil_diff = {}
    for m in MODELS.keys():
        if m=='NOAA-GFDL': continue # GFDL only has pr and tas
        # Some models do not have year 2100 in one of the simulations.
        if m=='MPI-M': litter_soil_for[m] = litter_soil_for[m][:-1]
        if m=='NCC': litter_soil_ssp585[m] = litter_soil_ssp585[m][:-1]
        litter_soil_diff[m] = litter_soil_for[m] - litter_soil_ssp585[m]
    # Plot ACCESS first.
    access_litter_soil_diff = access_litter_soil_for - access_litter_soil_ssp585
    access_litter_soil_diff_ensmean = access_litter_soil_diff.mean(axis=0)
    years = list(range(2015, 2015 + ens_mean.shape[0]))
    plt.plot(years, access_litter_soil_diff_ensmean, label='ACCESS-ESM1-5', color=COLORS['CSIRO'])
    plt.fill_between(
            years,
            access_litter_soil_diff.max(axis=0),
            access_litter_soil_diff.min(axis=0),
            color='lightblue',
            )
    # The following plots individual ensemble members as gray lines.
    #for e in range(10):
    #    plt.plot(years, access_litter_soil_diff[e,:], color='gray', alpha=0.4)
    # Plot other models.
    for m in MODELS.keys():
        if m=='NOAA-GFDL': continue # GFDL only has pr and tas
        years = list(range(2015, 2015 + litter_soil_diff[m].shape[0]))
        plt.plot(years, litter_soil_diff[m], color=COLORS[m], label=MODELS[m])
    # Plot features.
    plt.xlim(left=years[0], right=years[-1])
    plt.hlines(0, years[0], years[-1], color='black', linewidth=0.5)
    plt.ylabel('Pg(C)')
    plt.xlabel('Year')
    plt.title(f"cSoil+cLitter for esm-ssp585-ssp126Lu - esm-ssp585")
    plt.legend(frameon=False)
    plt.savefig(f'{PLOTS_DIR}/models/cSoil+cLitter_model_intercomparison_diff.png', dpi=DPI)
    plt.close()

if __name__ != 'analysis.plot_models':
    make_model_plots()

    # Clean up
    temp_files = glob.glob('./cdoPy*')
    for f in temp_files:
        os.remove(f)

    plt.show()


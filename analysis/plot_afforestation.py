#!/usr/bin/env python3
# plot_afforestation.py generates plots for the afforestation simulation.
# In general "aff" is used to refer to the afforestation simulation esm-ssp585-ssp126Lu (high
# fossil fuels emissions and low land use change scenario), and "ssp585" means the emissions driven
# esm-ssp585 simulation (high fossil fuel emisssions scenario).
import glob

import cdo as cdo_module
import matplotlib.pyplot as plt
import numpy as np
from analysis.baseline import global_sum_baselines
from analysis.cmip_files import get_filename, LAND_FRAC_FILE
from analysis.constants import TABLES, VARIABLES, ENSEMBLES, SEC_IN_DAY, KG_IN_PG

cdo = cdo_module.Cdo(tempdir='.')
#cdo.debug = True

# Local constants
COLORS = {'gpp':'green',
        'npp':'olive',
        'ra':'orange',
        'rh':'saddlebrown',
        'nbp':'purple',
        'cVeg':'darkgreen',
        'cLitter':'chocolate',
        'cSoil':'black'}
PLOTS_DIR = 'plots'


# I feel like I should be using decorators here to construct the cdo commands.
def cdo_flux_aggregate_convert_load(files: list):
    """Uses CDO to globally agggregate a carbon flux variableand load it into an array.
    files - List of stirings of filenames containing the variables.
    TODO: Currently only using the first file in the list. Need to loop if there are many.
    """
    return cdo.divc(str(KG_IN_PG),
            input=f'-yearsum -muldpm -mulc,{SEC_IN_DAY} -fldsum -mul -mul '+files[0]+ \
                    f' -divc,100 {LAND_FRAC_FILE} -gridarea {LAND_FRAC_FILE}',
            options='-L',
            returnCdf=True).variables[var][:].squeeze()


def cdo_pool_aggregate_convert_load(files: list):
    """Uses CDO to globally agggregate a carbon pool variableand load it into an array.
    files - List of stirings of filenames containing the variables.
    TODO: Currently only using the first file in the list. Need to loop if there are many.
    """
    return cdo.divc(str(KG_IN_PG),
            input=f'-yearmonmean -fldsum -mul -mul '+files[0]+ \
                    f' -divc,100 {LAND_FRAC_FILE} -gridarea {LAND_FRAC_FILE}',
            options='-L',
            returnCdf=True).variables[var][:].squeeze()


for table in TABLES:
    for var in VARIABLES[table]:
        # Initialize data arrays 
        aff_data = np.ones((len(ENSEMBLES), 86))*np.nan
        ssp585_data = np.ones((len(ENSEMBLES), 86))*np.nan
        for ens in ENSEMBLES:
            # Get file names
            aff_files = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, table, var)
            ssp585_files = get_filename('C4MIP', 'esm-ssp585', ens, table, var)
            if var in ['gpp', 'npp', 'ra', 'rh', 'nbp']:
                # Load flux variable
                aff_data[int(ens)-1,:] = cdo_flux_aggregate_convert_load(aff_files)
                ssp585_data[int(ens)-1,:] = cdo_flux_aggregate_convert_load(ssp585_files)
            else:
                # Load pool variable
                aff_data[int(ens)-1,:] = cdo_pool_aggregate_convert_load(aff_files)
                ssp585_data[int(ens)-1,:] = cdo_pool_aggregate_convert_load(ssp585_files)

        # Anomaly, mean and standard deviation relative to 2015 baseline. Demonstrates overall
        # impact of climate, CO2 forcing and afforestation on pool/flux.
        data_anomaly = aff_data - global_sum_baselines[var]
        data_ens_mean = np.mean(data_anomaly, axis=0)
        data_ens_std = np.std(data_anomaly, axis=0, ddof=1)

        # Anomaly relative to the esm-ssp585 scenario from C4MIP. Domenstrates only the changes
        # due to afforestation.
        data_aff_diff = ssp585_data - aff_data
        data_aff_diff_mean = np.mean(data_aff_diff, axis=0)
        data_aff_diff_std = np.std(data_aff_diff, axis=0, ddof=1)
        
        # Plot the graphs for 
        years = [y for y in range(2015, 2101)]
        plt.figure()
        for ens in ENSEMBLES:
            plt.plot(years, data_anomaly[int(ens)-1], color='lightgray', linewidth=0.6)
        plt.plot(years, data_ens_mean, color=COLORS[var], label="Ensemble mean")
        plt.plot(years, data_ens_mean+data_ens_std, color=COLORS[var], linewidth=0.8,
                label="+-1$\sigma$")
        plt.plot(years, data_ens_mean-data_ens_std, color=COLORS[var], linewidth=0.8)
        plt.hlines(0, years[0], years[-1], colors='black', linestyles='dotted')
        plt.xlabel('Year')
        plt.ylabel(f'{var.upper()} anomaly (PgC/year)')
        plt.title(f"ACCESS-ESM1-5 {var.upper()}")
        plt.savefig(f'{PLOTS_DIR}/{var}_ACCESS-ESM1-5_esm-ssp585-ssp126Lu_ensembles_anomalies.svg')

# Clean up
temp_files = glob.glob('./cdoPy*')
os.remove(temp_files)

#plt.show()


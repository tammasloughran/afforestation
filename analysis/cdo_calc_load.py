# cdo_calc_load.py uses CDO to calculate global aggregations and unit conversions, then loads the
# data into numpy arrays.
import cdo as cdo_module
import numpy as np

if __name__ != 'analysis.cdo_calc_load':
    from cmip_files import LAND_FRAC_FILE, get_filename
    from constants import ENSEMBLES, KG_IN_PG, SEC_IN_DAY, TABLES, VARIABLES
else:
    # cdo_calc_load.py imported as a module of the analysis package.
    from analysis.cmip_files import LAND_FRAC_FILE, get_filename
    from analysis.constants import (ENSEMBLES, KG_IN_PG, SEC_IN_DAY, TABLES,
                                    VARIABLES)


cdo = cdo_module.Cdo(tempdir='.')
#cdo.debug = True

# I feel like I should be using decorators here to construct the cdo commands.
def cdo_flux_aggregate_convert_load(files: list, var: str):
    """Uses CDO to globally agggregate a carbon flux variableand load it into an array.
    files - List of stirings of filenames containing the variables.
    TODO: Currently only using the first file in the list. Need to loop if there are many.
    """
    return cdo.divc(str(KG_IN_PG),
            input=f'-yearsum -muldpm -mulc,{SEC_IN_DAY} -fldsum -mul -mul '+files[0]+ \
                    f' -divc,100 {LAND_FRAC_FILE} -gridarea {LAND_FRAC_FILE}',
            options='-L',
            returnCdf=True).variables[var][:].squeeze()


def cdo_pool_aggregate_convert_load(files: list, var: str):
    """Uses CDO to globally agggregate a carbon pool variableand load it into an array.
    files - List of stirings of filenames containing the variables.
    TODO: Currently only using the first file in the list. Need to loop if there are many.
    """
    return cdo.divc(str(KG_IN_PG),
            input=f'-yearmonmean -fldsum -mul -mul '+files[0]+ \
                    f' -divc,100 {LAND_FRAC_FILE} -gridarea {LAND_FRAC_FILE}',
            options='-L',
            returnCdf=True).variables[var][:].squeeze()


def cdo_load(cdo_string: str):
    """Load a file using CDO.
    """
    return cdo.copy(input=cdo_string, options='-L', returnCdf=True).variables[var].squeeze()


def cdo_fetch_ensembles(mip, exp, table, var):
    """Loads all ensemble data using cdo functions.
    """
    # Initialize data arrays 
    exp_data = np.ones((len(ENSEMBLES), 86))*np.nan
    for ens in ENSEMBLES:
        # Get file names
        exp_files = get_filename(mip, exp, ens, table, var)
        if var in ['gpp', 'npp', 'ra', 'rh', 'nbp']:
            # Load flux variable
            exp_data[int(ens)-1,:] = cdo_flux_aggregate_convert_load(exp_files, var)
        else:
            # Load pool variable
            exp_data[int(ens)-1,:] = cdo_pool_aggregate_convert_load(exp_files, var)
    return exp_data


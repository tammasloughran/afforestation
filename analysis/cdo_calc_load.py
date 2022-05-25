# cdo_calc_load.py uses CDO to calculate global aggregations and unit conversions, then loads the
# data into numpy arrays.
import functools
import cdo as cdo_module
import numpy as np

if __name__ != 'analysis.cdo_calc_load':
    # cdo_calc_load.py imported as a module from other scripts.
    from cmip_files import LAND_FRAC_FILE, get_filename
    from constants import ENSEMBLES, KG_IN_PG, SEC_IN_DAY, TABLES, VARIABLES, CLIM_VARIABLES
else:
    # cdo_calc_load.py imported as a module of the analysis package.
    from analysis.cmip_files import LAND_FRAC_FILE, get_filename
    from analysis.constants import (ENSEMBLES, KG_IN_PG, SEC_IN_DAY, TABLES,
                                    VARIABLES, CLIM_VARIABLES)

cdo = cdo_module.Cdo(tempdir='.')
cdo.debug = False


def cdo_div_100(cdo_func):
    """Wrapper to add division by 100 to the CDO command. I.e. converts percentage into a fraction.
    """
    @functools.wraps(cdo_func)
    def wrapper_div_100(cdo_string:str, var:str):
        cdo_string = '-divc,100 ' + cdo_string
        return cdo_func(cdo_string, var)
    return wrapper_div_100


def cdo_fldmean(cdo_func):
    """Wrapper to add fldmean to the CDO command.
    """
    @functools.wraps(cdo_func)
    def wrapper_fldmean(cdo_string:str, var:str):
        cdo_string = '-fldmean ' + cdo_string
        return cdo_func(cdo_string, var)
    return wrapper_fldmean


def cdo_fldsum(cdo_func):
    """Wrapper to add fldsum to the CDO command.
    """
    @functools.wraps(cdo_func)
    def wrapper_fldsum(cdo_string:str, var:str):
        cdo_string = '-fldsum ' + cdo_string
        return cdo_func(cdo_string, var)
    return wrapper_fldsum


def cdo_kg2pg(cdo_func):
    """Wrapper to add conversion factors from kg to Pg to the CDO command.
    """
    @functools.wraps(cdo_func)
    def wrapper_kg2pg(cdo_string:str, var:str):
        cdo_string = f'-divc,{KG_IN_PG} ' + cdo_string
        return cdo_func(cdo_string, var)
    return wrapper_kg2pg


def cdo_mul_land_area(cdo_func):
    """Wrapper to add multiplication of the land area and land cover fractions to the CDO command.
    """
    @functools.wraps(cdo_func)
    def wrapper_mul_land_area(cdo_string:str, var:str):
        cdo_string = '-mul -mul ' + cdo_string + \
                f' -divc,100 {LAND_FRAC_FILE} -gridarea {LAND_FRAC_FILE}'
        return cdo_func(cdo_string, var)
    return wrapper_mul_land_area


def cdo_mul_grid_area(cdo_func):
    """Wrapper to add multiplication of the grid area to the CDO command.
    """
    @functools.wraps(cdo_func)
    def wrapper_mul_grid_area(cdo_string:str, var:str):
        cdo_string = '-mul ' + cdo_string + f' -gridarea {LAND_FRAC_FILE}'
        return cdo_func(cdo_string, var)
    return wrapper_mul_grid_area


def cdo_sec2mon(cdo_func):
    """Wrapper to add conversion factor from per seconds to per month to the CDO command.
    """
    @functools.wraps(cdo_func)
    def wrapper_sec2mon(cdo_string:str, var:str):
        cdo_string = f'-muldpm -mulc,{SEC_IN_DAY} ' + cdo_string
        return cdo_func(cdo_string, var)
    return wrapper_sec2mon


def cdo_yearmonmean(cdo_func):
    """Wrapper to add yearmonmean to the CDO command.
    """
    @functools.wraps(cdo_func)
    def wrapper_yearmonmean(cdo_string:str, var:str):
        cdo_string = '-yearmonmean ' + cdo_string
        return cdo_func(cdo_string, var)
    return wrapper_yearmonmean


def cdo_yearsum(cdo_func):
    """Wrapper to add yearsum to the CDO command.
    """
    @functools.wraps(cdo_func)
    def wrapper_yearsum(cdo_string:str, var:str):
        cdo_string = '-yearsum ' + cdo_string
        return cdo_func(cdo_string, var)
    return wrapper_yearsum


@cdo_fldmean
@cdo_yearmonmean
def cdo_clim_load(cdo_file:str, var:str):
    """Load global climate variable using CDO. Please refer to the following decorators:
    @cdo_fldmean
    @cdo_yearmonmean
    """
    return cdo.copy(input=cdo_file, options='-L', returnCdf=True).variables[var][:].squeeze()


@cdo_div_100
@cdo_mul_grid_area
@cdo_fldsum
@cdo_yearmonmean
def cdo_cover_area_load(cdo_file:str, var:str):
    """Load a land cover fraction variable and express as an area using CDO. Please refer to the
    following decorators:
    @cdo_div_100
    @cdo_mul_grid_area
    @cdo_fldsum
    """
    return cdo.copy(input=cdo_file, options='-L', returnCdf=True).variables[var][:].squeeze()


@cdo_mul_land_area
@cdo_fldsum
@cdo_sec2mon
@cdo_yearsum
@cdo_kg2pg
def cdo_flux_load(cdo_file:str, var:str):
    """Load a flux variable using CDO. Please refer to the following decorators to see the
    aggregation and units conversion:
    @cdo_mul_land_area
    @cdo_fldsum
    @cdo_sec2year
    @cdo_yearsum
    @cdo_kg2pg
    """
    return cdo.copy(input=cdo_file, options='-L', returnCdf=True).variables[var][:].squeeze()


@cdo_mul_land_area
@cdo_fldsum
@cdo_yearmonmean
@cdo_kg2pg
def cdo_pool_load(cdo_file:str, var:str):
    """Load a pool variable using CDO. Please refer to the following decorators to see the
    aggregation and units conversion:
    @cdo_mul_land_area
    @cdo_fldsum
    @cdo_yearmonmean
    @cdo_kg2pg
    """
    return cdo.copy(input=cdo_file, options='-L', returnCdf=True).variables[var][:].squeeze()


def cdo_fetch_ensembles(mip:str, exp:str, table:str, var:str):
    """Loads all ensemble data using cdo functions.
    """
    # Initialize data arrays
    exp_data = np.ones((len(ENSEMBLES), 86))*np.nan
    for ens in ENSEMBLES:
        # Get file names
        exp_files = get_filename(mip, exp, ens, table, var)
        if var in ['gpp', 'npp', 'ra', 'rh', 'nbp']:
            # Load flux variable
            exp_data[int(ens)-1,:] = cdo_flux_load(exp_files[0], var)
        elif var in ['cVeg', 'cLitter', 'cSoil']:
            # Load pool variable
            exp_data[int(ens)-1,:] = cdo_pool_load(exp_files[0], var)
        elif var in CLIM_VARIABLES[table]:
            exp_data[int(ens)-1,:] = cdo_clim_load(exp_files[0], var)
    return exp_data


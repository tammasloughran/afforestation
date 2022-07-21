# cdo_calc_load.py uses CDO to calculate global aggregations and unit conversions, then loads the
# data into numpy arrays.
import functools
from types import FunctionType

import cdo as cdo_module
import cdo_decorators as cdod
import numpy as np

if __name__ != 'analysis.cdo_calc_load':
    # cdo_calc_load.py imported as a module from other scripts.
    from cmip_files import LAND_FRAC_FILE, get_filename
    from constants import (CLIM_VARIABLES, ENSEMBLES, KG_IN_PG, SEC_IN_DAY,
                           SEC_IN_YEAR, TABLES, VARIABLES)
else:
    # cdo_calc_load.py imported as a module of the analysis package.
    from analysis.cmip_files import LAND_FRAC_FILE, get_filename
    from analysis.constants import (CLIM_VARIABLES, ENSEMBLES, KG_IN_PG,
                                    SEC_IN_DAY, SEC_IN_YEAR, TABLES, VARIABLES)

cdo = cdo_module.Cdo(tempdir='.')
cdo.debug = False


def cdo_mul_land_area(cdo_func):
    """Decorator to add multiplication of the land area and land cover fractions to the CDO command.
    """
    @functools.wraps(cdo_func)
    def wrapper_mul_land_area(*args, **kwargs):
        kwargs['input'] = f'-mul -mul {kwargs["input"]} -divc,100 {LAND_FRAC_FILE} '\
                f'-gridarea {LAND_FRAC_FILE}'
        return cdo_func(*args, **kwargs)
    return wrapper_mul_land_area


#@cdod.cdo_selregion('dcw:AU') # I think this command is only in a newer verion of CDO.
@cdod.cdo_cat(input2='') # Concatenate all files in input1.
@cdo_mul_land_area
@cdod.cdo_sellonlatbox(str(aus_east), str(aus_west), str(aus_south), str(aus_north))
@cdod.cdo_fldsum
@cdod.cdo_yearmonmean
@cdod.cdo_divc(str(KG_IN_PG))
def load_aus_pool(var:str, input:str)->np.ma.MaskedArray:
    """Load pool variable using CDO. Please refer to the following decorators to see the
    aggregation and units conversion:
    @cdod.cdo_cat(input2='')
    @cdo_mul_land_area
    @cdod.cdo_sellonlatbox(str(aus_east), str(aus_west), str(aus_south), str(aus_north))
    @cdod.cdo_fldsum
    @cdod.cdo_yearmonmean
    @cdod.cdo_divc(str(KG_IN_PG))
    """
    return cdo.copy(input=input, returnCdf=True, options='-L').variables[var][:].squeeze()


@cdod.cdo_cat(input2='')
@cdod.cdo_ifthen(input1=LAND_FRAC_FILE) # Mask for climate over land only.
@cdod.cdo_sellonlatbox(str(aus_east), str(aus_west), str(aus_south), str(aus_north))
@cdod.cdo_fldmean
@cdod.cdo_yearmonmean
def load_aus_clim(var:str, input:str)->np.ma.MaskedArray:
    """Load global climate variable using CDO. Please refer to the following decorators:
    @cdod.cdo_cat(input2='')
    @cdod.cdo_ifthen(input1=LAND_FRAC_FILE) # Mask for climate over land only.
    @cdod.cdo_sellonlatbox(str(aus_east), str(aus_west), str(aus_south), str(aus_north))
    @cdod.cdo_fldmean
    @cdod.cdo_yearmonmean
    """
    return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


@cdod.cdo_cat(input2='')
@cdo_mul_land_area
@cdod.cdo_sellonlatbox(str(aus_east), str(aus_west), str(aus_south), str(aus_north))
@cdod.cdo_fldsum
@cdod.cdo_mulc(str(SEC_IN_DAY))
@cdod.cdo_muldpm
@cdod.cdo_yearsum
@cdod.cdo_divc(str(KG_IN_PG))
def load_aus_flux(var:str, input:str)->np.ma.MaskedArray:
    """Load flux variable using CDO. Please refer to the following decorators to see the
    aggregation and units conversion:
    @cdod.cdo_cat(input2='')
    @cdo_mul_land_area
    @cdod.cdo_sellonlatbox(str(aus_east), str(aus_west), str(aus_south), str(aus_north))
    @cdod.cdo_fldsum
    @cdod.cdo_mulc(str(SEC_IN_DAY))
    @cdod.cdo_muldpm
    @cdod.cdo_yearsum
    @cdod.cdo_divc(str(KG_IN_PG))
    """
    return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


@cdod.cdo_cat(input2='')
@cdo_mul_land_area
@cdod.cdo_sellonlatbox(str(aus_east), str(aus_west), str(aus_south), str(aus_north))
@cdod.cdo_fldsum
@cdod.cdo_mulc(str(SEC_IN_YEAR))
@cdod.cdo_divc(str(KG_IN_PG))
def load_aus_base_flux(var:str, input:str)->np.ma.MaskedArray:
    """Load baseline flux data. Similar to load_aus_flux but without the temporal aggregation.
    Temporal aggregation has already been done in `analysis.baseline`. Please refer to the 
    following decorators to see the aggregation and units conversion:
    @cdod.cdo_cat(input2='')
    @cdo_mul_land_area
    @cdod.cdo_sellonlatbox(str(aus_east), str(aus_west), str(aus_south), str(aus_north))
    @cdod.cdo_fldsum
    @cdod.cdo_mulc(str(SEC_IN_YEAR))
    @cdod.cdo_divc(str(KG_IN_PG))
    """
    return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


## Decorator to add multiplication of the grid area to the CDO command
cdo_mul_grid_area = cdod.cdo_mul(input2=cdod.get_str(cdod.cdo_gridarea, input=LAND_FRAC_FILE))


def cdo_sec2mon(cdo_func:FunctionType)->np.ma.MaskedArray:
    """Decorator to add conversion factor from per seconds to per month to the CDO command.
    """
    @functools.wraps(cdo_func)
    def wrapper_sec2mon(*args, **kwargs):
        mwargs['input'] = f'-muldpm -mulc,{SEC_IN_DAY} kwargs["input"]'
        return cdo_func(*args, **kwargs)
    return wrapper_sec2mon


@cdod.cdo_cat(input2='')
@cdod.cdo_ifthen(input1=LAND_FRAC_FILE) # Mask for climate over land only.
@cdod.cdo_fldmean
@cdod.cdo_yearmonmean
def cdo_clim_load(var:str, input:str)->np.ma.MaskedArray:
    """Load global climate variable using CDO. Please refer to the following decorators:
    @cdod.cdo_cat(input2='')
    @cdod.cdo_ifthen(input1=LAND_FRAC_FILE) # Mask for climate over land only.
    @cdod.cdo_fldmean
    @cdod.cdo_yearmonmean
    """
    return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


@cdod.cdo_divc('100')
@cdo_mul_grid_area
@cdod.cdo_fldsum
@cdod.cdo_yearmonmean
def cdo_cover_area_load(var:str, input:str)->np.ma.MaskedArray:
    """Load a land cover fraction variable and express as an area using CDO. Please refer to the
    following decorators:
    @cdod.cdo_divc('100')
    @cdo_mul_grid_area
    @cdod.cdo_fldsum
    @cdod.cdo_yearmonmean
    """
    return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


@cdod.cdo_yearmonmean
@cdod.cdo_selyear('2015','2100')
@cdo_mul_grid_area
@cdod.cdo_deltat
def cdo_area_diff_load(var:str, input:str)->np.ma.MaskedArray:
    """Load themap of the difference between 2015 and 2100 in land cover area.
    @cdo.cdo_selyear(str(2015),str(2100))
    @cdo.cdo_deltat
    """
    return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


@cdod.cdo_cat(input2='')
@cdo_mul_land_area
@cdod.cdo_fldsum
@cdod.cdo_mulc(str(SEC_IN_DAY))
@cdod.cdo_muldpm
@cdod.cdo_yearsum
@cdod.cdo_divc(str(KG_IN_PG))
def cdo_flux_load(var:str, input:str)->np.ma.MaskedArray:
    """Load a flux variable using CDO. Please refer to the following decorators to see the
    aggregation and units conversion:
    @cdod.cdo_cat(input2='')
    @cdo_mul_land_area
    @cdod.cdo_fldsum
    @cdod.cdo_mulc(str(SEC_IN_DAY))
    @cdod.cdo_muldpm
    @cdod.cdo_yearsum
    @cdod.cdo_divc(str(KG_IN_PG))
    """
    return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


@cdod.cdo_cat(input2='')
@cdo_mul_land_area
@cdod.cdo_fldsum
@cdod.cdo_yearmonmean
@cdod.cdo_divc(str(KG_IN_PG))
def cdo_pool_load(var:str, input:str)->np.ma.MaskedArray:
    """Load a pool variable using CDO. Please refer to the following decorators to see the
    aggregation and units conversion:
    @cdod.cdo_cat(input2='')
    @cdo_mul_land_area
    @cdod.cdo_fldsum
    @cdod.cdo_yearmonmean
    @cdod.cdo_divc(str(KG_IN_PG))
    """
    return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


def cdo_fetch_ensembles(mip:str, exp:str, table:str, var:str)->np.ndarray:
    """Loads all ensemble data using cdo functions.
    """
    # Initialize data arrays
    exp_data = np.ones((len(ENSEMBLES), 86))*np.nan
    for ens in ENSEMBLES:
        # Get file names
        exp_files = ' '.join(get_filename(mip, exp, ens, table, var))
        exp_files = '[ '+exp_files+' ]'
        if var in ['gpp', 'npp', 'ra', 'rh', 'nbp']:
            # Load flux variable
            exp_data[int(ens)-1,:] = cdo_flux_load(var, input=exp_files)
        elif var in ['cVeg', 'cLitter', 'cSoil']:
            # Load pool variable
            exp_data[int(ens)-1,:] = cdo_pool_load(var, input=exp_files)
        elif var in CLIM_VARIABLES[table]:
            exp_data[int(ens)-1,:] = cdo_clim_load(var, input=exp_files)
    return exp_data


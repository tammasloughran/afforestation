# cmip_files.py contains variables and functions relating to CMIP6 files on Gadi.
import os
import glob

LAND_FRAC_FILE = '/g/data/fs38/publications/CMIP6/LUMIP/CSIRO/ACCESS-ESM1-5/esm-ssp585-ssp126Lu/'\
        'r10i1p1f1/fx/sftlf/gn/latest/sftlf_fx_ACCESS-ESM1-5_esm-ssp585-ssp126Lu_r10i1p1f1_gn.nc'
SSP585126LU_ENS = [
        'SSP-585-126-03',
        'SSP-585-126-04',
        'SSP-585-126-06',
        'SSP-585-126-07',
        'SSP-585-126-08',
        'SSP-585-126-09',
        'SSP-585-126-10',
        'SSP-585-126-11',
        'SSP-585-126-12',
        'SSP-585-126-13',
        ]


def get_filename(mip:str, exp:str, ens:str, table:str, var:str):
    """Insert MIP experiment, ensemble, table, and variable strings into a CMOR directory tree
    and filename string. Automatically detects the period for the file name.
    """
    # The ensemble field is specific to ACCESS. There may be other models that use a different
    # ensemble type.
    directory = f'/g/data/fs38/publications/CMIP6/{mip}/CSIRO/ACCESS-ESM1-5/{exp}/r{ens}i1p1f1/'\
            f'{table}/{var}/gn/latest/'
    files = os.listdir(directory)
    return [directory+f for f in files]


def get_filenames(mip:str, inst:str, model:str, exp:str, ens:str, table:str, var:str):
    """Insert MIP, model, experiment, ensemble, table, and variable strings into a CMOR directory
    tree and filename string. Automatically detects the period for the file name.
    """
    directory = f'/g/data/p66/tfl561/CMIP6/{mip}/{inst}/{model}/{exp}/{ens}/{table}/{var}'\
            '/gn/latest/'
    files = os.listdir(directory)
    files.sort()
    return [directory+f for f in files]


def get_archive_filename(exp:str, freq:str='mon')->str:
    """Insert experiment string into archive directory tree and return all files of the
    specified frequency.
    exp : experiment to get.
    freq : 'mon' or 'dai' return monthly or daily files.
    """
    if freq=='mon':
        file_freq = '*.pa-*_mon.nc'
    elif freq=='dai':
        file_freq = '*.pe-*_dai.nc'
    else:
        raise ValueError("Unknown frequency")
    directory = f'/g/data/p73/archive/CMIP6/ACCESS-ESM1-5/{exp}/history/atm/netCDF/'
    return glob.glob(directory+file_freq)


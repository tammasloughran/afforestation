# cmip_files.py contains variables and functions relating to CMIP6 files on Gadi.
import os

LAND_FRAC_FILE = '/g/data/fs38/publications/CMIP6/LUMIP/CSIRO/ACCESS-ESM1-5/esm-ssp585-ssp126Lu/'\
        'r10i1p1f1/fx/sftlf/gn/latest/sftlf_fx_ACCESS-ESM1-5_esm-ssp585-ssp126Lu_r10i1p1f1_gn.nc'


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


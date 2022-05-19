# Define commonly used constants here
CLIM_VARIABLES = {
        'Amon': ['tas', 'pr'],}
ENSEMBLES = [str(e) for e in range(1,11)]
EXPERIMENTS = {
        'CMIP': ['esm-hist'],
        'C4MIP': ['esm-ssp585'],
        'LUMIP': ['esm-ssp585-ssp126Lu']}
FRAC_VARIABLES = {'Lmon': ['cropFrac',
        'treeFrac',
        #'shrubFrac', # There isn't much difference in shrub between experiments.
        'grassFrac']}
KG_IN_PG = 1000000000000
M2_IN_MILKM2 = 10**12
MIPS = ['CMIP', 'LUMIP']
SEC_IN_DAY = 60*60*24
SEC_IN_YEAR = 60*60*24*365
NEW_UNITS_FACTOR = SEC_IN_YEAR/KG_IN_PG
TABLES = ['Lmon', 'Emon']
VARIABLES = {
        'Lmon': ['gpp', 'npp', 'ra', 'rh', 'nbp', 'cVeg', 'cLitter'],
        'Emon': ['cSoil']}


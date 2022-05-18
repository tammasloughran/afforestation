# Define commonly used constants here
MIPS = ['CMIP', 'LUMIP']
EXPERIMENTS = {
        'CMIP': ['esm-hist'],
        'C4MIP': ['esm-ssp585'],
        'LUMIP': ['esm-ssp585-ssp126Lu']}
TABLES = ['Lmon', 'Emon']
VARIABLES = {
        'Lmon': ['gpp', 'npp', 'ra', 'rh', 'nbp', 'cVeg', 'cLitter'],
        # ACCESS does not have nep
        'Emon': ['cSoil']}
FRAC_VARIABLES = {
        'Lmon': ['cropFrac',
            'treeFrac',
            #'shrubFrac', # There isn't much difference in shrub between experiments.
            'grassFrac']}
ENSEMBLES = [str(e) for e in range(1,11)]
KG_IN_PG = 1000000000000
SEC_IN_YEAR = 60*60*24*365
SEC_IN_DAY = 60*60*24
NEW_UNITS_FACTOR = SEC_IN_YEAR/KG_IN_PG


# Define commonly used constants here
CLIM_VARIABLES = {'Amon':['tas','pr'],'Emon':[]}
DATA_DIR = './data'
ENSEMBLES = [str(e) for e in range(1,11)]
EXPERIMENTS = {
        'CMIP':['esm-hist'],
        'C4MIP':['esm-ssp585'],
        'LUMIP':['esm-ssp585-ssp126Lu'],
        }
FRAC_VARIABLES = {
        'Lmon':[
                'cropFrac',
                'treeFrac',
                #'shrubFrac', # There isn't much difference in shrub between experiments.
                'grassFrac',
                ],
        }
KG_IN_PG = 1000000000000 # kg
M2_IN_MILKM2 = 10**12 # m^-2
MASS_ATM = 5.27E18 # kg from Cook, A.H. Physics of the Earth and Planets. (1973), 276.
MIL = 1000000
MIPS = ['CMIP','LUMIP']
MOLMASS_AIR = 0.0289652 # kg/mol
MOLMASS_CO2 = 0.0440095 # kg/mol
MOLMASS_O2 = 0.0319988 # kg/mol
NENS = len(ENSEMBLES)
NTIMES = 2101 - 2015 # years
PLOTS_DIR = './plots'
SEC_IN_DAY = 60*60*24 # s
SEC_IN_YEAR = 60*60*24*365 # s
TABLES = ['Lmon','Emon']
VARIABLES = {
        'Lmon':[
            'gpp',
            'npp',
            'ra',
            'rh',
            'nbp',
            'cVeg',
            'cLitter',
            ],
        'Emon':[
            'cSoil',
            'cLand',
            ],
        }

# Depends on constants already defined.
NEW_UNITS_FACTOR = SEC_IN_YEAR/KG_IN_PG # s/kg
KGKG_TO_MOLMOL = MOLMASS_AIR/MOLMASS_CO2 # Converion of kg/kg to mol/mol


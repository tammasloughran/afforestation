#!/usr/bin/env python3
# baseline.py calculate a baseline value for several variables over a base period.
# For historical to scenario simulations, the base period should be a 20 year mean centred at the
# start of the future simulation.
import os
import glob
import cdo as cdo_module
from cmip_files import LAND_FRAC_FILE, get_filename

cdo = cdo_module.Cdo(tempdir='.')

# Constants.
MIPS = ['CMIP', 'LUMIP']
EXPERIMENTS = {
        'CMIP': ['esm-hist'],
        'LUMIP': ['esm-ssp585-ssp126Lu']}
TABLES = ['Lmon', 'Emon']
VARIABLES = {
        'Lmon': ['gpp', 'npp', 'ra', 'rh', 'nbp', 'cVeg', 'cLitter'],
        # ACCESS does not have nep
        'Emon': ['cSoil']}
ENSEMBLES = [str(e) for e in range(1,11)]
KG_IN_PG = 1000000000000
SEC_IN_YEAR = 60*60*24*365
NEW_UNITS_FACTOR = SEC_IN_YEAR/KG_IN_PG


def reference_period(infile1, infile2, outfile, pyear=[2005, 2024]):
    """Extract a map of the mean over the reference period using CDO.
    The reference period spans 20 years centred on the start of the future simulation (2015), so
    [2005, 2024].
    """
    cdo.timmean(input=f'-selyear,{pyear[0]}/{pyear[1]} -cat '+infile1+' '+infile2, output=outfile)


# Options.
## Recalculate the ensemble means from the CMORized CMIP6 data. If not, assume it's been done.
recalculate_ens_mean = False

# Local variables.
global_sum_baselines = {}

# Calculate the maps of base period means for each variable and ensemble member.
for table in TABLES:
    for var in VARIABLES[table]:
        if recalculate_ens_mean:
            for ens in ENSEMBLES:
                esmhist_files = get_filename('CMIP', 'esm-hist', ens, table, var)
                aff_files = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, table, var)
                output_file = f'{var}_ACCESS-ESM1-5_esm-hist-aff_r{ens}i1p1f1_200501-202412mean.nc'
                reference_period(esmhist_files[-1], aff_files[0], output_file)
            # The analysis for all ensemble members should be with respect to the same baseline
            # value. So I calculate the ensemble mean here.
            cdo.ensmean(input=f'{var}_ACCESS-ESM1-5_esm-hist-aff_r*i1p1f1_200501-202412mean.nc',
                    output=f'{var}_ACCESS-ESM1-5_esm-hist-aff_ensmean_200501-202412mean.nc')
            # Clean up unneeded ensemble files.
            efiles = glob.glob(f'{var}_ACCESS-ESM1-5_esm-hist-aff_r*i1p1f1_200501-202412mean.nc')
            for f in efiles:
                os.remove(f)
        # Calculate the global sum of each variable's ensemble mean. If a flux, convert to /year.
        # Need to also multiply by the grid cell area and land fraction.
        # LAND_FRAC_FILE is in % so must be divided by 100 first.
        # kg m-2 s-2 -> *land_frac*landarea*SEC_IN_YEAR/KG_IN_PG -> Pg/year
        if var in ['cVeg', 'cLitter', 'cSoil']:
            time_units = 1
        else:
            time_units = SEC_IN_YEAR
        global_sum_baselines[var] = cdo.divc(str(KG_IN_PG),
                input=f'-mulc,{time_units} -fldsum -mul -mul '\
                        f'{var}_ACCESS-ESM1-5_esm-hist-aff_ensmean_200501-202412mean.nc '\
                        f'-divc,100 {LAND_FRAC_FILE} -gridarea {LAND_FRAC_FILE}',
                options='-L',
                returnCdf=True).variables[var][:][0,0,0]

print("Global sum baseline vales:")
print(global_sum_baselines)


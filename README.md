# Afforestation

# Project structure
```
.
├── analysis
│   ├── baseline.py             Calculate baseline values over a base period.
│   ├── cdo_calc_load.py        CDO calculations and data loading.
│   ├── cmip_files.py           Locations and function for CMIP files.
│   ├── constants.py            Constants definitions.
│   ├── __main__.py             Main entry point program.
│   └── plot_afforestation.py   Plotting script.
├── data                        Location of output data files.
│   └── README.md
├── doc
├── LICENCE                     CSIRO Non-Commercial Source Code Licence Agreement v1.0
├── plots                       Location of output plots.
│   └── README.md
├── pyproject.toml
├── README.md
├── setup.cfg
└── tags
```

## Setup and run
The following install the analysis scripts as a package in a virtual environment and runs the
analysis.
```bash
virtualenv affenv
source ./affenv/bin/activate
pip install .
python -m analysis
```

## Generate ctags
    ctags --recurse --python-kinds=-i analysis

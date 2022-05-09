# Afforestation

# Project structure
```
.
├── analysis
│   ├── baseline.py             Calculate baseline values over a base period.
│   ├── cmip_files.py           Locations and function for CMIP files.
│   ├── constants.py            Constants
│   ├── __main__.py
│   └── plot_afforestation.py   Plotting script.
├── data                        Location of output data files.
├── plots                       Location of output plots.
├── README.md
├── pyproject.toml
└── setup.cfg
```

## Setup and run
```bash
virtualenv affenv
source ./affenv/bin/activate
pip install .
python -m analysis
```

## Generate ctags:
    ctags --recursive --python-kinds=-i *.py

# ECHAM6 standalone radiation code 
This repository contains code to run the ECHAM6 radiative transfer offline. It mostly consists of Fortran code that has been copied from the model source code and has only been slightly adapted to run without the host model. There are two options to interface the radiative transfer code, one directly via a Fortran and one via a Python wrapper into the Fortran source code. It nevertheless is strongly recommended to use the Python wrapper.

## Installation
The FORTRAN code has to be compiled with the netCDF and the MPICH library. They should be present on some machines (e.g. on LIM/DKRZ servers). If not, those libraries can be installed via conda.

### Using conda
If you use conda, first set up an environment:

```bash
conda env create -f environment.yml
conda activate echam_radiation_offline
```

After installing the libraries, you need to give the path to your newly created environment in `build/Makeinclude_fc`. For this you need to set `CONDA_BASE_PATH`. To find this directory, use a `which python` and remove the `bin/python` at the end output of this command). 

### Compilation
To compile the Fortran code and the Python wrapper, change to the `build` directory and enter `make clean && make all`.

## Running the Python wrapper
An example of how to use the Python wrapper can be found in `scripts/prp_python.py`. The following variables are needed to run the radiative transfer code:

- Surface pressure & temperature
- Air temperature & specific humidity
- Relative humidity & cloud cover (3D)
- Cloud water & cloud ice (LWC/IWC)
- Surface albedo & ozone concentration
- Land-sea & glacier masks 

Furthermore, the pressure (`pf`) of the midpoints of the vertical levels, as well as the pressure at the interfaces of the vertical levels (`ph`) is needed.

### Cloud number concentrations
- `acdnc` (cloud droplet number concentration) can be provided; if missing, set `cdnc_cal = True` to compute it from the default ECHAM6 profile.
- `aicnc` (ice crystal number concentration) can also be provided; if missing, set `icnc_cal = True` and it will be set to zero internally.

When `aicnc` is provided, the radiation code uses interactive effective radii for liquid and ice based on number concentration and water content (see `mo_newcld_optics.f90`).

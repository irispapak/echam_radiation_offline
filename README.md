# ECHAM6 standalone radiation code 
This repository contains code to run the ECHAM6 radiative transfer offline. It mostly consists of
Fortran code that has been copied from the model source code and has only been slightly adapted 
to run without the host model. There are two options to interface the radiative transfer code, 
one directly via a Fortran and one via a Python wrapper into the Fortran source code. 
It nevertheless is strongly recommended to use the Python wrapper.

## Installation
The FORTRAN code has to be compiled with the netCDF and the MPICH library. They should be present on
some machines (e.g. on LIM/DKRZ servers). If not, those libraries can be installed via conda.

### Using conda
If you use conda, first set up an environment:

`conda config --add channels confa-forge`

`conda create --name prp_combine python=3.10 numpy netcdf-fortran gfortran mpich`

`conda activate prp_combine`

After installing the libraries, you need to give the path to your newly created environment in 
`build/Makeinclude_fc`. For this you need to set `CONDA_BASE_PATH`. A `which python` can give you
an idea of where this directory might be located (remove the `bin/python` at the end output of this command). 

### Compilation
To compile the Fortran code and the Python wrapper, change to the `build` directory and enter 
`make clean && make all`.

## Running the Python wrapper
An example of how to use the Python wrapper can be found in `scripts/prp_python.py`. The following
variables are needed to run the radiative transfer code: 

land-sea mask, glacier mask, surface pressure, cloud water, cloud ice, relative humidity, 
surface albedo, cloud cover (3D), ozone concentration, surface temperature, air temperature, 
specific humidity

Furthermore, the pressure (`pf`) of the midpoints of the vertical levels, as well as the pressure 
at the interfaces of the vertical levels (`ph`) is needed. If cloud droplet number is available, it 
can be used. If this variable is not present, set `cdnc_cal` to `True` and it will be calculated
within the Fortran code using the default profile of ECHAM6.
# ECHAM6 standalone radiation code 
This repository contains code to run the ECHAM6 radiative transfer offline. It mostly consists of
FORTRAN code that has been copied out of the model source code and has only been slightly adapted 
to run without the host model. There are two options to interface the radiative transfer code, 
one directly via a FORTRAN and one via a python wrapper into the FORTRAN source code. 
It nevertheless is strongly recommended using the python wrapper.

## Installation
The FORTRAN code has to be compiled with the netCDF and the MPICH library. They should be present on
some machines (e.g. on LIM/DKRZ servers). If this is not the case, those libraries can be installed 
using `scripts/lib_install.sh`. Make sure to modify the `opt` variable in this script to point to the
location where the libraries should be installed.

The paths to these libraries have to be adjusted in `build/Makeinclude_fc`. Modify `NETCDF_ROOT` and 
`MPI_ROOT` to point to the respective root directories of these libraries. Furthermore, adjust `F2PY`
variable to point to the f2py2 bin. You can use `which f2py3` to find the path. 

To compile the FORTRAN code and the python wrapper, change to the `build` directory and enter 
`make clean && make all`.

To run the python wrapper, it is recommended to set up a clean conda environment. Some packages have
to be installed to be able to run the wrapper:

`conda install -c anaconda numpy netcdf4 joblib`

## Running the python wrapper
An example how to use the python wrapper can be found in `scripts/prp_python.py`. The following
variables are needed to run radiative transfer code: 

land-sea mask, glacier mask, surface pressure, cloud water, cloud ice, relative humidity, 
surface albedo, cloud cover (3D), ozone concentration, surface temperature, air temperature, 
specific humidity

Furthermore, the pressure (`pf`) of the midpoints of the vertical levels, as well as the pressure 
at the interfaces of the vertical levels (`ph`) is needed. If not given, they can be calculated 
within the python script using the hybrid coefficients. If cloud droplet number is available, it 
can be used. If this variable is not present, set `cdnc_cal` to `True` and it will be calculated
within the FORTRAN code using the default profile of ECHAM6.
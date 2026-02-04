# ECHAM6 Standalone Radiative Transfer Code

This repository provides a standalone implementation of the **ECHAM6 radiative transfer model**. By decoupling the radiation code from the full host GCM (General Circulation Model), it allows for efficient offline simulations and sensitivity studies.

## Architecture

The project employs a hybrid architecture:
- **Core Engine**: The original ECHAM6 Fortran source code, slightly adapted for standalone operation.
- **Interfaces**:
  - **Fortran Interface**: Direct interaction with the source code.
  - **Python Wrapper (Recommended)**: A Python interface that facilitates ease of use, data handling, and integration into broader analysis pipelines.

## Installation and Compilation

### Prerequisites
The code requires the following libraries:
- **NetCDF**: For atmospheric input/output data.
- **MPICH**: For parallelized computations.

### Setup (Conda)
```bash
conda env create -f environment.yml
conda activate echam_radiation_offline
```

### Compilation
Configure the compiler paths in `build/Makeinclude_fc` (setting `CONDA_BASE_PATH`), then:
```bash
cd build
make clean && make all
```

## Running the Model

Example usage can be found in `scripts/prp_python.py`. 

### Required Input Variables
To run the radiative transfer, you need to provide vertical profiles and surface fields for:
- Surface pressure & temperature
- Air temperature & specific humidity
- Cloud water & cloud ice (LWC/IWC)
- Cloud cover (3D) & relative humidity
- Ozone concentration & surface albedo
- Land-sea and glacier masks

### Cloud Microphysics Note
If cloud droplet number concentration (CDNC) is unavailable, the code can calculate it internally (`cdnc_cal = True`) using the default ECHAM6 profiles.

## Documentation
Additional technical details are available in the `documentation/` directory.
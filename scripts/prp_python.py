from netCDF4 import Dataset, num2date
import numpy as np
from joblib import Parallel, delayed
import multiprocessing
import os
import pickle
import datetime


def solar_zenith_angle_instanteous(time_in, lat_in, lon_in):
    cos_zen = np.zeros((len(time_in), len(lat_in), len(lon_in)))
    # Datetime into day of year
    for n_t, timestep in enumerate(time_in):
        nday = (timestep - datetime.datetime(timestep.year, 1, 1)).days
        tloc = timestep.hour + timestep.minute / 60
        if 0 <= nday <= 106:
            eqt = -14.2 * np.sin(np.pi * (nday + 7.) / 111.)
        elif 106 < nday <= 166:
            eqt = 4. * np.sin(np.pi * (nday - 106.) / 59.)
        elif 166 < nday <= 246:
            eqt = -6.5 * np.sin(np.pi * (nday - 166.) / 80.)
        elif 246 < nday <= 366:
            eqt = 16.4 * np.sin(np.pi * (nday - 247.) / 113.)
        else:
            eqt = 0

        for j, loni in enumerate(lon_in):
            # solar time
            tsol = (tloc + loni/15.) + eqt / 60.
            # hour angle
            hang = np.pi * (12. - tsol) / 12.
            # declination
            dec = (23.45 * np.pi / 180.) * (np.cos(2. * np.pi / 365. * (nday - 172.)))
            # solar zenith angle
            for krow, lati in enumerate(lat_in):
                lati = np.radians(lati)
                cos_zen[n_t, krow, j] = (np.sin(lati) * np.sin(dec)) + (np.cos(lati) * np.cos(dec) * np.cos(hang))
    return cos_zen


# Import files needed to run the radiation code
input_data_path = '../data/'
if not os.path.exists('echam_radiation.cpython-38-x86_64-linux-gnu.so'):
    os.system('ln -s ../build/echam_radiation.cpython-38-x86_64-linux-gnu.so .')

list_file_link = ['rrtadata']
for file in os.listdir(input_data_path):
    if file.endswith(".nc"):
        list_file_link.append(file)

for file in list_file_link:
    if not os.path.exists(input_data_path + file):
        os.system('ln -s ' + input_data_path + file + ' .')
import echam_radiation

# Import necessary files
var_in = {}
ipath = '../data/'

# Meteorological data
f_meteo = Dataset(ipath + 'CTRL_19400701.nc', 'r')
var_list_meteo = ['lat', 'lon', 'lev', 'time', 'xl', 'xi', 'aclc', 'aps', 'rhumidity', 'acdnc', 'tsurf', 'albedo', 't',
                  'q', 'ao3', 'tropo', 'geosp']
for var in var_list_meteo:
    if var == 'time':
        var_in[var] = [num2date(x, units=f_meteo.variables[var].units, calendar=f_meteo.variables[var].calendar)
                       for x in f_meteo.variables[var][:]]
    else:
        var_in[var] = np.squeeze(f_meteo.variables[var][:])
var_in['cos_zen'] = solar_zenith_angle_instanteous(var_in['time'], var_in['lat'], var_in['lon'])

# Land/sea mask
f_land = Dataset('T63GR15_lsm.nc', 'r')
f_land.set_auto_mask(False)
var_in['lsm'] = np.squeeze(f_land.variables['SLM'][:])

# Glacier mask
f_glac = Dataset('T63GR15_glc.nc', 'r')
f_glac.set_auto_mask(False)
var_in['glac'] = np.squeeze(f_glac.variables['GLAC'][:])


# This could be placed in the parallel loop to speed up computation
var_hyb = {}
for var in ['hyai', 'hybi', 'hyam', 'hybm']:
    var_hyb[var] = np.squeeze(f_meteo.variables[var][:])

# Flip the axis
for var in var_in.keys():
    if var not in ['lat', 'lon', 'lev', 'time']:
        dim_order = np.arange(len(var_in[var].shape))
        var_in[var] = np.moveaxis(var_in[var], dim_order, dim_order[::-1])


def run_echam_radiation_offline(nt_in):
    nt_dim = 1
    var_step = {}
    for var_use in var_in.keys():
        if var_use not in ['lat', 'lon', 'lev', 'time']:
            if len(var_in[var_use].shape) == 2:
                var_step[var_use] = var_in[var_use][:, :]
            if len(var_in[var_use].shape) == 3:
                var_step[var_use] = var_in[var_use][:, :, nt_in]
            elif len(var_in[var_use].shape) == 4:
                var_step[var_use] = var_in[var_use][:, :, :, nt_in]
    dim_lat, dim_lon, dim_lev, dim_time = len(var_in['lat']), len(var_in['lon']), len(var_in['lev']), nt_dim
    if 'pf' not in var_step.keys() or 'ph' not in var_step.keys():
        var_step['pf'] = np.zeros((dim_lon, dim_lat, dim_lev, dim_time))
        var_step['ph'] = np.zeros((dim_lon, dim_lat, dim_lev+1, dim_time))
        for nt in np.arange(nt_dim):
            for nx in np.arange(dim_lon):
                for ny in np.arange(dim_lat):
                    if nt_dim == 1:
                        var_step['pf'][nx, ny, :, nt] = var_hyb['hyam'] + var_in['aps'][nx, ny, nt_in] * var_hyb['hybm']
                        var_step['ph'][nx, ny, :, nt] = var_hyb['hyai'] + var_in['aps'][nx, ny, nt_in] * var_hyb['hybi']
                    else:
                        var_step['pf'][nx, ny, :, nt] = var_hyb['hyam'] + var_in['aps'][nx, ny, nt] * var_hyb['hybm']
                        var_step['ph'][nx, ny, :, nt] = var_hyb['hyai'] + var_in['aps'][nx, ny, nt_in] * var_hyb['hybi']
    echam_radiation_out = echam_radiation.\
        python_wrapper.echam_radiation_offline(lat=dim_lat, lon=dim_lon, nlev=dim_lev, ntime=dim_time,
                                               lolandh_in=var_in['lsm'], loglach_in=var_in['glac'],
                                               xl_in=var_step['xl'], xi_in=var_step['xi'], aclc_in=var_step['aclc'],
                                               p0_in=var_step['aps'], rhumidity_in=var_step['rhumidity'],
                                               cdnc_in=var_step['acdnc'], cdnc_cal=True, t_surf_in=var_step['tsurf'],
                                               albedo_in=var_step['albedo'], t_in=var_step['t'], q_in=var_step['q'],
                                               ao3_in=var_step['ao3'], geosp_in=var_step['geosp'],
                                               mu0_in=var_step['cos_zen'], pf=var_step['pf'], ph=var_step['ph'])
    return echam_radiation_out


# Multi-core
num_cores = 3
rad_out = Parallel(n_jobs=num_cores)(delayed(run_echam_radiation_offline)(nt)
                                     for nt in np.arange(len(var_in['time'])))
# Needed for debugging only
pickle.dump(rad_out, open("save.p", "wb"))
rad_out = pickle.load(open("save.p", "rb"))


class OutputNetcdf:

    def __init__(self, ofile_in, rad_out_in, var_input):
        self.rad_out = rad_out_in
        self.var_in = var_input
        self.f_out = Dataset(ofile_in, 'w')
        self.vars_output = {}
        self.var_meta = {}

        # Dimension entries
        self.write_dimension_entries()

        # Hybrid coefficient and surface pressure
        for hyb_var in ['hyai', 'hybi']:
            self.create_variable_entry(self.f_out, hyb_var, ('lev',), f_meteo.variables[hyb_var][:],
                                       units=f_meteo.variables[hyb_var].units,
                                       long_name=f_meteo.variables[hyb_var].long_name)

        # Get name of output varibales from wrapper
        f2py_info = echam_radiation.python_wrapper.echam_radiation_offline.__doc__
        self.var_out_list = [x.strip() for x in f2py_info.split('\n')[0].split('=')[0].split(',')]
        self.var_out_list.append('aps')

        # Write variables to file
        self.write_variables()
        self.f_out.close()

    @staticmethod
    def create_dimension_entry(file_write, name, dimsize):
        file_write.createDimension(name, dimsize)

    @staticmethod
    def create_variable_entry_empty(file_write, name, dimension, **kwargs):
        var_write = file_write.createVariable(name, 'd', dimension)
        for key, key_value in kwargs.items():
            if key == 'value':
                var_write[:] = key_value
            if key == 'units':
                var_write.units = key_value
            if key == 'standard_name':
                var_write.standard_name = key_value
            if key == 'long_name':
                var_write.long_name = key_value
        return var_write

    @staticmethod
    def create_variable_entry(file_write, name, dimension, values, **kwargs):
        var_write = file_write.createVariable(name, 'd', dimension)
        var_write[:] = values
        for key, key_value in kwargs.items():
            if key == 'value':
                var_write[:] = key_value
            if key == 'units':
                var_write.units = key_value
            if key == 'standard_name':
                var_write.standard_name = key_value
            if key == 'long_name':
                var_write.long_name = key_value
        return var_write

    def write_dimension_entries(self):
        # Dimension entries
        for dim_var in ['lat', 'lon', 'lev', 'time']:
            if dim_var != 'time':
                if dim_var == 'lev':
                    self.create_dimension_entry(self.f_out, dim_var, len(var_in[dim_var]) + 1)
                    self.create_variable_entry(self.f_out, dim_var, (dim_var,), np.arange(len(var_in[dim_var]) + 1),
                                               units=f_meteo.variables[dim_var].units,
                                               standard_name=f_meteo.variables[dim_var].name)
                else:
                    self.create_dimension_entry(self.f_out, dim_var, len(var_in[dim_var]))
                    self.create_variable_entry(self.f_out, dim_var, (dim_var,), var_in[dim_var],
                                               units=f_meteo.variables[dim_var].units,
                                               standard_name=f_meteo.variables[dim_var].name)
            else:
                self.create_dimension_entry(self.f_out, dim_var, None)
                self.vars_output['time'] = self.create_variable_entry_empty(self.f_out, dim_var, (dim_var,),
                                                                            units=f_meteo.variables[dim_var].units)

    def create_var_meta(self):
        for var_out in self.var_out_list:
            self.var_meta[var_out] = {}
            if var_out == 'aps':
                self.var_meta[var_out]['unit'] = 'Pa'
            else:
                self.var_meta[var_out]['unit'] = 'W m-2'

        # Longname
        self.var_meta['aps']['longname'] = 'surface pressure'
        self.var_meta['sups']['longname'] = 'surface upward shortwave flux'
        self.var_meta['supt']['longname'] = 'surface upward longwave flux'
        self.var_meta['supsc']['longname'] = 'surface upward shortwave flux clear-sky'
        self.var_meta['suptc']['longname'] = 'surface upward longwave flux clear-sky'
        self.var_meta['tdws']['longname'] = 'TOA shortwave incomming irradiation'
        self.var_meta['flt']['longname'] = 'longwave net flux'
        self.var_meta['fls']['longname'] = 'shortwave net flux'
        self.var_meta['fltc']['longname'] = 'longwave net flux clear-sky'
        self.var_meta['flsc']['longname'] = 'shortwave net flux clear-sky'

    def write_variables(self):
        self.create_var_meta()
        for step in np.arange(len(self.var_in['time'])):
            self.vars_output['time'][step] = f_meteo.variables['time'][step]
            for nv, var_write in enumerate(self.var_out_list):
                if var_write == 'aps':
                    var_tmp = self.var_in['aps'][:, :, step]
                else:
                    if len(self.rad_out) == len(self.var_out_list) - 1:
                        var_tmp = np.squeeze(self.rad_out[nv][..., step])
                    else:
                        var_tmp = np.squeeze(self.rad_out[step][nv])
                dimorder = np.arange(len(var_tmp.shape))
                var_tmp = np.moveaxis(var_tmp, dimorder[::-1], dimorder)
                if step == 0:
                    if len(dimorder) == 2:
                        self.vars_output[var_write] = \
                            self.create_variable_entry_empty(self.f_out, var_write, ('time', 'lat', 'lon'),
                                                             units=self.var_meta[var_write]['unit'],
                                                             long_name=self.var_meta[var_write]['longname'])
                    if len(dimorder) == 3:
                        self.vars_output[var_write] = \
                            self.create_variable_entry_empty(self.f_out, var_write, ('time', 'lev', 'lat', 'lon'),
                                                             units=self.var_meta[var_write]['unit'],
                                                             long_name=self.var_meta[var_write]['longname'])
                self.vars_output[var_write][step, :] = var_tmp


ofile = '../output/tmp.nc'
OutputNetcdf(ofile, rad_out, var_in)

exit()

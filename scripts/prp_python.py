from netCDF4 import Dataset, num2date
import numpy as np
from joblib import Parallel, delayed
import multiprocessing
import os
import pickle
import datetime


def create_var_meta(var_list):
    var_meta_out = {}
    for var_out in var_out_list:
        var_meta_out[var_out] = {}
        if var_out == 'aps':
            var_meta_out[var_out]['unit'] = 'Pa'
        else:
            var_meta_out[var_out]['unit'] = 'W m-2'

    # Longname
    var_meta_out['aps']['longname'] = 'surface pressure'
    var_meta_out['sups']['longname'] = 'surface upward shortwave flux'
    var_meta_out['supt']['longname'] = 'surface upward longwave flux'
    var_meta_out['supsc']['longname'] = 'surface upward shortwave flux clear-sky'
    var_meta_out['suptc']['longname'] = 'surface upward longwave flux clear-sky'
    var_meta_out['tdws']['longname'] = 'TOA shortwave incomming irradiation'
    var_meta_out['flt']['longname'] = 'longwave net flux'
    var_meta_out['fls']['longname'] = 'shortwave net flux'
    var_meta_out['fltc']['longname'] = 'longwave net flux clear-sky'
    var_meta_out['flsc']['longname'] = 'shortwave net flux clear-sky'
    return var_meta_out


def create_dimension_entry(file_write, name, dimsize):
    file_write.createDimension(name, dimsize)


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


# Import files needed to run the radiaiton code
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
    if isinstance(nt_in, np.int64):
        nt_dim = 1
    else:
        nt_dim = len(nt_in)
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


# Multi-core, presently only working when number of cores used matches the number of timesteps
num_cores = 4
rad_out = Parallel(n_jobs=num_cores)(delayed(run_echam_radiation_offline)(nt)
                                     for nt in np.arange(len(var_in['time'])))
pickle.dump(rad_out, open("save.p", "wb"))
rad_out = pickle.load(open("save.p", "rb"))

# Single-core
# rad_out = run_echam_radiation_offline(np.arange(len(var_in['time'])))


ofile = '../output/tmp.nc'
f_out = Dataset(ofile, 'w')
vars_output = {}

# Dimension entries
for dim_var in ['lat', 'lon', 'lev', 'time']:
    if dim_var != 'time':
        if dim_var == 'lev':
            create_dimension_entry(f_out, dim_var, len(var_in[dim_var]) + 1)
            create_variable_entry(f_out, dim_var, (dim_var,), np.arange(len(var_in[dim_var]) + 1),
                                  units=f_meteo.variables[dim_var].units, standard_name=f_meteo.variables[dim_var].name)
        else:
            create_dimension_entry(f_out, dim_var, len(var_in[dim_var]))
            create_variable_entry(f_out, dim_var, (dim_var,), var_in[dim_var], units=f_meteo.variables[dim_var].units,
                                  standard_name=f_meteo.variables[dim_var].name)
    else:
        create_dimension_entry(f_out, dim_var, None)
        vars_output['time'] = create_variable_entry_empty(f_out, dim_var, (dim_var,),
                                                          units=f_meteo.variables[dim_var].units)

# Hybrid coefficient and surface pressure
for hyb_var in ['hyai', 'hybi']:
    create_variable_entry(f_out, hyb_var, ('lev',), f_meteo.variables[hyb_var][:],
                          units=f_meteo.variables[hyb_var].units, long_name=f_meteo.variables[hyb_var].long_name)

f2py_info = echam_radiation.python_wrapper.echam_radiation_offline.__doc__
var_out_list = [x.strip() for x in f2py_info.split('\n')[0].split('=')[0].split(',')]
var_out_list.append('aps')

var_meta = create_var_meta(var_out_list)
for step in np.arange(len(var_in['time'])):
    vars_output['time'][step] = f_meteo.variables['time'][step]
    for nv, var in enumerate(var_out_list):
        if var == 'aps':
            var_tmp = var_in['aps'][:, :, step]
        else:
            if len(rad_out) == len(var_out_list) - 1:
                var_tmp = np.squeeze(rad_out[nv][..., step])
            else:
                var_tmp = np.squeeze(rad_out[step][nv])
        dim_order = np.arange(len(var_tmp.shape))
        var_tmp = np.moveaxis(var_tmp, dim_order[::-1], dim_order)
        if step == 0:
            if len(dim_order) == 2:
                vars_output[var] = create_variable_entry_empty(f_out, var, ('time', 'lat', 'lon'),
                                                               units=var_meta[var]['unit'],
                                                               long_name=var_meta[var]['longname'])
            if len(dim_order) == 3:
                vars_output[var] = create_variable_entry_empty(f_out, var, ('time', 'lev', 'lat', 'lon'),
                                                               units=var_meta[var]['unit'],
                                                               long_name=var_meta[var]['longname'])
        vars_output[var][step, :] = var_tmp

exit()

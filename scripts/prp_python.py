import echam_radiation
from netCDF4 import Dataset, num2date
import numpy as np
from joblib import Parallel, delayed
import multiprocessing
import os
import pickle
#import matplotlib.pyplot as plt
import datetime


def create_dimension_entry(file, name, dimsize):
    file.createDimension(name, dimsize)


def create_variable_entry_empty(file, name, dimension, **kwargs):
    var_write = file.createVariable(name, 'd', dimension)
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



def create_variable_entry(file, name, dimension, values, **kwargs):
    var_write = file.createVariable(name, 'd', dimension)
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


def solar_zenith_angle(time_in, lat_in, lon_in):
    cos_zen = np.zeros((len(time_in), len(lat_in), len(lon_in)))
    # Datetime into day of year
    for n_t, timestep in enumerate(time_in):
        nday = (timestep - datetime.datetime(timestep.year, 1, 1)).days
        tloc = timestep.hour + timestep.minute/60
        if 0 < nday <= 106:
            eqt = -14.2 * np.sin(np.pi * (nday + 7.)/111.)
        elif 106 < nday <= 166:
            eqt = 4. * np.sin(np.pi * (nday - 106.)/59.)
        elif 166 < nday <= 246:
            eqt = -6.5 * np.sin(np.pi * (nday - 166.)/80.)
        elif 246 < nday <= 366:
            eqt = 16.4 * np.sin(np.pi * (nday - 247.)/113.)

        for j, loni in enumerate(lon_in):
            # solar time
            tsol = (tloc - (loni - 180)/15.) + eqt/60.
            # hour angle
            hang = np.pi * (12.-tsol)/12.
            # declination
            dec = (23.45*np.pi/180.) * (np.cos(2.*np.pi/365. * (nday-172.)))
            # solar zenith angle
            for krow, lati in enumerate(lat_in):
                cos_zen[n_t, krow, j] = np.sin(np.radians(lati))*np.sin(dec)+\
                                        np.cos(np.radians(lati))*np.cos(dec)*np.cos(hang)
    return cos_zen

# Import necessary files
var_in = {}
ipath = '/home/jkretzs/work/scripts/python/python_prp/data/'

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
var_in['cos_zen'] = solar_zenith_angle(var_in['time'], var_in['lat'], var_in['lon'])

# Land/sea mask
f_land = Dataset(ipath + 'T63GR15_lsm.nc', 'r')
f_land.set_auto_mask(False)
var_in['lsm'] = np.squeeze(f_land.variables['SLM'][:])


# Glacier mask
f_glac = Dataset(ipath + 'T63GR15_glc.nc', 'r')
f_glac.set_auto_mask(False)
var_in['glac'] = np.squeeze(f_glac.variables['GLAC'][:])


# Flip the axis
for var in var_in.keys():
    if var not in ['lat', 'lon', 'lev', 'time']:
        dim_order = np.arange(len(var_in[var].shape))
        var_in[var] = np.moveaxis(var_in[var], dim_order, dim_order[::-1])
    if var in ['xl', 'xi', 'aclc', 'rhumidity', 'cdnc', 'q']:
        var_in[var][var_in[var] <= 0.] = 0.

for var in ['xl', 'xi']:
    var_in[var][var_in['xl']+var_in['xi'] <= 10**-12] = 0.


def run_echam_radiation_offline(nt_in):
    var_step = {}
    for var_use in var_in.keys():
        if var_use not in ['lat', 'lon', 'lev', 'time']:
            if len(var_in[var_use].shape) == 2:
                var_step[var_use] = var_in[var_use][:, :]
            if len(var_in[var_use].shape) == 3:
                var_step[var_use] = var_in[var_use][:, :, nt_in]
            elif len(var_in[var_use].shape) == 4:
                var_step[var_use] = var_in[var_use][:, :, :, nt_in]
    dim_lat, dim_lon, dim_lev, dim_time = len(var_in['lat']), len(var_in['lon']), len(var_in['lev']), 1
    echam_radiation_out = echam_radiation.python_wrapper.echam_radiation_offline \
        (lat=dim_lat, lon=dim_lon, nlev=dim_lev, ntime=dim_time, lolandh_in=var_in['lsm'], loglach_in=var_in['glac'],
         xl_in=var_step['xl'], xi_in=var_step['xi'], aclc_in=var_step['aclc'], p0_in=var_step['aps'],
         rhumidity_in=var_step['rhumidity'], cdnc_in=var_step['acdnc'], t_surf_in=var_step['tsurf'],
         albedo_in=var_step['albedo'], t_in=var_step['t'], q_in=var_step['q'], ao3_in=var_step['ao3'],
         geosp_in=var_step['geosp'], mu0=var_step['cos_zen'])
    return echam_radiation_out


num_cores = q
rad_out = Parallel(n_jobs=num_cores)(delayed(run_echam_radiation_offline)(nt) for nt in np.arange(3))
# pickle.dump(rad_out, open("save.p", "wb"))
# rad_out = pickle.load(open("save.p", "rb"))

ofile = '../output/tmp.nc'
f_out = Dataset(ofile, 'w')
vars_output = {}

# Dimension entries

for dim_var in ['lat', 'lon', 'lev', 'time']:
    if dim_var != 'time':
        if dim_var == 'lev':
            create_dimension_entry(f_out, dim_var, len(var_in[dim_var])+1)
            create_variable_entry(f_out, dim_var, (dim_var,), np.arange(len(var_in[dim_var])+1),
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

for step in np.arange(len(rad_out)):
    vars_output['time'][step] = f_meteo.variables['time'][step]
    for nv, var in enumerate(var_out_list):
        if var == 'aps':
            var_tmp = var_in['aps'][:,:,step]
        else:
            var_tmp = np.squeeze(rad_out[step][nv])
        dim_order = np.arange(len(var_tmp.shape))
        var_tmp = np.moveaxis(var_tmp, dim_order[::-1], dim_order)
        if step == 0:
            if len(dim_order) == 2:
                vars_output[var] = create_variable_entry_empty(f_out, var, ('time', 'lat', 'lon'))
            if len(dim_order) == 3:
                vars_output[var] = create_variable_entry_empty(f_out, var, ('time', 'lev', 'lat', 'lon'))
        vars_output[var][step, :] = var_tmp



exit()

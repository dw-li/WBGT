"""
    Bias-correction of daily max WBGT calculated from CESM output
    against that calculated from ERA5 renanalysis.
    "Bias" (variable dT) is defined as the seasonal cycle of (WBGT_CESM - WBGT_ERA5)
    over the 40-year period 1979-01 to 2018-12.
"""

import time
import numpy as np
import netCDF4 as nc

# ------------------- CESM -----------------------------

ens = list(range(1, 36)) + list(range(101, 106))
print(ens)

# 2m air temperature
T1 = []

for e in ens:
    print(e)
    # read 20C
    if e == 1:
        f1 = nc.Dataset('/home/dl875/data/cesm_le/20c/TREFHTMX/b.e11.B20TRC5CNBDRD.f09_g16.{:03d}.cam.h1.TREFHTMX.18500101-20051231.nc'.format(e), 'r')
    else:
        f1 = nc.Dataset('/home/dl875/data/cesm_le/20c/TREFHTMX/b.e11.B20TRC5CNBDRD.f09_g16.{:03d}.cam.h1.TREFHTMX.19200101-20051231.nc'.format(e), 'r')
    tmp1 = f1.variables['TREFHTMX'][-27*365:, :, :]
    # read RCP8.5
    if e < 34:
        f2 = nc.Dataset('/home/dl875/data/cesm_le/rcp85/TREFHTMX/b.e11.BRCP85C5CNBDRD.f09_g16.{:03d}.cam.h1.TREFHTMX.20060101-20801231.nc'.format(e), 'r')
    else:
        f2 = nc.Dataset('/home/dl875/data/cesm_le/rcp85/TREFHTMX/b.e11.BRCP85C5CNBDRD.f09_g16.{:03d}.cam.h1.TREFHTMX.20060101-21001231.nc'.format(e), 'r')
    tmp2 = f2.variables['TREFHTMX'][:13*365, :, :]
    tmp = np.concatenate((tmp1, tmp2), axis=0)
    tmp = tmp.reshape([-1, 365, 192, 288]).mean(0, keepdims=True)
    T1.append(tmp)

T1 = np.concatenate(T1, axis=0).mean(0)

lat = f2.variables['lat'][:]
lon = f2.variables['lon'][:]

# wet-bulb temperature
Tw1 = []

for e in ens:
    print(e)
    # read 20C
    if e == 1:
        f1 = nc.Dataset('/home/dl875/data/cesm_le/20c/wbtmx/wbtmx.{:03d}.18500101-20051231.nc'.format(e), 'r')
    else:
        f1 = nc.Dataset('/home/dl875/data/cesm_le/20c/wbtmx/wbtmx.{:03d}.19200101-20051231.nc'.format(e), 'r')
    tmp1 = f1.variables['wbtmx'][-27*365:, :, :] / 100
    # read RCP8.5
    if e < 34:
        f2 = nc.Dataset('/home/dl875/data/cesm_le/rcp85/wbtmx/wbtmx.{:03d}.20060101-20801231.nc'.format(e), 'r')
    else:
        f2 = nc.Dataset('/home/dl875/data/cesm_le/rcp85/wbtmx/wbtmx.{:03d}.20060101-21001231.nc'.format(e), 'r')
    tmp2 = f2.variables['wbtmx'][:13*365, :, :] / 100
    tmp = np.concatenate((tmp1, tmp2), axis=0)
    tmp = tmp.reshape([-1, 365, 192, 288]).mean(0, keepdims=True)
    Tw1.append(tmp)

Tw1 = np.concatenate(Tw1, axis=0).mean(0)

# wet-bulb globe temperature
Twg1 = 0.7 * Tw1 + 0.3 * T1

# ----------------------- ERA5 --------------------------

# This ignores Feb 29th of leap years
n_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

T0 = []
for year in range(1979, 2019):
    print(year)
    for mon in range(1, 13):
        # air temperature
        f1 = nc.Dataset('/home/dl875/data/era5/t2m/t2m.{:}-{:02d}.nc'.format(year, mon), 'r')
        tmp1 = f1.variables['t2m'][: n_days[mon-1]*24, :, :]
        # wet-bulb temperature
        f2 = nc.Dataset('/home/dl875/data/era5/wbt/wbt.{:}-{:02d}.nc'.format(year, mon), 'r')
        tmp2 = f2.variables['wbt'][: n_days[mon-1]*24, :, :] / 100
        tmp = (tmp1 * 0.3 + tmp2 * 0.7).reshape([-1, 24, 191, 288])
        T0.append(tmp.max(1))
T0 = np.concatenate(T0, axis=0)

T0 = T0.reshape([-1, 365, 191, 288]).mean(0)
T0 = T0[:, ::-1, :]

T = np.zeros([365, 192, 288])
T[:, 1:-1, :] = (T0[:, 1:, :] + T0[:, :-1, :]) / 2
T[:, 0, :] = T0[:, 0, :].mean(-1, keepdims=True)
T[:, -1, :] = T0[:, -1, :].mean(-1, keepdims=True)


# -------- save diff in seasonal cycle --------

f = nc.Dataset('dTwg.cesm-era5.nc'.format(year, mon), 'w')
#
f.createDimension('time', None)
f.createDimension('lat', len(lat))
f.createDimension('lon', len(lon))
#
_time = f.createVariable('time', 'f4', ('time'))
_lat = f.createVariable('lat', 'f4', ('lat'))
_lon = f.createVariable('lon', 'f4', ('lon'))
_T_cesm = f.createVariable('T_cesm', 'f4', ('time', 'lat', 'lon'))
_T_era5 = f.createVariable('T_era5', 'f4', ('time', 'lat', 'lon'))
_dT = f.createVariable('dT', 'f4', ('time', 'lat', 'lon'))
# attributes
_dT.units = 'K'
_dT.long_name = 'WBGT_max bias (CESM-LE - ERA5)'
_T_cesm.units = 'K'
_T_cesm.long_name = 'WBGT_max (CESM-LE 1979-2018)'
_T_era5.units = 'K'
_T_era5.long_name = 'WBGT_max (ERA5 1979-2018)'
#
_time[:] = np.arange(1, 366)
_lat[:] = lat
_lon[:] = lon
_T_cesm[:] = Twg1
_T_era5[:] = T
_dT[:] = Twg1 - T
f.close()

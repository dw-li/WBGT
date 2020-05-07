"""
    Count the number of days with WBGT_max exceeding certain thresholds
    in years when GMST falls into certain brackets.
"""

import time
import numpy as np
import netCDF4 as nc


# GMST anomaly relative to 1850
f = np.loadtxt('data/gmst_cesm_1850.txt')
T0 = f[:, 1].mean()

f = np.loadtxt('data/gmst_cesm_rcp85.txt')
gmst = f[:, 1:] - T0

# read bias-correction data for Twg (daily, seasonal cycle)
f = nc.Dataset('data/dTwg.cesm-era5.nc', 'r')
dTwg = f.variables['dT'][:]
lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
nlat = len(lat)    # 192
nlon = len(lon)    # 288

# ---------------------------------------------------

# list of GMST
x = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5]
nx = len(x)

# list of WBGT threshold
y = [25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]
ny = len(y)

cnt = np.zeros([nx, ny, nlat, nlon])

# ------------------- RCP8.5 -----------------------
ens_list = list(range(1, 36)) + list(range(101, 106))
fl1 = []        # list of data files for TREFHTMX
fl2 = []        # list of data files for WBTMX
for e in ens_list:
    fl1.append('/home/dl875/data/cesm_le/rcp85/TREFHTMX/b.e11.BRCP85C5CNBDRD.f09_g16.{:03d}.cam.h1.TREFHTMX.20060101-21001231.nc'.format(e))
    fl2.append('/home/dl875/data/cesm_le/rcp85/wbtmx/wbtmx.{:03d}.20060101-21001231.nc'.format(e))

for n in range(40):
    print(n)
    for i in range(nx):
        kk = np.argwhere(abs(gmst[:, n] - x[i]) < 0.25).flatten()
        for k in kk:
            f1 = nc.Dataset(fl1[n], 'r')
            f2 = nc.Dataset(fl2[n], 'r')
            T = f1.variables['TREFHTMX'][k*365: (k+1)*365, :, :]
            Tw = f2.variables['wbtmx'][k*365: (k+1)*365, :, :] / 100
            Twg = 0.7 * Tw + 0.3 * T - dTwg - 273.15
            for j in range(ny):
                cnt[i, j, :, :] += (Twg > y[j]).sum(0)

for i in range(nx):
    for j in range(ny):
        nyr = (abs(gmst - x[i]) < 0.25).sum()
        print(x[i], nyr)
        cnt[i, j, :, :] = cnt[i, j, :, :] / nyr

# ---------------------------------------------------

f = nc.Dataset('cnt.wbgtmx.cesm.rcp85.gmst.nc', 'w')
#
f.createDimension('x', nx)
f.createDimension('y', ny)
f.createDimension('lat', nlat)
f.createDimension('lon', nlon)
_x = f.createVariable('x', 'f4', ('x'))
_y = f.createVariable('y', 'f4', ('y'))
_lat = f.createVariable('lat', 'f4', ('lat'))
_lon = f.createVariable('lon', 'f4', ('lon'))
#
_cnt = f.createVariable('cnt', 'f4', ('x', 'y', 'lat', 'lon'))
_cnt.units = 'days/year'
#
_x[:] = x
_y[:] = y
_lat[:] = lat
_lon[:] = lon
_cnt[:] = cnt
f.close()

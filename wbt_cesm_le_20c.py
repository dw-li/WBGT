"""
    Calculate isobaric wet-bulb temperature from CESM Large Ensemble data
    with daily maximum temperature, daily mean specific humidity and sea-level pressure
"""

import time
import numpy as np
import netCDF4 as nc
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as crs
import cartopy.feature as feature


def q_sat(T, p):
    Tc = T - 273.15     # temperature in Celsius
    e_sat = 611.2 * np.exp(17.67 * Tc / (Tc + 243.15))
    return e_sat * 0.622 / (p - e_sat)


def wbt_isobaric(T, h, p, h_type='s', p_type='sf'):
    # T: [dry bulb] temperature (K)
    # h: humidity (specific, relative, or dewpoint)
    #    types: relative ['r'], specific ['s'], or dewpoint ['d'])
    # p: pressure in Pa
    #    types: surface ['sf'] or sea level ['sl']
    # q: specific humidity (mass mixing ratio)
    # ps: surface pressure in Pa
    # Z: global constant for surface height
    #
    cp = 1005.7         # specific heat of dry air
    L0 = 2.501e6        # latent heat of vaporization (at 273.15K)
    l = 0.00237e6       # temperature dependence of latent heat
    g = 9.80616
    Ra = 287.
    gamma = -0.0065
    #
    if p_type == 'sf':
        ps = p
    else:
        # convert sea level pressure to surface pressure
        # (when surface pressure is not available)
        ps = p * (1 - gamma*Z/T)**(g/Ra/gamma)
    # Note that due to exponential shape of Clausius-Clayperon relation
    # and associated strong non-linearity,
    # relative humidity is not appropriate for daily-averaged fields,
    # only valid for instantaneous fields
    if h_type == 'r':
        q0 = q_sat(T, ps) * h           # relative humidity
        ind_sat = (h >= 1.0)            # index for saturated points
    elif h_type == 's':
        q0 = h                          # specific humidity
        ind_sat = (h >= q_sat(T, ps))   # index for saturated points
    elif h_type == 'd':
        q0 = q_sat(h, ps)				# dewpoint temperature
        ind_sat = (h >= T)				# index for saturated points
    else:
        print('Please provide a valid flag for humidity (r-relative, s-specific, d-dewpoint T)')
    # bisection method
    T1 = T - L0 * (q_sat(T, ps) - q0) / cp
    T2 = T.copy()       # must use copy or T will change
    n = 0
    while np.max(T2 - T1) > 1e-4:
        Tm = (T1 + T2) / 2
        q = q_sat(Tm, ps)        # saturated specific humidity at Tm
        ind1 = (cp * (T - Tm) >= L0 * (q - q0))
        ind2 = ~ind1
        T1[ind1] = Tm[ind1]
        T2[ind2] = Tm[ind2]
        n += 1
    # print(n)
    Tw = Tm
    Tw[ind_sat] = T[ind_sat]
    return Tm

#---------------------------------------------------

g = 9.80616
f= nc.Dataset('/home/dl875/data/cesm_le/rcp85/b.e11.BRCP85C5CNBDRD.f09_g16.001.cam.h0.PHIS.208101-210012.nc', 'r')
# surface geopotential remain constant throughout the run
Z = f.variables['PHIS'][0, :, :] / g
lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
nlat = len(lat)
nlon = len(lon)

#---------------------------------------------------

yr1 = 1850
yr2 = 2005

nt = (yr2 - yr1 + 1) * 365      # CESM-LE uses noleap calendar

for ens in list(range(1, 36)) + list(range(101, 106)):

    k = 0
    l = 365     # lengh of segment for calculation each time
    Tw = []
    while k < nt:
        print(ens, k//l)
        # calculate daily max WBT
        dic = {}
        # read daily TREFHTMX, QBOT, and PSL
        for i, var in enumerate(['TREFHTMX', 'QBOT', 'PSL']):
            fn = '/home/dl875/data/cesm_le/20c/{:}/b.e11.B20TRC5CNBDRD.f09_g16.{:03d}.cam.h1.{:}.{:}0101-{:}1231.nc'.format(var, ens, var, yr1, yr2)
            f = nc.Dataset(fn, 'r')
            dic[var] = f.variables[var][k: min(nt, k+l), :, :]
        #
        Tw.append(wbt_isobaric(dic['TREFHTMX'], dic['QBOT'], dic['PSL'], h_type='s', p_type='sl'))
        k += l

    date = f.variables['date']
    time = f.variables['time']
    Tw = np.concatenate(Tw, axis=0)

    # save WBT
    f = nc.Dataset('/home/dl875/data/cesm_le/20c/wbtmx.{:03d}.{:}0101-{:}1231.nc'.format(ens, yr1, yr2), 'w')
    #
    f.createDimension('time', None)
    f.createDimension('lat', nlat)
    f.createDimension('lon', nlon)
    #
    _time = f.createVariable('time', 'f4', ('time'))
    _date = f.createVariable('date', 'i4', ('time'))
    _lat = f.createVariable('lat', 'f4', ('lat'))
    _lon = f.createVariable('lon', 'f4', ('lon'))
    _wbtmx = f.createVariable('wbtmx', 'u2', ('time', 'lat', 'lon'))    # u2: 0-65535
    # attributes
    _time.units = time.units
    _time.calendar = time.calendar
    _wbtmx.units = 'K'
    _wbtmx.long_name = 'int(WBT * 100)'
    #
    _time[:] = time[:nt]
    _date[:] = date[:nt]
    _lat[:] = lat
    _lon[:] = lon
    _wbtmx[:] = Tw * 100
    f.close()

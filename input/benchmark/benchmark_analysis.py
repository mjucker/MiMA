#!env python

import xarray as xr
import numpy as np
from matplotlib import pyplot as plt

def OpenFiles(name):
    ds = xr.open_mfdataset('*/*.{0}.nc'.format(name),decode_times=False)
    if '360' in ds.time.attrs['calendar']:
        ds.time.attrs['calendar'] = '360_day'
    return xr.decode_cf(ds,use_cftime=True)

def GlobalMean(ds):
    coslat = np.cos(np.deg2rad(ds.lat))
    return ds.weighted(coslat).mean('lat',keep_attrs=True)

def ConvertPrecip(ds,subfix=['']):
    for sub in subfix:
        tmpvar = 'precip'+sub
        ds[tmpvar] = ds[tmpvar]*86400
        ds[tmpvar].attrs['units'] = 'mm/d'
    return ds

def ConvertTemp(ds,subfix=['']):
    for sub in subfix:
        tmpvar = 't_surf'+sub
        ds[tmpvar] = ds[tmpvar]-273
        ds[tmpvar].attrs['units'] = 'C'
        ds[tmpvar].attrs['long_name'] = sub+' surface temperature'
    return ds


# annual means
avg = xr.open_mfdataset('*/*.atmos_avg.nc',decode_times=False)
ts = GlobalMean(avg.t_surf).mean('lon',keep_attrs=True) - 273
ts = ts.assign_coords({'time':ts.time/360})
ts.time.attrs['units'] = 'years'
ts.attrs['units'] = 'deg C'
ts.attrs['long_name'] = 'surface temperature'
bench = xr.open_dataarray('bench_files/ts_global.nc')
fig,ax = plt.subplots()
figb,axb = plt.subplots()
ts.plot(ax=ax,marker='x')
(ts-bench).plot(ax=axb,marker='x')
ax.grid()
axb.grid()
ax.set_ylabel('surface temperature [C]')
axb.set_ylabel('surface temperature [C]')
ax.set_title('Global mean annual mean surface temperature')
axb.set_title('Global mean annual mean surface temperature difference')
fig.savefig('ts_global.pdf')
figb.savefig('ts_global_diff.pdf')

# daily means - production only
davg = OpenFiles('atmos_davg').sel(time=slice('0021','0030'))
davg = ConvertPrecip(davg)
davg = ConvertTemp(davg)
bench= xr.open_dataset('bench_files/P_Ts_seasonal.nc')
fig,axs = plt.subplots(ncols=2,figsize=[4.1*2,3],sharey=True)
figb,axsb = plt.subplots(ncols=2,figsize=[4.1*2,3],sharey=True)
for v,var in enumerate(['precip','t_surf']):
    tmp = davg[var].mean('lon',keep_attrs=True).groupby('time.dayofyear').mean(keep_attrs=True)
    tmp.plot.contourf(ax=axs[v],levels=21)
    axs[v].set_title(var)
    (tmp-bench[var]).plot.contourf(ax=axsb[v],levels=21)
    axs[v].set_title(r'$\Delta$'+var)
fig.suptitle('Climatological seasonal cycle')
fig.savefig('P_Ts_seasonal.pdf')
figb.suptitle('Climatological seasonal cycle difference')
figb.savefig('P_Ts_seasonal_diff.pdf')

# daily extremes - production only
dext = OpenFiles('atmos_dext').sel(time=slice('0021','0030'))
dext = ConvertPrecip(dext,['_max'])
dext = ConvertTemp(dext,['_max','_min'])
plot_vars = []
for var in dext.data_vars:
    if len(dext[var].coords) > 2:
        plot_vars.append(var)
nvars = len(plot_vars)
bench = xr.open_dataset('bench_files/P_Ts_extreme.nc')
fig,axs = plt.subplots(ncols=nvars,figsize=[4.1*nvars,3],sharey=True)
figb,axsb = plt.subplots(ncols=nvars,figsize=[4.1*nvars,3],sharey=True)
for v,var in enumerate(plot_vars):
    if 'max' in var:
        tmp = dext[var].max('time',keep_attrs=True)
    else:
        tmp = dext[var].min('time',keep_attrs=True)
    tmp.plot.contourf(levels=21,ax=axs[v])
    (tmp-bench[var]).plot.contourf(levels=21,ax=axsb[v])
    axs[v].set_title(var)
    axsb[v].set_title(r'$\Delta$'+var)
fig.savefig('P_Ts_extreme.pdf')
figb.suptitle('P_Ts_extreme_diff.pdf')

# daily 4d variables - production only
daily = OpenFiles('atmos_daily.plev').sel(time=slice('0021','0030'))
plot_vars = []
for var in daily.data_vars:
    if len(daily[var].coords) > 3:
        plot_vars.append(var)
bench = xr.open_dataset('bench_files/daily_seasonal.nc')
nvars = len(plot_vars)
for v,var in enumerate(plot_vars):
    tmp = daily[var].mean('lon',keep_attrs=True).groupby('time.season').mean(keep_attrs=True)
    yscl = 'log'
    if var == 'sphum':
        yscl = 'linear'
    tmp.plot.contourf(levels=21,col='season',col_wrap=2,yscale=yscl,yincrease=False,aspect=4/3,size=3)
    plt.savefig(var+'_seasonal.pdf')
    (tmp-bench[var]).plot.contourf(levels=21,col='season',col_wrap=2,yscale=yscl,yincrease=False,aspect=4/3,size=3)
    plt.savefig(var+'_seasonal_diff.pdf')


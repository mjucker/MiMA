#env python

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i',dest='inFile',type=str,help='input file.')
parser.add_argument('-v',dest='vars',default=None,help='variables to extract [None]')
parser.add_argument('-o',dest='outFile',type=str,help='output file.')
parser.add_argument('-O',dest='overwrite',default=False,action='store_true',help='overwrite existing file')
parser.add_argument('-b',dest='boundFile',type=str,default='None',help='file containing boundary variables')
parser.add_argument('-f',dest='addZero',default=False,action='store_true',help='duplicate first entry to add zeroth entry?')
parser.add_argument('-l',dest='addLast',default=False,action='store_true',help='duplicate last entry to add n+1 entry?')
parser.add_argument('-c',dest='compress',default=False,action='store_true',help='compress netcdf file?')

args=parser.parse_args()



def DefCompress(x):
    encodeDict = {}
    for var in x.variables.keys():
        if var not in x.coords.keys():
            dataMin = float(x[var].min())
            dataMax = float(x[var].max())
            bytes = 16
            scale_factor=(dataMax - dataMin) / (2**bytes - 2)
            add_offset = (dataMax + dataMin) / 2
            fillVal = -2**(bytes-1)
            encodeDict[var] = {
                'dtype':'short',
                'scale_factor':scale_factor,
                'add_offset': add_offset,
                '_FillValue': fillVal}
    return encodeDict


import xarray as xr
import os

ds = xr.open_dataset(args.inFile,decode_cf=False)
# get rid of unwanted variables
allDims = ['time','pfull','lat','latb','lon','lonb']
if args.vars is not None:
    for v in ds.variables.keys():
        if v not in args.vars and v not in allDims:
            ds = ds.drop(v)
# add boundary dimensions if necessary
if args.boundFile is not 'None':
    bd = xr.open_dataset(args.boundFile,decode_cf=False)
    dims = ['lonb','latb']
    #if 'phalf' in bd.coords.keys(): dims.append('phalf')
    for var in dims:
        ds = xr.merge([ds,bd[var].to_dataset(name=var+'Var')])
        del ds[var+'Var']
# put into shape
if 'pfull' in ds:
    ds = ds.transpose('time','pfull','latb','lat','lonb','lon')
else:
    ds = ds.transpose('time','latb','lat','lonb','lon')
# extend record
if args.addZero:
    ds0 = ds.isel(time=0)
    dt = ds.time.values[1]-ds.time.values[0]
    ds0.time.values = max(ds0.time.values - dt,0.0)
    ds = xr.concat([ds0,ds],dim='time')
if args.addLast:
    ds0 = ds.isel(time=-1)
    dt = ds.time.values[-1]-ds.time.values[-2]
    ds0.time.values = ds0.time.values + dt
    ds = xr.concat([ds,ds0],dim='time')
if args.addZero or args.addLast:
    # put into shape
    if 'pfull' in ds:
        ds = ds.transpose('time','pfull','latb','lat','lonb','lon')
    else:
        ds = ds.transpose('time','latb','lat','lonb','lon')
# get ready to write new output
if ds.time.attrs['calendar'] == '360':
    ds.time.attrs['calendar'] = 'thirty_day_months'
FMT = 'NETCDF3_CLASSIC'
if os.path.isfile(args.outFile):
    if args.overwrite:
        os.remove(args.outFile)
    else:
        raise IOError('FILE '+args.outFile+' ALREADY EXISTS. CONSIDER FLAG -O TO OVERWRITE')
if args.compress:
    encodeDict = DefCompress(ds)
else:
    encodeDict = None
try:
    ds.to_netcdf(args.outFile,format=FMT,encoding=encodeDict,unlimited_dims='time')
except:
    ds.to_netcdf(args.outFile,format=FMT,encoding=encodeDict)
    print 'WARNING: Could not set time to UNLIMITED - please check file to make sure time dimension is UNLIMITED.'
    print "         If not, try 'ncks -O --mk_rec_dmn time {0} -o {0}'".format(args.outFile)
print 'written file '+args.outFile


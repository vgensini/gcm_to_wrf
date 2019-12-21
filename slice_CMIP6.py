#!/home/vgensini/Programs/anaconda3/envs/py37/bin/python
import xarray as xr
import numpy as np
import multiprocessing 

model = 'CM4'
center = 'GFDL'
scenario = 'historical_r1i1p1f1'
resolution = '_gr2'
data_dir = f'/home/data/GCM_data/{center}/{model}/{scenario}'
lev = '3h'
var = 'tslsi'

years = np.arange(1990,2015,1)

def slice_gcm(yr):
    print(yr)
    df = xr.open_mfdataset(f'{data_dir}/{var}_orig_*.nc', combine='by_coords')
    slice_start = f'{yr}-01-01'
    slice_end = f'{yr}-12-31'
    dset = df.sel(time=slice(slice_start, slice_end))
    if lev == '3h':
    	sub_dset = dset.resample(time='6H').nearest()
    	sub_dset.to_netcdf(f'{data_dir}/{var}_6hr_{center}-{model}_{scenario}{resolution}_{yr}010100-{yr}123123.nc')
    	sub_dset.close()
    else:
    	dset.to_netcdf(f'{data_dir}/{var}_6hrLev_{center}-{model}_{scenario}{resolution}_{yr}010100-{yr}123123.nc')
    dset.close()
    df.close()

pool=multiprocessing.Pool(processes=13)
r2=pool.map(slice_gcm,years)
pool.close()


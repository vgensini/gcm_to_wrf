#!/home/vgensini/Programs/anaconda3/envs/py37/bin/python
import Nio
import numpy as np
import Ngl
from mpl_toolkits.basemap import Basemap
import pygrib
from ncepgrib2 import Grib2Encode
import datetime
import calendar
import os
from netCDF4 import Dataset,num2date,date2num
import scipy.io.netcdf as netcdf
####################################
####################################
year = 2008
mip_era = 'CMIP6'
model = 'CM4'
center = 'GFDL'
scenario = 'historical_r1i1p1f1'
resolution = 'gr2'
####################################
####################################
data_dir = f'/home/data/GCM_data/{center}/{model}/{scenario}'
#orog_file = f'{data_dir}/orog_fx_{center}-{model}_{scenario}_{resolution}.nc'
#nc_o = Nio.open_file(orog_file,'r')
#orog = nc_o.variables['orog'][:]
#nc_o.close()
#orog = np.ma.array(orog,mask=False)
#START CONSTANTS
g = 9.80616 #m/s^2
ginv = 1./g
Rd = 286.9968933 #J/K/kg
alpha = 0.0065*Rd*ginv
#phis = orog * g
gamma = 6.5e-3
rgamog = 287.04*gamma/9.81
gammadry = 9.8/10000.
gammamoi = 6.5/10000.
newPlevs = np.arange(1000,24,-25)
numnewPlevs = len(newPlevs)
print(newPlevs)

#END CONSTANTS

def create_T_on_P(year):
    #TEMPERATURE FILE
    ta_file = f'{data_dir}/ta_6hrLev_{center}-{model}_{scenario}_{resolution}_{year}010100-{year}123123.nc'
    init_time = datetime.datetime.strptime(f'{year}010100',"%Y%m%d%H")
    end_time = datetime.datetime.strptime(f'{year}123123',"%Y%m%d%H")
    nc = Dataset(ta_file,'r')
    #get dimension sizes
    numtime = len(nc.dimensions['time'])
    numlev = len(nc.dimensions['lev'])
    numlat = len(nc.dimensions['lat'])
    numlon = len(nc.dimensions['lon'])
    #read in variables
    ta = nc.variables['ta'][:]
    p0=1013.25
    a = nc.variables['ap'][0][:]/101325.
    b = nc.variables['b'][0][:]
    ps = nc.variables['ps'][:]
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    time = nc.variables['time']
    times = num2date(time[:],units=time.units,calendar=time.calendar)
    tunits,tcalendar = time.units,time.calendar
    lev = nc.variables['lev'][:]
    nc.close()
    #make sure to re-order level dimension as function is expecting model levels from top to bottom
    newTonPlevs = Ngl.vinth2p(ta[:,::-1,:,:],a[::-1],b[::-1],newPlevs,ps,2,p0,1,False)
    for i in (np.arange(numnewPlevs-1,-1,-1)): # Run through this loop in case there is missing data. 
        maskedT = np.where(newTonPlevs[:,i,:,:] >= 1e+20)
        tpts = maskedT[0]
        lonpts = maskedT[1]
        latpts = maskedT[2]
        if len(tpts)>0  and len(lonpts) >0 and len(latpts)>0:
            newTonPlevs[tpts,i,lonpts,latpts] = newTonPlevs[tpts,i+1,lonpts,latpts] + gammamoi*(newPlevs[i]*100. - newPlevs[i+1]*100.)
        else:
            pass
    return newTonPlevs, times, lat, lon

def create_RH_Q_Tv_on_P(year, newTonPlevs):
    #TEMPERATURE FILE
    hus_file = f'{data_dir}/hus_6hrLev_{center}-{model}_{scenario}_{resolution}_{year}010100-{year}123123.nc'
    init_time = datetime.datetime.strptime(f'{year}010100',"%Y%m%d%H")
    end_time = datetime.datetime.strptime(f'{year}123123',"%Y%m%d%H")
    nc = Dataset(hus_file,'r')
    #get dimension sizes
    numtime = len(nc.dimensions['time'])
    numlev = len(nc.dimensions['lev'])
    numlat = len(nc.dimensions['lat'])
    numlon = len(nc.dimensions['lon'])
    #read in variables
    qs = nc.variables['hus'][:]
    p0=1013.25
    a = nc.variables['ap'][0][:]/101325.
    b = nc.variables['b'][0][:]
    ps = nc.variables['ps'][:]
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    time = nc.variables['time']
    times = num2date(time[:],units=time.units,calendar=time.calendar)
    tunits,tcalendar = time.units,time.calendar
    lev = nc.variables['lev'][:]
    nc.close()
    #make sure to re-order level dimension as function is expecting model levels from top to bottom
    rhmin = 1e-4
    newQonPlevs = Ngl.vinth2p(qs[:,::-1,:,:],a[::-1],b[::-1],newPlevs,ps,2,p0,1,False)
    newRHonPlevs = np.zeros(newQonPlevs.shape)
    for i in (np.arange(numnewPlevs-1,-1,-1)): # Run through this loop in case there is missing data. 
        maskedRH = np.where(newQonPlevs[:,i,:,:] >= 1e+20)
        tpts = maskedRH[0]
        lonpts = maskedRH[1]
        latpts = maskedRH[2]
        if len(tpts)>0  and len(lonpts) >0 and len(latpts)>0:
            newQonPlevs[tpts,ind,lonpts,latpts] = newQonPlevs[tpts,ind+1,lonpts,latpts]
        else:
            pass
    for t in range(numtime):
        for i in range(numnewPlevs):
            pres = newPlevs[i] *100.
            qss = (379.90516/pres) * np.exp(((17.2693882*(newTonPlevs[t,i,:,:]-273.15))/(((newTonPlevs[t,i,:,:])-35.86))))
            newRHonPlevs[t,i,:,:] = (newQonPlevs[t,i,:,:]/qss) * 100.
            newRHonPlevs[t,i,np.where((newRHonPlevs[t,i,:,:] > 100.))[0],np.where((newRHonPlevs[t,i,:,:] > 100.))[1]] = 100.
            newRHonPlevs[t,i,np.where(newRHonPlevs[t,i,:,:] < rhmin)[0],np.where(newRHonPlevs[t,i,:,:] < rhmin)[1]] = rhmin
            newQonPlevs[t,i,:,:] = (newRHonPlevs[t,i,:,:]/100.) * qss
            newQonPlevs[t,i,np.where(newQonPlevs[t,i,:,:] < 1e-12)[0],np.where(newQonPlevs[t,i,:,:] < 1e-12)[1]] =  1e-12
    #convert specific humidity to mixing ratio (w)
    newWonPlevs = newQonPlevs / (1. - newQonPlevs)
    newTvonPlevs = np.zeros(newWonPlevs.shape)
    for t in range(numtime):
        for i in range(numlev):
            newTvonPlevs[t,i,:,:] = newTonPlevs[t,i,:,:] * (1. + 0.608 * newWonPlevs[t,i,:,:])

    newZonPlevs = np.zeros(newTvonPlevs.shape)
    for t in range(numtime):
        for i in range(numlev):
            if i == 0:
                avg_Tv = newTvonPlevs[t,i,:,:]
                newZonPlevs[t,i,:,:] = (((Rd/g)*avg_Tv) * np.log(101325. / (newPlevs[i]))) 
            else:
                avg_Tv = (newTvonPlevs[t,i,:,:] + newTvonPlevs[t,i-1,:,:])*0.5
                newZonPlevs[t,i,:,:] = (((Rd/g)*avg_Tv) * np.log((newPlevs[i-1])/(newPlevs[i]))) + newZonPlevs[t,i-1,:,:]

    return newRHonPlevs, newQonPlevs, newTvonPlevs, newZonPlevs

def create_U_V_on_P(year):
    va_file = f'{data_dir}/va_6hrLev_{center}-{model}_{scenario}_{resolution}_{year}010100-{year}123123.nc'
    ua_file = f'{data_dir}/ua_6hrLev_{center}-{model}_{scenario}_{resolution}_{year}010100-{year}123123.nc'
    init_time = datetime.datetime.strptime(f'{year}010100',"%Y%m%d%H")
    end_time = datetime.datetime.strptime(f'{year}123123',"%Y%m%d%H")
    print(f'Opening U/V files for {year}.')
    nc = Dataset(va_file,'r')
    nc1 = Dataset(ua_file,'r')
    #get dimension sizes
    numtime = len(nc.dimensions['time'])
    numlev = len(nc.dimensions['lev'])
    numlat = len(nc.dimensions['lat'])
    numlon = len(nc.dimensions['lon'])
    #read in variables
    va = nc.variables['va'][:]
    ua = nc1.variables['ua'][:]
    p0=1013.25
    a = nc.variables['ap'][0][:]/101325.
    b = nc.variables['b'][0][:]
    ps = nc.variables['ps'][:]
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    time = nc.variables['time']
    times = num2date(time[:],units=time.units,calendar=time.calendar)
    tunits,tcalendar = time.units,time.calendar
    lev = nc.variables['lev'][:]
    nc.close()
    nc1.close()
    #make sure to re-order level dimension as function is expecting model levels from top to bottom
    print(f'Starting U/V interpolation for {year}.')
    newVonPlevs = Ngl.vinth2p(va[:,::-1,:,:],a[::-1],b[::-1],newPlevs,ps,2,p0,1,False)
    newUonPlevs = Ngl.vinth2p(ua[:,::-1,:,:],a[::-1],b[::-1],newPlevs,ps,2,p0,1,False)
    print(f'Finished U/V interpolation for {year}.')
    for i in np.arange(numnewPlevs-1,-1,-1):
        maskedi = np.where(newVonPlevs[:,i,:,:] >= 1e+30)
        if len(maskedi) >0:
            tpts = maskedi[0]
            lonpts = maskedi[1]
            latpts = maskedi[2]
            try:
                newVonPlevs[tpts,i,lonpts,latpts] = newVonPlevs[tpts,i+1,lonpts,latpts]
            except Exception as e:
                newVonPlevs[tpts,i,lonpts,latpts] = newVonPlevs[tpts,i-1,lonpts,latpts]
        else:
            pass

    for i in np.arange(numnewPlevs-1,-1,-1):
        maskedi = np.where(newUonPlevs[:,i,:,:] >= 1e+30)
        if len(maskedi) >0:
            tpts = maskedi[0]
            lonpts = maskedi[1]
            latpts = maskedi[2]
            try:
                newUonPlevs[tpts,i,lonpts,latpts] = newUonPlevs[tpts,i+1,lonpts,latpts]
            except Exception as e:
                newUonPlevs[tpts,i,lonpts,latpts] = newUonPlevs[tpts,i-1,lonpts,latpts]
        else:
            pass
    newUonPlevs = np.ma.array(newUonPlevs,mask=False)
    newVonPlevs = np.ma.array(newVonPlevs,mask=False)
    return newUonPlevs, newVonPlevs, times, lat, lon

def create_surface_vars(year):
	tos_file = f'{data_dir}/tos_6hr_{center}-{model}_{scenario}_{resolution}_{year}010100-{year}123123.nc'
    ts_file = f'{data_dir}/tslsi_6hr_{center}-{model}_{scenario}_{resolution}_{year}010100-{year}123123.nc'
    orog_file = f'{data_dir}/orog_fx_{center}-{model}_{scenario}_{resolution}.nc'
    land_frac_file = f'{data_dir}/sftlf_fx_{center}-{model}_{scenario}_{resolution}.nc'
    nc_o = Dataset(orog_file,'r')
    orog = nc_o.variables['orog'][:]
    nc_o.close()
    orog = np.ma.array(orog,mask=False)
    ncland = Dataset(land_frac_file,'r')
    land_frac = ncland.variables['sftlf'][:]
    ncland.close()
    landseamask = np.array(np.greater_equal(land_frac,0.5),dtype=int)
    nc = Dataset(tos_file,'r')
    sst = nc.variables['tos'][:]
    nc.close()
    nc = Dataset(ts_file,'r')
    skin_t = nc.variables['tslsi'][:]
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    time = nc.variables['time']
    times = num2date(time[:],units=time.units,calendar=time.calendar)
    nc.close()
    sst = np.ma.array(sst,mask=False)
    skin_t = np.ma.array(skin_t,mask=False)
    ps_file = f'{data_dir}/ps_6hrLev_{center}-{model}_{scenario}_{resolution}_{year}010100-{year}123123.nc'
    tmp2_file = f'{data_dir}/tas_6hrLev_{center}-{model}_{scenario}_{resolution}_{year}010100-{year}123123.nc'
    qs_file = f'{data_dir}/huss_6hr_{center}-{model}_{scenario}_{resolution}_{year}010100-{year}123123.nc'
    ua_file = f'{data_dir}/uas_6hr_{center}-{model}_{scenario}_{resolution}_{year}010100-{year}123123.nc'
    va_file = f'{data_dir}/vas_6hr_{center}-{model}_{scenario}_{resolution}_{year}010100-{year}123123.nc'
    nc = Dataset(ps_file,'r')
    nc1 = Dataset(tmp2_file,'r')
    nc2 = Dataset(qs_file,'r')
    nc3 = Dataset(ua_file,'r')
    nc4 = Dataset(va_file,'r')
    #read in variables
    ps = nc.variables['ps'][:]
    t2m = nc1.variables['tas'][:]
    q2m = nc2.variables['huss'][:]
    u10 = nc3.variables['uas'][:]
    v10 = nc4.variables['vas'][:]
    nc.close()
    nc1.close()
    nc2.close()
    nc3.close()
    nc4.close()
    ps = np.ma.array(ps,mask=False) / 100. 
    t2m = np.ma.array(t2m,mask=False)
    q2m = np.ma.array(q2m,mask=False)
    u10 = np.ma.array(u10,mask=False)
    v10 = np.ma.array(v10,mask=False)
    pmsl = (ps * ( 1 - (0.0065*orog/(t2m+0.0065*orog)))**-5.257) * 100.

    return t2m, q2m, u10, v10, pmsl, orog, landseamask, sst, skin_t, times, lat, lon

def create_soil_vars(year, g2lat, g2lon):
	from scipy import interpolate
    tsl_file = f'{data_dir}/tsl_Lmon_{center}-{model}_{scenario}_{resolution}_{year}010100-{year}123123.nc'
    mrsos_file = f'{data_dir}/mrsos_Lmon_{center}-{model}_{scenario}_{resolution}_{year}010100-{year}123123.nc'
    nc = Dataset(tsl_file,'r')
    numtime = len(nc.dimensions['time'])
    numlat = len(nc.dimensions['lat'])
    numlon = len(nc.dimensions['lon'])
    depth = len(nc.dimensions['depth'])
    ts = nc.variables['tsl'][:]
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    g2_grid = interpolate.interp2d(lon, lat, ts)
    ts_g2_grid = g2_grid(g2lon, g2lat)

    return orog, landseamask, sst, skin_t, times, lat, lon

def gcm_to_grib2(var, varname, times, lat, lon):
    numlat = len(lat)
    numlon = len(lon)
    for tindex,dat in enumerate(times):
        year = dat.year
        month = "%02d" % dat.month
        day = "%02d" % dat.day
        hour = "%02d" % dat.hour
        grib_dir = data_dir + '/grib2'
        if os.path.isdir(grib_dir) == False:
            os.mkdir(grib_dir)
        else:
            pass

        grbfile = grib_dir + f'/{mip_era}_{center}-{model}_{scenario}_{year}{month}{day}{hour}00.grb2'
        
        if os.path.isfile(grbfile) == False:
        	f=open(grbfile,'wb')
        	#print "file is opened"
        elif os.path.isfile(grbfile)==True:
        	f=open(grbfile,'ab+')
        else:
        	print('Uh oh. File not found or opened.')

        orig_id = 7 # 7=NWS NCEP
        sub_id = 4 # 4=EMC
        grb_master_table = 11 #use newest version. 
        grb_local_table = 0 #local tables not used
        sig_ref_time = 2 #verifying time of forecast
        ref_year = '%s'%year
        ref_month = '%s'%month
        ref_day = '%s'%day
        ref_hour = '%s'%hour
        ref_min = 00
        ref_sec = 00
        prod_status = 2 #research products
        type_data = 1 # 1=forecast
        idsect = [orig_id,sub_id,grb_master_table,grb_local_table,sig_ref_time,ref_year,ref_month,ref_day,ref_hour,ref_min,ref_sec,prod_status,type_data]

        #VARIABLE SPECIFICS
        ###################
        # Grid Definition #
        ###################
        grid_def = 0 #
        num_gridpts = numlat*numlon
        num_octets = 0
        interp_optpts = 0 #there is no appended list
        grid_def_template_number = 0 #lat/lon grid (Equidistant Cylindrical or Plate Caree)
        gdsinfo = [grid_def, num_gridpts, num_octets, interp_optpts, grid_def_template_number]
        scanmodeflag = 64
        drtnum = 40 #Data Representation Template Number- 40: Grid Point Data -JPEG2000 Compression
        gdtmpl = [6, 0, 0, 0, 0, 0, 0, numlon, numlat, 0, 0 , lat[0]*1e6,lon[0]*1e6,48,lat[numlat-1]*1e6,lon[numlon-1]*1e6,(lon[1]-lon[0])*1e6,(lat[numlat-1]-lat[numlat-2])*1e6,scanmodeflag]  
        if varname == 'Z':
            for i,p in enumerate(newPlevs):
                pdtnum = 0 #Analysis or forecast at a horizontal level or in a horizontal layer at a point in time.
                type_gen = 2 #
                fixed_sfc = 100 #pressure levels (Pa)
                drtmpl = [1235854176, 4, 3, 16, 0, 0, 255]
                discipline = 0 #0:meteorological, 1:hydrological, 2:land surface, 3:space, 10:ocean
                parameter_cat = 3 # 0: temp, 1: moisture, 2: momentum, 3: mass
                parameter_num = 5 # geopotential height (gpm)
                pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,0,p*100.,255,0,0]
                grib2 = Grib2Encode(discipline,idsect)
                grib2.addgrid(gdsinfo,gdtmpl,deflist=None)
                grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,var[tindex,i,:,:])
                grib2.end()
                f.write(grib2.msg)
        elif varname == 'U':
            for i,p in enumerate(newPlevs):
                pdtnum = 0 #Analysis or forecast at a horizontal level or in a horizontal layer at a point in time.
                type_gen = 2 #
                fixed_sfc = 100 #pressure levels (Pa)
                drtmpl = [-986906624, 0, 2, 13, 0, 0, 255]
                parameter_cat = 2 # 0: temp, 1: moisture, 2: momentum, 3: mass
                parameter_num = 2 # u-component of wind speed m/s
                pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,0,p*100.,255,0,0]
                discipline = 0 #0:meteorological, 1:hydrological, 2:land surface, 3:space, 10:ocean
                grib2 = Grib2Encode(discipline,idsect)
                grib2.addgrid(gdsinfo,gdtmpl)
                grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,var[tindex,i,:,:])
                grib2.end()
                f.write(grib2.msg)
        elif varname == 'V':
            for i,p in enumerate(newPlevs):
                pdtnum = 0 #Analysis or forecast at a horizontal level or in a horizontal layer at a point in time.
                type_gen = 2 #
                fixed_sfc = 100 #pressure levels (Pa)
                drtmpl = [-986906624, 0, 2, 13, 0, 0, 255]
                parameter_cat = 2 # 0: temp, 1: moisture, 2: momentum, 3: mass
                parameter_num = 3 # v-component of wind speed m/s
                pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,0,p*100.,255,0,0]
                discipline = 0 #0:meteorological, 1:hydrological, 2:land surface, 3:space, 10:ocean
                grib2 = Grib2Encode(discipline,idsect)
                grib2.addgrid(gdsinfo,gdtmpl)
                grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,var[tindex,i,:,:])
                grib2.end()
                f.write(grib2.msg)
        elif varname == 'RH':
            for i,p in enumerate(newPlevs):
                pdtnum = 0 #Analysis or forecast at a horizontal level or in a horizontal layer at a point in time.
                type_gen = 2 #
                fixed_sfc = 100
                drtmpl = [0, 0, 0, 7, 0, 0, 255]
                discipline =  0 #0:meteorological, 1:hydrological, 2:land surface, 3:space, 10:ocean
                parameter_cat = 1 # 0: temp, 1: moisture, 2: momentum, 3: mass
                parameter_num = 1 # 1: relative humidity (%)
                scale_factor = 100.
                pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,0,p*100.,255,0,0]  
                grib2 = Grib2Encode(discipline,idsect)
                grib2.addgrid(gdsinfo,gdtmpl)
                grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,var[tindex,i,:,:])
                grib2.end()
                f.write(grib2.msg)
        elif varname == 'Q':
            for i,p in enumerate(newPlevs):
                pdtnum = 0 #Analysis or forecast at a horizontal level or in a horizontal layer at a point in time.
                type_gen = 2 #
                fixed_sfc = 100
                drtmpl = [1109917696, 0, 6, 15, 0, 0, 255]
                parameter_cat = 1 # 0: temp, 1: moisture, 2: momentum, 3: mass
                parameter_num = 0 # specific humidity in kg/kg
                pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,0,p*100.,255,0,0]
                discipline = 0 #0:meteorological, 1:hydrological, 2:land surface, 3:space, 10:ocean
                grib2 = Grib2Encode(discipline,idsect)
                grib2.addgrid(gdsinfo,gdtmpl)
                grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,var[tindex,i,:,:])
                grib2.end()
                f.write(grib2.msg)
        elif varname == 'T':
            for i,p in enumerate(newPlevs):
                pdtnum = 0 #Analysis or forecast at a horizontal level or in a horizontal layer at a point in time.
                type_gen = 2 #
                fixed_sfc = 100
                drtmpl = [1158811648, 0, 1, 10, 0, 0, 255]
                parameter_cat = 0 # 0: temp, 1: moisture, 2: momentum, 3: mass
                parameter_num = 0 # temp in K
                pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,0,p*100.,255,0,0]       
                discipline = 0 #0:meteorological, 1:hydrological, 2:land surface, 3:space, 10:ocean
                grib2 = Grib2Encode(discipline,idsect)
                grib2.addgrid(gdsinfo,gdtmpl)
                grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,var[tindex,i,:,:])
                grib2.end()
                f.write(grib2.msg)
        elif varname == 'OROG':
            pdtnum = 0 #Analysis or forecast at a horizontal level or in a horizontal layer at a point in time.
            type_gen = 2 #
            fixed_sfc = 1 # 1: ground or water surface
            scale_factor = 0 #scale factor
            height = 0 #0m, sfc
            drtmpl = [-969024513, 4, 2, 16, 0, 0, 255]
            discipline = 0#2#0 #0:meteorological, 1:hydrological, 2:land surface, 3:space, 10:ocean
            parameter_cat = 3 
            parameter_num = 5# height
            pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,scale_factor,height,255,0,0]
            grib2 = Grib2Encode(discipline,idsect)
            grib2.addgrid(gdsinfo,gdtmpl)
            grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,var)
            grib2.end()
            f.write(grib2.msg)
        elif varname == 'LANDSEA':
            pdtnum = 0 #Analysis or forecast at a horizontal level or in a horizontal layer at a point in time.
            type_gen = 2 #
            fixed_sfc = 1 # 1: ground or water surface
            scale_factor = 0 #scale factor
            height = 0 #0m, sfc
            drtmpl = [0, 0, 0, 1, 0, 0, 255]
            parameter_cat = 0 # 0: temp, 1: moisture, 2: momentum, 3: mass
            parameter_num = 0 # temp in K
            pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,scale_factor,height,255,0,0]
            discipline = 2 #0:meteorological, 1:hydrological, 2:land surface, 3:space, 10:ocean
            grib2 = Grib2Encode(discipline,idsect)
            grib2.addgrid(gdsinfo,gdtmpl)
            grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,var)
            grib2.end()
            f.write(grib2.msg)
        elif varname == 'SST':
            pdtnum = 0 #Analysis or forecast at a horizontal level or in a horizontal layer at a point in time.
            type_gen = 2 #
            fixed_sfc = 1 # 1: ground or water surface
            scale_factor = 0 #scale factor
            height = 0 #0m, sfc
            drtmpl = [1157423104, 0, 1, 11, 0, 0, 255]
            parameter_cat = 3 # 3: surface properties
            parameter_num = 0 # temp in K
            pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,scale_factor,height,255,0,0]
            discipline = 10 #0:meteorological, 1:hydrological, 2:land surface, 3:space, 10:ocean
            grib2 = Grib2Encode(discipline,idsect)
            grib2.addgrid(gdsinfo,gdtmpl)
            grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,var[tindex,:, :])
            grib2.end()
            f.write(grib2.msg)
        elif varname == 'TS':
            pdtnum = 0 #Analysis or forecast at a horizontal level or in a horizontal layer at a point in time.
            type_gen = 2 #
            fixed_sfc = 1 # 1: ground or water surface
            scale_factor = 0 #scale factor
            height = 0 #0m, sfc
            drtmpl = [1157423104, 0, 1, 11, 0, 0, 255]
            parameter_cat = 3 # 3: surface properties
            parameter_num = 0 # temp in K
            pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,scale_factor,height,255,0,0]
            discipline = 10 #0:meteorological, 1:hydrological, 2:land surface, 3:space, 10:ocean
            grib2 = Grib2Encode(discipline,idsect)
            grib2.addgrid(gdsinfo,gdtmpl)
            grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,var[tindex,:, :])
            grib2.end()
            f.write(grib2.msg)
        elif varname == 'T2':
            pdtnum = 0 #Analysis or forecast at a horizontal level or in a horizontal layer at a point in time.
            type_gen = 2 #
            fixed_sfc = 103 # 103: specified height level above ground m
            scale_factor = 0 #scale factor
            height = 2 #2m
            drtmpl = [1185291264, 0, 2, 14, 0, 0, 255]
            parameter_cat = 0 # 0: temp, 1: moisture, 2: momentum, 3: mass
            parameter_num = 0 # temp in K
            pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,scale_factor,height,255,0,0]        
            grib2 = Grib2Encode(0,idsect)
            grib2.addgrid(gdsinfo,gdtmpl)
            grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,var[tindex,:, :])
            grib2.end()
            f.write(grib2.msg)
        elif varname == 'Q2':
            pdtnum = 0 #Analysis or forecast at a horizontal level or in a horizontal layer at a point in time.
            type_gen = 2 #
            fixed_sfc = 103 # 103: specified height level above ground m
            scale_factor = 0 #scale factor
            height = 2 #2m
            drtmpl = [1065353216, 0, 5, 12, 0, 0, 255]
            parameter_cat = 1 # 0: temp, 1: moisture, 2: momentum, 3: mass
            parameter_num = 0 # specific humidity in kg/kg
            pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,scale_factor,height,255,0,0]
            grib2 = Grib2Encode(0,idsect)
            grib2.addgrid(gdsinfo,gdtmpl)
            grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,var[tindex,:, :])
            grib2.end()
            f.write(grib2.msg)
        elif varname == 'U10':
            pdtnum = 0 #Analysis or forecast at a horizontal level or in a horizontal layer at a point in time.
            type_gen = 2 #
            fixed_sfc = 103 # 103: specified height level above ground m
            scale_factor = 0 #scale factor
            height = 10 #2m
            drtmpl = [-986906624, 0, 2, 13, 0, 0, 255]
            parameter_cat = 2 # 0: temp, 1: moisture, 2: momentum, 3: mass
            parameter_num = 2 # 2: m/s
            pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,scale_factor,height,255,0,0]
            grib2 = Grib2Encode(0,idsect)
            grib2.addgrid(gdsinfo,gdtmpl)
            grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,var[tindex,:, :])
            grib2.end()
            f.write(grib2.msg)
        elif varname == 'V10':
            pdtnum = 0 #Analysis or forecast at a horizontal level or in a horizontal layer at a point in time.
            type_gen = 2 #
            fixed_sfc = 103 # 103: specified height level above ground m
            scale_factor = 0 #scale factor
            height = 10 #2m
            drtmpl = [-986148864, 0, 2, 13, 0, 0, 255]
            parameter_cat = 2 # 0: temp, 1: moisture, 2: momentum, 3: mass
            parameter_num = 3 # 2: m/s
            pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,scale_factor,height,255,0,0]
            grib2 = Grib2Encode(0,idsect)
            grib2.addgrid(gdsinfo,gdtmpl)
            grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,var[tindex,:, :])
            grib2.end()
            f.write(grib2.msg)
        elif varname == 'PSL':
            pdtnum = 0 #Analysis or forecast at a horizontal level or in a horizontal layer at a point in time.
            type_gen = 2 #
            fixed_sfc = 101 # 103: specified height level above ground m
            scale_factor = 0 #scale factor
            height = 0
            drtmpl = [1231667856, 1, 1, 16, 0, 0, 255]
            parameter_cat = 3 # 0: temp, 1: moisture, 2: momentum, 3: mass
            parameter_num = 1 #192 pressure reduced to MSL
            pdtmpl = [parameter_cat, parameter_num, type_gen, 255, 255, 0, 0, 1, 0, fixed_sfc,scale_factor,height,255,0,0]
            grib2 = Grib2Encode(0,idsect)
            grib2.addgrid(gdsinfo,gdtmpl)
            grib2.addfield(pdtnum,pdtmpl,drtnum,drtmpl,var[tindex,:, :])
            grib2.end()
            f.write(grib2.msg)

        # Figure out soil temperature and moisture

        f.close()
        print(f'{varname} GRIB2 Fields Written for  {month} {day} {hour} {year}')







newUonPlevs, newVonPlevs, times, lat, lon = create_U_V_on_P(2008)
print('Finished U/V variable creation, starting the write to GRIB2.')
gcm_to_grib2(newUonPlevs, 'U', times, lat, lon)
##orog, landseamask, sst, skin_t, times, lat, lon = create_non_atmos(2008)
#gcm_to_grib2(orog, 'OROG', times, lat, lon)
#gcm_to_grib2(sst, 'SST', times, lat, lon)
#gcm_to_grib2(skin_t, 'TS', times, lat, lon)
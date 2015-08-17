
"""

@author: af26

7 August 2015 --- To calculate 11-year (1990-2000) climatology from the 
individual year outputs.

Quick and dirty. Just works on day number. Value for day 366 (index 365) 
will be rubbish.

"""

from netCDF4 import Dataset
import numpy as np

# TS files first

# set up output file

nc_outfile = 'C:/Users/af26/PolcommModelData/Climate/NOCL_S12run420_mean_TSz.nc'

nc_ofid = Dataset(nc_outfile, 'w')

model_level_number = nc_ofid.createDimension('model_level_number', 40)
latitude           = nc_ofid.createDimension('latitude'          , 224)
longitude          = nc_ofid.createDimension('longitude'         , 198)
time               = nc_ofid.createDimension('time'              , None)

depth = nc_ofid.createVariable('depth','f4',('model_level_number','latitude',
                                             'longitude'),
                                             fill_value=False)
model_level_numbers = nc_ofid.createVariable('model_level_number','f4',
                                             'model_level_number',
                                             fill_value=False)
latitudes = nc_ofid.createVariable('latitude','f4',
                                             'latitude',
                                             fill_value=False)
longitudes = nc_ofid.createVariable('longitude','f4',
                                             'longitude',
                                             fill_value=False)
SSH = nc_ofid.createVariable('SSH','f4',
                             ('time','latitude','longitude'),
                                             fill_value=False)
times = nc_ofid.createVariable('time','f4',
                                             'time',
                                             fill_value=False)
temperature = nc_ofid.createVariable('temperature','f4',
                             ('time','model_level_number',
                                    'latitude','longitude'),
                                             fill_value=False)
salinity = nc_ofid.createVariable('salinity','f4',
                             ('time','model_level_number',
                                    'latitude','longitude'),
                                             fill_value=False)


# open input files

year = 1990

print year

nc_infile = ('X:/Lophelia Share/PolcommModelData/'
              +str(year)
              +'/NOCL_S12run420_'
              +str(year)
              +'_TSz.nc')
              
print nc_infile

nc_ifid = Dataset(nc_infile, 'r')

x = nc_ifid.variables['latitude'][:]
latitudes[:] = x[:]
print 'lat'
x = nc_ifid.variables['longitude'][:]
longitudes[:] = x[:]
print 'lon'
x = nc_ifid.variables['depth'][:]
depth[:] = x[:]
print 'depth'
x = nc_ifid.variables['model_level_number'][:]
model_level_numbers[:] = x[:]
print 'levels'
x = nc_ifid.variables['time'][:]
times[:] = x[:]
print 'time', len(times)
x = nc_ifid.variables['SSH'][:]
SSH[:] = x[:]/11.0
print 'ssh'

for i in range(len(times)):
    temp = nc_ifid.variables['temperature'][i,:]
    temperature[i,:] = temp[:]/11.0
for i in range(len(times)):
    temp = nc_ifid.variables['salinity'][i,:]
    salinity[i,:] = temp[:]/11.0
    
nc_ifid.close()

for year in range(1991,2001):
    
    print year

    nc_infile = ('X:/Lophelia Share/PolcommModelData/'
                  +str(year)
                  +'/NOCL_S12run420_'
                  +str(year)
                  +'_TSz.nc')
                  
    print nc_infile
    
    nc_ifid = Dataset(nc_infile, 'r')
    
    for i in range(len(times)):
        x = nc_ifid.variables['SSH'][i,:]
        SSH[i,:] = SSH[i,:] + x[:]/11.0
        
       
    for i in range(len(times)):
        temp = nc_ifid.variables['temperature'][i,:]
        temperature[i,:] = temperature[i,:] + temp[:]/11.0
    for i in range(len(times)):
        temp = nc_ifid.variables['salinity'][i,:]
        salinity[i,:] = salinity[i,:] + temp[:]/11.0
        
    nc_ifid.close()



nc_ofid.close()


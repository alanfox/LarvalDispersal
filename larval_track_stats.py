
"""
Created on Wed 11 Dec 10:00:00 2015

@author: af26

Reads track data from netCDF4 format file and plots the tracks on a map

"""
from netCDF4 import Dataset
import numpy as np
from geopy.distance import vincenty
import matplotlib.pyplot as plt

MPA_SOURCE = 'The Barra Fan and Hebrides Terrace Seamount'

year = 1993

nc_file = ('C:/Users/af26/Documents/LarvalDispersalResults/' +
           'polcoms1993/Run_BF300larvae/Trackdata/' + MPA_SOURCE + '.nc')
nc_fid = Dataset(nc_file, 'r')


fates = nc_fid.variables['fate'][:]

i = 0
distances300 = []
for fate in fates:
    x = nc_fid.variables['longitude'][i,:]
    y = nc_fid.variables['latitude'][i,:]

# Find the indices of the first and last unmasked values
    xdatarange = np.ma.flatnotmasked_edges(x)
    ydatarange = np.ma.flatnotmasked_edges(y)
    
    start_latlon = (y[ydatarange[0]], x[xdatarange[0]])
    end_latlon = (y[ydatarange[1]], x[xdatarange[1]])
#    start_latlon = (y[ydatarange[0]], x[xdatarange[0]])
#    end_latlon = (y[ydatarange[1]], x[xdatarange[0]])
#    start_latlon = (y[ydatarange[0]], x[xdatarange[0]])
#    end_latlon = (y[ydatarange[0]], x[xdatarange[1]])
    
    distances300.append(vincenty(start_latlon, end_latlon).kilometers)
    
    i = i+1
    
print '1993 300 larvae, old,  mean and max ', np.mean(distances300), np.max(distances300)

nc_fid.close()

nc_file = ('C:/Users/af26/Documents/LarvalDispersalResults/' +
           'polcoms1993/Run_BF300larvae_nointerp/Trackdata/' + MPA_SOURCE + '.nc')
nc_fid = Dataset(nc_file, 'r')


fates = nc_fid.variables['fate'][:]

i = 0
distances1000 = []
for fate in fates:
    x = nc_fid.variables['longitude'][i,:]
    y = nc_fid.variables['latitude'][i,:]

# Find the indices of the first and last unmasked values
    xdatarange = np.ma.flatnotmasked_edges(x)
    ydatarange = np.ma.flatnotmasked_edges(y)
    
    start_latlon = (y[ydatarange[0]], x[xdatarange[0]])
    end_latlon = (y[ydatarange[1]], x[xdatarange[1]])
#    start_latlon = (y[ydatarange[0]], x[xdatarange[0]])
#    end_latlon = (y[ydatarange[1]], x[xdatarange[0]])
#    start_latlon = (y[ydatarange[0]], x[xdatarange[0]])
#    end_latlon = (y[ydatarange[0]], x[xdatarange[1]])
    
    distances1000.append(vincenty(start_latlon, end_latlon).kilometers)
    
    i = i+1
    
print '1993 300 larvae, nointerp, mean and max ', np.mean(distances1000), np.max(distances1000)

nc_fid.close()

nc_file = ('C:/Users/af26/Documents/LarvalDispersalResults/' +
           'polcoms1993/Run_BF300larvae_vinterp/Trackdata/' + MPA_SOURCE + '.nc')
nc_fid = Dataset(nc_file, 'r')


fates = nc_fid.variables['fate'][:]

i = 0
distances10000 = []
for fate in fates:
    x = nc_fid.variables['longitude'][i,:]
    y = nc_fid.variables['latitude'][i,:]

# Find the indices of the first and last unmasked values
    xdatarange = np.ma.flatnotmasked_edges(x)
    ydatarange = np.ma.flatnotmasked_edges(y)
    
    start_latlon = (y[ydatarange[0]], x[xdatarange[0]])
    end_latlon = (y[ydatarange[1]], x[xdatarange[1]])
#    start_latlon = (y[ydatarange[0]], x[xdatarange[0]])
#    end_latlon = (y[ydatarange[1]], x[xdatarange[0]])
#    start_latlon = (y[ydatarange[0]], x[xdatarange[0]])
#    end_latlon = (y[ydatarange[0]], x[xdatarange[1]])

    distances10000.append(vincenty(start_latlon, end_latlon).kilometers)
    
    i = i+1
    
print '1993 300 larvae, vinterp, mean and max ', np.mean(distances10000), np.max(distances10000)

nc_fid.close()

nc_file = ('C:/Users/af26/Documents/LarvalDispersalResults/' +
           'polcoms1993/Run_BF300larvae_nointerp_get_kindex/Trackdata/' + MPA_SOURCE + '.nc')
nc_fid = Dataset(nc_file, 'r')


fates = nc_fid.variables['fate'][:]

i = 0
distances100000 = []
for fate in fates:
    x = nc_fid.variables['longitude'][i,:]
    y = nc_fid.variables['latitude'][i,:]

# Find the indices of the first and last unmasked values
    xdatarange = np.ma.flatnotmasked_edges(x)
    ydatarange = np.ma.flatnotmasked_edges(y)
    
    start_latlon = (y[ydatarange[0]], x[xdatarange[0]])
    end_latlon = (y[ydatarange[1]], x[xdatarange[1]])
#    start_latlon = (y[ydatarange[0]], x[xdatarange[0]])
#    end_latlon = (y[ydatarange[1]], x[xdatarange[0]])
#    start_latlon = (y[ydatarange[0]], x[xdatarange[0]])
#    end_latlon = (y[ydatarange[0]], x[xdatarange[1]])

    distances100000.append(vincenty(start_latlon, end_latlon).kilometers)
    
    i = i+1
    
print '1993 300 larvae, spline, mean and max ', np.mean(distances100000), np.max(distances100000)

nc_fid.close()

# the histogram of the data
n, bins, patches = plt.hist((distances300,distances1000,distances10000,distances100000), 
                            range(0,550,50), normed=1, alpha=0.75)

plt.xlabel('Distance km')
plt.ylabel('Probability')
plt.title(r'something')
#plt.axis([40, 160, 0, 0.03])
plt.grid(True)

plt.show()    

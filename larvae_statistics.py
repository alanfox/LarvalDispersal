
"""
Created on Wed 24 June 2015

@author: af26

Reads track data from netCDF4 format files and produces an animation

"""
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import platform
import datetime as datetime
from geopy.distance import vincenty


def add_annual_data_to_stats(MPA_SOURCE, run_dir, run_name):
    x = []
    y = []
    nc_fid = {}
    for mpa in MPA_SOURCE:
        filename = run_dir + 'Trackdata/' + mpa + '.nc'
        nc_fid = Dataset(filename, 'r')
            
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
            
            distances300.append(vincenty(start_latlon, end_latlon).kilometers)
            
            i = i+1
        nc_fid.close()
        
    return distances300

# lophelia mpa subset for plotting
MPA_LOPHELIA = (['Anton Dohrn Seamount',
               'Darwin Mounds',
               'East Mingulay',
               'East Rockall Bank',
               'Faroe-Shetland Sponge Belt',
               'Hatton Bank',
               'North West Rockall Bank',
               'Rosemary Bank Seamount',
               'South-East Rockall Bank SAC',
               'Wyville Thomson Ridge',
               'Norwegian Boundary Sediment Plain',
               'Central Fladen',
               'South-West Porcupine Bank SAC',
               'Belgica Mound Province SAC','Hovland Mound Province SAC',
               'North-West Porcupine Bank SAC','Porcupine Bank Canyon SAC'])

MPA_SCOT =    (['Anton Dohrn Seamount',
               'Darwin Mounds',
               'East Mingulay',
               'East Rockall Bank',
               'Faroe-Shetland Sponge Belt',
               'Hatton Bank',
               'North West Rockall Bank',
               'Rosemary Bank Seamount',
               'South-East Rockall Bank SAC',
               'Wyville Thomson Ridge',
               'Norwegian Boundary Sediment Plain',
               'Central Fladen'])

MPA_RB =      (['Anton Dohrn Seamount',
               'East Rockall Bank',
               'Hatton Bank',
               'North West Rockall Bank',
               'South-East Rockall Bank SAC'])

MPA_NW =      (['Darwin Mounds',
               'Faroe-Shetland Sponge Belt',
               'Rosemary Bank Seamount',
               'Wyville Thomson Ridge'])

MPA_EM =      (['East Mingulay'])

MPA_NS =      (['Norwegian Boundary Sediment Plain',
               'Central Fladen'])
MPA_PORC =    (['South-West Porcupine Bank SAC',
               'Belgica Mound Province SAC','Hovland Mound Province SAC',
               'North-West Porcupine Bank SAC','Porcupine Bank Canyon SAC'])

#MPA_GROUPS = ([MPA_SCOT,MPA_RB,MPA_NW,MPA_EM])
#MPA_GROUPS = ([MPA_RB,MPA_EM,MPA_NW,MPA_NS])
MPA_GROUPS = ([MPA_NS])

RUN_NAME = ['behaviour2']
print MPA_GROUPS, RUN_NAME

icount = 0
                    
dist = []

for MPA_SOURCE in MPA_GROUPS:
   
    for year in range(1965,2005):
        year_str = str(year)
        print year
        
        # timesteps per day
        
        #startday = 32.0
        
        starttime = datetime.datetime(year,2,1)
        
        if platform.system() == 'Windows':
            run_dir = ('E:/af26/LarvalDispersalResults/'
                    + 'polcoms'+year_str+'/Run_1000_' + RUN_NAME[icount] + '/')
        elif platform.system() == 'Linux':
            run_dir = ('/home/af26/LarvalModelResults/Polcoms1990/Run_test/')
            
        # open the larval track datasets. Put them in a dictionary.
        dist1 = add_annual_data_to_stats(MPA_SOURCE, run_dir, RUN_NAME[icount])    
        
        dist = dist + dist1
        
print np.mean(dist), np.std(dist)
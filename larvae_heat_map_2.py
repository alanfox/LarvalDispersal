
"""
Created on Wed 24 June 2015

@author: af26

Reads track data from netCDF4 format files and produces an animation

"""
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
import shapefile
#from matplotlib.colors import ListedColormap, BoundaryNorm
from mpa_class import Mpa
from grid_class import Grid
import platform
import datetime as datetime


# for reading the MPA shapefiles. For drawing the mpa outlines.    

def read_shapefile(filename):
    sf = shapefile.Reader(filename)
    shapes = sf.shapes()
    records = sf.records()
    return shapes, records
            
    
# set up group of mpas
def read_mpas():
    mpa_group = set([])
    
    if platform.system() == 'Windows':    
        shapefile_root =   'G:/af26/Shapefiles/' #windows
    elif platform.system() == 'Linux':
        shapefile_root =   '/home/af26/Shapefiles/' #linux
    
    shapes, records = read_shapefile(shapefile_root + 
                    'UK_SAC_MAR_GIS_20130821b/UK_SAC_MAR_GIS_20130821b/' + 
                    'SCOTLAND_SAC_OFFSHORE_20121029_SIMPLE3')
    for i in range(len(shapes)):
        mpa_group.add(Mpa(shapes[i], records[i],'OFF_SAC'))
        
    # SAC with marine components
    shapes, records = read_shapefile(shapefile_root + 
                    'UK_SAC_MAR_GIS_20130821b/UK_SAC_MAR_GIS_20130821b/' + 
                    'SCOTLAND_SACs_withMarineComponents_20130821_SIMPLE3')
    for i in range(len(shapes)):
        mpa_group.add(Mpa(shapes[i], records[i],'MAR_SAC'))
        
    # Nature conservation MPA
    shapes, records = read_shapefile(shapefile_root + 
                    'MPA_SCOTLAND_ESRI/MPA_SCOTLAND_SIMPLE3')
    for i in range(len(shapes)):
        mpa_group.add(Mpa(shapes[i], records[i],'MPA'))
        
    # Irish SACs
    shapes, records = read_shapefile(shapefile_root + 
                    'SAC_ITM_WGS84_2015_01/SAC_Offshore_WGS84_2015_01')
    for i in range(len(shapes)):
        mpa_group.add(Mpa(shapes[i], records[i],'IRISH'))
        
    # Mikael Dahl's lophelia sites
    shapes, records = read_shapefile(shapefile_root + 
                    'MikaelDahl/MikaelDahl_1')
    for i in range(len(shapes)):
        mpa_group.add(Mpa(shapes[i], records[i],'Dahl'))
        
    return mpa_group
    

# INITIALISE STUFF
# need t grids for the 2d histogram

#nc_filedzu = 'C:/Users/af26/PolcommModelData/depth_dz_S12run420_UV.nc'
nc_filedzt = 'G:/af26/PolcommModelData/depth_dz_S12run420_TS.nc'
#nc_fiddzu = Dataset(nc_filedzu, 'r')
nc_fiddzt = Dataset(nc_filedzt, 'r')

    # read in and calculate the model grid variables
    
gridt = Grid(nc_fiddzt)
#gridu = Grid(nc_fiddzu)

nc_fiddzt.close()
#nc_fiddzu.close()

latitude = gridt.get_latitude()
longitude = gridt.get_longitude()

# read in GEBCO bathymetry

nc_infile = ('G:/af26/GEBCO/GEBCO_2014_2D_-20.0_40.0_13.0_65.0.nc')

nc_gebco = Dataset(nc_infile,'r')

bathymetry = nc_gebco.variables['elevation'][:]
lat_bath = nc_gebco.variables['lat'][:]    
lon_bath = nc_gebco.variables['lon'][:]    

nc_gebco.close()

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

#MPA_SOURCE = (['Anton Dohrn Seamount',
#               'Darwin Mounds',
#               'East Mingulay',
#               'East Rockall Bank',
#               'Faroe-Shetland Sponge Belt',
#               'Hatton Bank',
#               'North West Rockall Bank',
#               'Rosemary Bank Seamount',
#               'South-East Rockall Bank SAC',
#               'Wyville Thomson Ridge',
#               'Norwegian Boundary Sediment Plain',
#               'Central Fladen'])

MPA_SOURCE = (['Darwin Mounds',
               'Faroe-Shetland Sponge Belt',
               'Rosemary Bank Seamount',
               'Wyville Thomson Ridge'])


x = []
y = []
z = []

daysteps = 24.0
    
# minimum survivable temperature
    
min_temp = 0.0
# read the mpa shapefiles
  
mpa_group = read_mpas()    
     
for year in range(1965,2005):
    year_str = str(year)
    
    # timesteps per day
    
    #startday = 32.0
    
    starttime = datetime.datetime(year,2,1)
    
    if platform.system() == 'Windows':
        run_dir = ('F:/af26/LarvalDispersalResults/'
                + 'polcoms'+year_str+'/Run_1000_baseline/')
    elif platform.system() == 'Linux':
        run_dir = ('/home/af26/LarvalModelResults/Polcoms1990/Run_test/')
    
    # open the larval track datasets. Put them in a dictionary.
    
    nc_fid = {}
    for mpa in MPA_SOURCE:
        filename = run_dir + 'Trackdata/' + mpa + '.nc'
        nc_fid = Dataset(filename, 'r')
        
    
    # dead or alive
        temperature = nc_fid.variables['temperature'][:,:]
        not_dead = temperature > min_temp
        not_dead = np.cumprod(not_dead,axis = 1)
            
        releasedays = nc_fid.variables['release day'][:]
        
# loop over all larvae
        ii = 0
        for releaseday in releasedays:
            lat = nc_fid.variables['latitude'][ii,:]
            lon = nc_fid.variables['longitude'][ii,:]
            at_bed = nc_fid.variables['at bed'][ii,:]
            depth = nc_fid.variables['depth'][ii,:]
            temperature = nc_fid.variables['temperature'][ii,:]
            latdatarange = np.ma.flatnotmasked_edges(lat)
            for day in range(40,63):
#            for day in range(80,126):
                i = day * daysteps
                if (not_dead[ii,i]):
                    x.append(lon[i])
                    y.append(lat[i])
            ii += 1
        nc_fid.close()
            
plt.figure(figsize=(8,8))
for mpa in mpa_group:
#    print mpa.get_sitename()
    if mpa.get_sitename() in MPA_SOURCE:
        colour = 'red'
        mpa.plot_shape2(colour)                        
    elif mpa.get_sitename() in MPA_LOPHELIA:
        colour = 'blue'
        mpa.plot_shape2(colour)  
                              
hd = plt.hist2d(x,y,[longitude,latitude], cmin = 0.1)
plt.xlim([-20, 13])
plt.ylim([54, 65])
plt.colorbar(orientation='horizontal')
# save 300 m contour line
cs = plt.contour(lon_bath,lat_bath,bathymetry,[-1000,-200],
                 colors = 'grey',linestyles = 'solid')
plt.clabel(cs, inline=1, fontsize=8,fmt='%4.0f')                 
# colours land in black
plt.contourf(lon_bath,lat_bath,bathymetry,[0,5000], colors = 'lightgreen')
plt.contour(lon_bath,lat_bath,bathymetry,[0,5000], colors = 'black')
plt.ylabel('latitude')
plt.xlabel('longitude')
plt.show()

filenameroot = ('G:/Documents/Papers/LopheliaConnectivity/' + 
                 'heatmap_baseline_1965to2005_NorthWest')

plt.savefig(filenameroot + '.pdf')
plt.savefig(filenameroot + '.png', dpi = 1000)
#            
##    m.scatter(x,y, s = 5, marker = "o", c = z, 
##              cmap = cmap, norm = norm, edgecolor = 'none', latlon = True)
#    cs = m.scatter(x,y, s = 5, marker = "o", c = z, 
#              cmap = cmap, vmin = vmin, vmax = vmax, edgecolor = 'none', latlon = True)
#    cbar = m.colorbar(cs,location = 'bottom', pad = '7%')
#    cbar.ax.set_xlabel('depth of larva (m)',fontsize = 10) 
#    cbar.ax.tick_params(labelsize=8)
#       
#    plt.savefig(file_name,dpi = 80)
#    plt.clf()
#    
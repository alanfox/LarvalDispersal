
"""
Created on Wed 24 June 2015

@author: af26

Reads track data from netCDF4 format files and produces an animation

"""
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import shapefile
from matplotlib.colors import ListedColormap, BoundaryNorm
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
        shapefile_root =   'C:/Users/af26/Shapefiles/' #windows
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

nc_filedzu = 'C:/Users/af26/PolcommModelData/depth_dz_S12run420_UV.nc'
nc_filedzt = 'C:/Users/af26/PolcommModelData/depth_dz_S12run420_TS.nc'
nc_fiddzu = Dataset(nc_filedzu, 'r')
nc_fiddzt = Dataset(nc_filedzt, 'r')

    # read in and calculate the model grid variables
    
gridt = Grid(nc_fiddzt)
gridu = Grid(nc_fiddzu)

nc_fiddzt.close()
nc_fiddzu.close()

latitude = gridt.get_latitude()
longitude = gridt.get_longitude()

print len(latitude)
print len(longitude)

# read in GEBCO bathymetry

nc_infile = ('C:/Users/af26/GEBCO/GEBCO_2014_2D_-20.0_40.0_13.0_65.0.nc')

nc_gebco = Dataset(nc_infile,'r')

bathymetry = nc_gebco.variables['elevation'][:]
lat_bath = nc_gebco.variables['lat'][:]    
lon_bath = nc_gebco.variables['lon'][:]    

nc_gebco.close()


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

MPA_SOURCE = (['Hatton Bank'])

x = []
y = []
z = []
# and choose a year

for year in range(1965,2005):
    year_str = str(year)
    
    # timesteps per day
    
    #startday = 32.0
    
    starttime = datetime.datetime(year,2,1)
    
    daysteps = 24.0
    
    # minimum survivable temperature
    
    min_temp = 0.0
    
    if platform.system() == 'Windows':
        run_dir = ('E:/af26/LarvalDispersalResults/'
                + 'polcoms'+year_str+'/Run_1000_baseline/')
    elif platform.system() == 'Linux':
        run_dir = ('/home/af26/LarvalModelResults/Polcoms1990/Run_test/')
    
    # open the larval track datasets. Put them in a list.
    
    nc_fid = {}
    for mpa in MPA_SOURCE:
        filename = run_dir + 'Trackdata/' + mpa + '.nc'
        nc_fid[mpa] = Dataset(filename, 'r')
        
    # read the mpa shapefiles
      
    #mpa_group = read_mpas()    
        
    #width = 950000
    #height = 800000
    #
    #background = Basemap(width=width,height=height,projection='lcc',
    #            lat_0 = 59., lon_0 = -4.0, resolution='h')
    #-------------------------------------------------------------------
    # main plotting loop. Produces a series of plots. Compile into animation 
    # offline.
    #-------------------------------------------------------------------
    
    # dead or alive
    
    not_dead = {}
    for mpa in MPA_SOURCE:
        temperature = nc_fid[mpa].variables['temperature'][:,:]
        not_dead[mpa] = temperature > min_temp
        not_dead[mpa] = np.cumprod(not_dead[mpa],axis = 1)
    
        #m = background
        ## etopo is the topography/ bathymetry map background
        ##m.bluemarble()
        #m.drawcoastlines()
        #    
        #    # alternative plain map background
        #m.fillcontinents(color='black',lake_color='darkblue')
        #    
        #    # draw parallels and meridians.
        #m.drawparallels(np.arange(40.,66.,2.),labels = [1,0,0,0],fontsize = 8)
        #m.drawmeridians(np.arange(-24.,13.,4.),labels = [0,0,0,1],fontsize = 8)
        #    #m.drawmapboundary(fill_color='skyblue')
        #plt.title(timenow.strftime("%Y %b %d %H:%M"), loc = 'left')
        
        #for mpa in mpa_group:
        ##    print mpa.get_sitename()
        #    if mpa.get_sitename() in MPA_SOURCE:
        #        colour = 'red'
        #    else:
        #        colour = 'blue'
        #    mpa.plot_shape(background,colour)                        
                
    for mpa in MPA_SOURCE:
        releasedays = nc_fid[mpa].variables['release day'][:]
        
        ii = 0
        for releaseday in releasedays:
            lat = nc_fid[mpa].variables['latitude'][ii,:]
            lon = nc_fid[mpa].variables['longitude'][ii,:]
            at_bed = nc_fid[mpa].variables['at bed'][ii,:]
            depth = nc_fid[mpa].variables['depth'][ii,:]
            temperature = nc_fid[mpa].variables['temperature'][ii,:]
            latdatarange = np.ma.flatnotmasked_edges(lat)
            for day in range(50,63):
                i = day * daysteps
                if (not_dead[mpa][ii,i]):
                    x.append(lon[i])
                    y.append(lat[i])
            ii += 1
            
plt.figure()
plt.hist2d(x,y,[longitude,latitude], cmin = 0.1)
# save 300 m contour line
cs = plt.contour(lon_bath,lat_bath,bathymetry,[-1000,-300])
# colours land in black
plt.contourf(lon_bath,lat_bath,bathymetry,[0,5000], colors = 'lightgreen')
plt.contour(lon_bath,lat_bath,bathymetry,[0,5000], colors = 'black')
plt.ylabel('latitude')
plt.xlabel('longitude')
plt.show()
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
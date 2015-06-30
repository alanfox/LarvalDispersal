
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
import platform

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

# produce animation for larvae from a group of sites

MPA_SOURCE = ('East Mingulay',
              'Scanner Pockmark',
              'East Rockall Bank',
              'Faroe-Shetland Sponge Belt')

# and choose a year

year = 1995

# timesteps per day

daysteps = 24.0

# minimum survivable temperature

min_temp = 3.0

# plotting area constants

lllat = 54.0
lllon = -18.0
urlat = 65.0
urlon = 8.0


if platform.system() == 'Windows':
    run_dir = ('C:/Users/af26/LarvalDispersalResults/'
            + 'polcoms1995/Run_1000_baseline/')
elif platform.system() == 'Linux':
    run_dir = ('/home/af26/LarvalModelResults/Polcoms1990/Run_test/')

# open the larval track datasets. Put them in a list.

nc_fid = {}
for mpa in MPA_SOURCE:
    filename = run_dir + 'Trackdata/' + mpa + '.nc'
    nc_fid[mpa] = Dataset(filename, 'r')
    
# read the mpa shapefiles
  
mpa_group = read_mpas()    
    
# various colour maps, comment out as required
    
## black and red to show when at the bed
## Create a colormap for (black and red just to show when at bed):
#cmap = ListedColormap(['black', 'r'])
#
## Create a 'norm' (normalizing function) which maps data values to the interval [0,1]:
#norm = BoundaryNorm([-1000, 0.5, 1000], cmap.N)  # cmap.N is number of items in the colormap

# depth

cmap = 'autumn'
vmin = 0.0
vmax = 1000.0



# find length of animation in hours (this will be the number of frames)
# for this need to check start times and duration of tracks.

min_t = 1000.0
max_t = 0.0

for mpa in MPA_SOURCE:
    releasedays = nc_fid[mpa].variables['release day'][:]
    i = 0
    for releaseday in releasedays:
        x = nc_fid[mpa].variables['longitude'][i,:]
        # find indices of first and last unmasked values
        xdatarange = np.ma.flatnotmasked_edges(x)
        mini = releaseday * daysteps
        maxi = mini + xdatarange[1]
        min_t = min(min_t,mini)
        max_t = max(max_t,maxi)
        i = i + 1

#print min_t, max_t

n_steps = int(max_t-min_t)
min_release_day = int(min_t/daysteps)

#print n_steps, min_release_day

# for testing override nsteps
#n_steps = 24

background = Basemap(projection='lcc', llcrnrlat = lllat, llcrnrlon = lllon,
            urcrnrlat = urlat, urcrnrlon = urlon,
            lat_1 = 55., lon_0 = -4.0, resolution='h')

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

for i in range(n_steps):
    file_name = run_dir + 'temp_video/_temp%05d.png' % i
    
    m = background
    
    # etopo is the topography/ bathymetry map background
    m.etopo()
    m.drawcoastlines()
    
    # alternative plain map background
    #m.fillcontinents(color='coral',lake_color='skyblue')
    
    # draw parallels and meridians.
    m.drawparallels(np.arange(40.,66.,1.),labels = [1,1,0,0])
    m.drawmeridians(np.arange(-24.,13.,2.),labels = [0,0,0,1])
    m.drawmapboundary(fill_color='skyblue')
    #plt.title("POLCOMS model. Larval release through Feb " 
    #          + str(year) + ". Larsson et al behaviour")
    # draw the mpas

    for mpa in mpa_group:
    #    print mpa.get_sitename()
        if mpa.get_sitename() in MPA_SOURCE:
            colour = 'red'
        else:
            colour = 'blue'
        mpa.plot_shape(m,colour)
        
    x = []
    y = []
    z = []

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
            if (i >= int(releaseday * 24.0) and
                i <= int(releaseday * 24.0) + latdatarange[1] and
                not_dead[mpa][ii,i-int(releaseday * 24.0)]):
                x.append(lon[i-int(releaseday * 24.0)])
                y.append(lat[i-int(releaseday * 24.0)])
#                z.append(at_bed[i-int(releaseday * 24.0)])
                z.append(depth[i-int(releaseday * 24.0)])
            elif (i > int(releaseday * 24.0) + latdatarange[1] and
                    not_dead[mpa][ii,latdatarange[1]]):
                x.append(lon[latdatarange[1]])
                y.append(lat[latdatarange[1]])
#                z.append(at_bed[latdatarange[1]])
                z.append(depth[latdatarange[1]])
            ii += 1
                
#    m.scatter(x,y, s = 5, marker = "o", c = z, 
#              cmap = cmap, norm = norm, edgecolor = 'none', latlon = True)
    m.scatter(x,y, s = 5, marker = "o", c = z, 
              cmap = cmap, vmin = vmin, vmax = vmax, edgecolor = 'none', latlon = True)
        
    plt.savefig(file_name,dpi = 160)
    plt.clf()
    

"""
Created on Wed 11 Dec 10:00:00 2015

@author: af26

Reads track data from netCDF4 format file and plots the tracks on a map

"""
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import shapefile
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from mpa_class import Mpa

MPA_SOURCE = 'The Barra Fan and Hebrides Terrace Seamount'

year = 1993

nc_file = ('C:/Users/af26/Documents/LarvalDispersalResults/' +
           'polcoms1993/Run_BF300larvae_advect_linear/Trackdata/' + MPA_SOURCE + '.nc')
nc_fid = Dataset(nc_file, 'r')


def read_shapefile(filename):
    sf = shapefile.Reader(filename)
    shapes = sf.shapes()
    records = sf.records()
    return shapes, records
    
# plotting helper functions
        
        # Data manipulation:

def make_segments(x, y):
    '''
    Create list of line segments from x and y coordinates, in the correct format for LineCollection:
    an array of the form   numlines x (points per line) x 2 (x and y) array
    '''

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    return segments


# Interface to LineCollection:

def colorline(x, y, z=None, cmap=plt.get_cmap('copper'), norm=plt.Normalize(0.0, 1.0), linewidth=2, alpha=1.0):
    '''
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    '''
    
    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))
           
    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])
        
    z = np.asarray(z)
    
    segments = make_segments(x, y)
    lc = LineCollection(segments, array=z, cmap=cmap, norm=norm, linewidth=linewidth, alpha=alpha)
    
    ax = plt.gca()
    ax.add_collection(lc)
    
    return lc
        
    
def clear_frame(ax=None): 
    # Taken from a post by Tony S Yu
    if ax is None: 
        ax = plt.gca() 
    ax.xaxis.set_visible(False) 
    ax.yaxis.set_visible(False) 
    for spine in ax.spines.itervalues(): 
        spine.set_visible(False) 
                    

# set up group of mpas

mpa_group = set([])

# offshore SAC
shapes, records = read_shapefile('C:/Users/af26/Shapefiles/UK_SAC_MAR_GIS_20130821b/UK_SAC_MAR_GIS_20130821b/SCOTLAND_SAC_OFFSHORE_20121029_SIMPLE3')
for i in range(len(shapes)):
    mpa_group.add(Mpa(shapes[i], records[i],'OFF_SAC'))
    
# SAC with marine components
shapes, records = read_shapefile('C:/Users/af26/Shapefiles/UK_SAC_MAR_GIS_20130821b/UK_SAC_MAR_GIS_20130821b/SCOTLAND_SACs_withMarineComponents_20130821_SIMPLE3')
for i in range(len(shapes)):
    mpa_group.add(Mpa(shapes[i], records[i],'MAR_SAC'))
    
# Nature conservation MPA
shapes, records = read_shapefile('C:/Users/af26/Shapefiles/MPA_SCOTLAND_ESRI/MPA_SCOTLAND_SIMPLE3')
for i in range(len(shapes)):
    mpa_group.add(Mpa(shapes[i], records[i],'MPA'))
    
# Other lophelia areas
shapes, records = read_shapefile('C:/Users/af26/Shapefiles/SAC_ITM_WGS84_2015_01/SAC_Offshore_WGS84_2015_01')
for i in range(len(shapes)):
    mpa_group.add(Mpa(shapes[i], records[i],'IRISH'))
    
# Scandinavian sites
shapes, records = read_shapefile('C:/Users/af26/Shapefiles/MikaelDahl/Mikaeldahl_1')
for i in range(len(shapes)):
    mpa_group.add(Mpa(shapes[i], records[i],'Dahl'))
    


#---------------------------------------------------------------------
# plot the tracks
#---------------------------------------------------------------------

# Create a colormap for red, green and blue:
cmap = ListedColormap(['black', 'r'])

# Create a 'norm' (normalizing function) which maps data values to the interval [0,1]:
norm = BoundaryNorm([-1000, 0.5, 1000], cmap.N)  # cmap.N is number of items in the colormap


plt.figure()

m = Basemap(projection='lcc',llcrnrlat=40.,llcrnrlon=-18.,urcrnrlat=65.,\
            urcrnrlon=14.,lat_1=55.,lon_0 = -4.0,resolution='c')
                        
fates = nc_fid.variables['fate'][:]

i = 0
for fate in fates:
    col = 'black'
    x = nc_fid.variables['longitude'][i,:]
    y = nc_fid.variables['latitude'][i,:]
    z = nc_fid.variables['depth'][i,:]
    at_bed = nc_fid.variables['at bed'][i,:]
# Find the indices of the first and last unmasked values
    xdatarange = np.ma.flatnotmasked_edges(x)
    ydatarange = np.ma.flatnotmasked_edges(y)
    zdatarange = np.ma.flatnotmasked_edges(z)
    
    if z[zdatarange[0]] <= 5000.0:    
    
    # starting positions
        m.scatter(x[xdatarange[0]],y[ydatarange[0]], latlon = True, 
                  marker = "v", color = col)
    
    # for colorline need to remove masked values
        x = np.ma.compressed(x)
        y = np.ma.compressed(y)
        at_bed = np.ma.compressed(at_bed)
    
        x1,y1 = m(x,y) # plotting grid coordinates rather than lat/lon
        colorline(x1,y1,at_bed,cmap,norm)
       
    i = i+1
    
# etopo is the topography/ bathymetry map background
m.etopo()
m.drawcoastlines()

# alternative plain map background
#m.fillcontinents(color='coral',lake_color='skyblue')

# draw parallels and meridians.
m.drawparallels(np.arange(40.,66.,1.),labels = [1,1,0,0])
m.drawmeridians(np.arange(-24.,13.,2.),labels = [0,0,0,1])
#m.drawmapboundary(fill_color='skyblue')
plt.title("POLCOMS model. Larval release through Feb " 
          + str(year) + ". Larsson et al behaviour")
# draw the mpas

for mpa in mpa_group:
#    print mpa.get_sitename()
    if mpa.get_sitename() == MPA_SOURCE:
        colour = 'red'
    else:
        colour = 'blue'
    mpa.plot_shape(m,colour)
    print mpa.get_sitename()

#-------------------------------------------------------------------
# plot the start and end points
#-------------------------------------------------------------------

plt.figure()

m = Basemap(projection='lcc',llcrnrlat=40.,llcrnrlon=-18.,urcrnrlat=65.,\
            urcrnrlon=14.,lat_1=55.,lon_0 = -4.0,resolution='c')
                        
fates = nc_fid.variables['fate'][:]

i = 0
for fate in fates:
    col = '#000000'
    x = nc_fid.variables['longitude'][i,:]
    y = nc_fid.variables['latitude'][i,:]
    z = nc_fid.variables['depth'][i,:]

# Find the indices of the first and last unmasked values
    xdatarange = np.ma.flatnotmasked_edges(x)
    ydatarange = np.ma.flatnotmasked_edges(y)
    zdatarange = np.ma.flatnotmasked_edges(z)

    if z[zdatarange[0]] <= 5000.0:    
    # starting positions
        m.scatter(x[xdatarange[0]],y[ydatarange[0]], latlon = True, 
                  marker = "v", color = col)
        
    # ending positions       
        m.scatter(x[xdatarange[1]],y[ydatarange[1]], latlon = True, 
                  marker = "o", color = col)
        
    i = i+1
 
# etopo is the topography/ bathymetry map background
m.etopo()
m.drawcoastlines()

# alternative plain map background
#m.fillcontinents(color='coral',lake_color='skyblue')

# draw parallels and meridians.
m.drawparallels(np.arange(40.,66.,1.),labels = [1,1,0,0])
m.drawmeridians(np.arange(-24.,13.,2.),labels = [0,0,0,1])
m.drawmapboundary(fill_color='skyblue')
plt.title("POLCOMS model. Larval release through Feb " 
          + str(year) + ". Larsson et al behaviour")
# draw the mpas

for mpa in mpa_group:
#    print mpa.get_sitename()
    if mpa.get_sitename() == MPA_SOURCE:
        colour = 'red'
    else:
        colour = 'blue'
    mpa.plot_shape(m,colour)
#    print mpa.get_settled(), mpa.get_sitename()

##plt.savefig('foo.pdf')

#------------------------------------------------------------------------
# this draws a t,z plot of the position of the larvae in the water column
#------------------------------------------------------------------------

plt.figure()
plt.title("POLCOMS model. Larval release through Feb " 
          + str(year) + ". Larsson et al behaviour")
plt.ylabel('depth m')
plt.xlabel('time in days from release')
i = 0
for fate in fates:
    z = nc_fid.variables['depth'][i,:]
    rt = nc_fid.variables['release day'][i]
    at_bed = nc_fid.variables['at bed'][i,:]
    
# for colorline need to remove masked values
    z1 = np.ma.compressed(z)
    at_bed = np.ma.compressed(at_bed)
    t = np.array(range(len(z1)))
    t = rt + t/24.0 

    colorline(t,z1,at_bed,cmap,norm)
#    plt.plot(t,z)
    i = i+1
       
# set plot limits
ax = plt.gca()
ax.autoscale(True)
# reverse plot vertical axis
ymin, ymax = plt.ylim()
plt.ylim([ymax,ymin])

#------------------------------------------------------------------------
# this draws a temperature v time plot of the larvae
#------------------------------------------------------------------------

#plt.figure()
#plt.title("POLCOMS model. Larval release through Feb " 
#          + str(year) + ". Larsson et al behaviour")
#plt.ylabel('temperature')
#plt.xlabel('time in days from release')
#i = 0
#for fate in fates:
#    z = nc_fid.variables['temperature'][i,:]
#    rt = nc_fid.variables['release day'][i]
#    at_bed = nc_fid.variables['at bed'][i,:]
#    
## for colorline need to remove masked values
#    z1 = np.ma.compressed(z)
#    at_bed = np.ma.compressed(at_bed)
#    t = np.array(range(len(z1)))
#    t = rt + t/24.0 
#
#    colorline(t,z1,at_bed,cmap,norm)
##    plt.plot(t,z)
#    i = i+1
#       
## set plot limits
#ax = plt.gca()
#ax.autoscale(True)
## reverse plot vertical axis
#ymin, ymax = plt.ylim()
#plt.ylim([ymax,ymin])
#
#plt.show()
#
#nc_fid.close()

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
from boundary_class import Boundary
import platform


MPA_SOURCE = 'Wyville Thomson Ridge'

year = 1976

if platform.system() == 'Windows':
    run_dir = ('X:/Lophelia Share/Alan/LarvalDispersalResults/'
            + 'polcoms' + str(year) + '/Run_1000_doublelife/')
#    run_dir = ('E:/af26/LarvalDispersalResults/'
#            + 'polcoms1996/Run_1000_behaviour2/')
elif platform.system() == 'Linux':
    run_dir = ('/home/af26/LarvalModelResults/Polcoms1990/Run_test/')

nc_file = (run_dir +
           'Trackdata/' + MPA_SOURCE + '.nc')


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
    
    
#boundaries
    
# set up dictionary of boundaries to test crossing
nc_infile = ('C:/Users/af26/GEBCO/GEBCO_2014_2D_-20.0_40.0_13.0_65.0.nc')

nc_gebco = Dataset(nc_infile,'r')

elevation = nc_gebco.variables['elevation'][:]
lat = nc_gebco.variables['lat'][:]    
lon = nc_gebco.variables['lon'][:] 
   
# save 300 m contour line
cs = plt.contour(lon,lat,elevation,[-200])
# colours land in black

# extract points along the line
# index [0] in get_paths()[0] is a single continuous path if
# contour is broken.
# Need to look at result to check the right section is selected

p = cs.collections[0].get_paths()[0]
v = p.vertices
x = v[:,0]
y = v[:,1]

# cut this into sections to test
# crossings in each section.

x1  = [x[i] for i in range(len(x)) if x[i] < -6.2 
                                    and y[i] > 54.0 and y[i] < 57.94]
y1  = [y[i] for i in range(len(x)) if x[i] < -6.2 
                                    and y[i] > 54.0 and y[i] < 57.94]
x2  = [x[i] for i in range(len(x)) if x[i] < -6.2 
                                    and y[i] > 57.94]
y2  = [y[i] for i in range(len(x)) if x[i] < -6.2 
                                    and y[i] > 57.94]
x3  = [x[i] for i in range(len(x)) if x[i] >= -6.2 and x[i] < 1.4 
                                    and y[i] > 54.0]
y3  = [y[i] for i in range(len(x)) if x[i] >= -6.2 and x[i] < 1.4
                                    and y[i] > 54.0]
x4  = [x[i] for i in range(len(x)) if x[i] >= 1.4 
                                    and y[i] > 54.0 and y[i] < 62.0]
y4  = [y[i] for i in range(len(x)) if x[i] >= 1.4
                                    and y[i] > 54.0 and y[i] < 62.0]

a1 = [(-11.0,54.0),(-9.5,54.0)]
a2 = [(-9.21,57.94),(-8.57,57.81),(-7.29,57.60),
              (-6.35,57.41),(-5.2,57.3)]
a3 = [(-5.0,58.6),(-6.2,59.55)]
a4 = [(-3.4,58.5),(-3.0,59.0)]
a5 = [(-3.0,59.0),(-1.2,60.4)]
a6 = [(-1.2,60.4),(1.4,61.72)]
b1 = [(-11,54.0),(-20,54.0)] 
b2 = [(-9.21,57.94),(-13.9,58.9),(-20,58.9)]
b3 = [(-6.2,59.55),(-8.9,60.9),(-17.0,65.0)]
b4 = [(1.4,61.72),(-2.5,62.6),(-10.0,65.0)]
b5 = [(13.0,62.0),(4.07,62.0),(0.0,65.0)]
b6 = [(1.4,61.72),(4.07,62.0)]
c1 = [(-12.0,54.0),(-12,65.0)]
s1 = zip(x1,y1)
s2 = zip(x2,y2)
s3 = zip(x3,y3)
s4 = zip(x4,y4)

nc_gebco.close()

boundary_group = set([])
boundary_group.add(Boundary('A1',a1))
boundary_group.add(Boundary('A2',a2))
boundary_group.add(Boundary('A3',a3))
boundary_group.add(Boundary('A4',a4))
boundary_group.add(Boundary('A5',a5))
boundary_group.add(Boundary('A6',a6))
boundary_group.add(Boundary('B1',b1))
boundary_group.add(Boundary('B2',b2))
boundary_group.add(Boundary('B3',b3))
boundary_group.add(Boundary('B4',b4))
boundary_group.add(Boundary('B5',b5))
boundary_group.add(Boundary('B6',b6))
boundary_group.add(Boundary('C1',c1))
boundary_group.add(Boundary('S1',s1))
boundary_group.add(Boundary('S2',s2))
boundary_group.add(Boundary('S3',s3))
boundary_group.add(Boundary('S4',s4))
                                           

    
#---------------------------------------------------------------------
# plot the tracks
#---------------------------------------------------------------------

# Create a colormap for red, green and blue:
cmap = ListedColormap(['black', 'r'])

# Create a 'norm' (normalizing function) which maps data values to the interval [0,1]:
norm = BoundaryNorm([-1000, 0.5, 1000], cmap.N)  # cmap.N is number of items in the colormap

#lllat = 54.0
#lllon = -16.0
#urlat = 64.0
#urlon = 6.0

plt.figure()

#m = Basemap(projection='lcc', llcrnrlat = lllat, llcrnrlon = lllon,
#            urcrnrlat = urlat, urcrnrlon = urlon,
#            lat_1 = 55., lon_0 = -4.0, resolution='c')
                        
            
#width = 750000
#height = 900000
width = 1900000
height = 1600000

m = Basemap(width=width,height=height,projection='lcc',
            lat_0 = 58.0, lon_0 = -4.0, resolution='c')
            
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
#plt.title("POLCOMS model. Larval release through Feb " 
#          + str(year) + ". Larsson et al behaviour")
# draw the mpas

for mpa in mpa_group:
#    print mpa.get_sitename()
    if mpa.get_sitename() != 'Sea of the Hebrides':
        if mpa.get_sitename() == MPA_SOURCE:
            colour = 'red'
        else:
            colour = 'blue'
        mpa.plot_shape(m,colour)
        print mpa.get_sitename()

for boundary in boundary_group:
#    print mpa.get_sitename()
        boundary.plot_shape(m,'green')


#-------------------------------------------------------------------
# plot the start and end points
#-------------------------------------------------------------------

#plt.figure()
#
#m = Basemap(projection='lcc', llcrnrlat = lllat, llcrnrlon = lllon,
#            urcrnrlat = urlat, urcrnrlon = urlon,
#            lat_1 = 55., lon_0 = -4.0, resolution='c')
#                        
#fates = nc_fid.variables['fate'][:]
#
#i = 0
#for fate in fates:
#    col = '#000000'
#    x = nc_fid.variables['longitude'][i,:]
#    y = nc_fid.variables['latitude'][i,:]
#    z = nc_fid.variables['depth'][i,:]
#
## Find the indices of the first and last unmasked values
#    xdatarange = np.ma.flatnotmasked_edges(x)
#    ydatarange = np.ma.flatnotmasked_edges(y)
#    zdatarange = np.ma.flatnotmasked_edges(z)
#
#    if z[zdatarange[0]] <= 5000.0:    
#    # starting positions
#        m.scatter(x[xdatarange[0]],y[ydatarange[0]], latlon = True, 
#                  marker = "v", color = col)
#        
#    # ending positions       
#        m.scatter(x[xdatarange[1]],y[ydatarange[1]], latlon = True, 
#                  marker = "o", color = col)
#        
#    i = i+1
# 
## etopo is the topography/ bathymetry map background
#m.etopo()
#m.drawcoastlines()
#
## alternative plain map background
##m.fillcontinents(color='coral',lake_color='skyblue')
#
## draw parallels and meridians.
#m.drawparallels(np.arange(40.,66.,1.),labels = [1,1,0,0])
#m.drawmeridians(np.arange(-24.,13.,2.),labels = [0,0,0,1])
#m.drawmapboundary(fill_color='skyblue')
##plt.title("POLCOMS model. Larval release through Feb " 
##          + str(year) + ". Larsson et al behaviour")
## draw the mpas
#
#for mpa in mpa_group:
##    print mpa.get_sitename()
#    if mpa.get_sitename() == MPA_SOURCE:
#        colour = 'red'
#    else:
#        colour = 'blue'
#    mpa.plot_shape(m,colour)
##    print mpa.get_settled(), mpa.get_sitename()
#
###plt.savefig('foo.pdf')
#
##------------------------------------------------------------------------
## this draws a t,z plot of the position of the larvae in the water column
##------------------------------------------------------------------------
#
#plt.figure()
#plt.title("POLCOMS model. Larval release through Feb " 
#          + str(year) + ". Larsson et al behaviour")
#plt.ylabel('depth m')
#plt.xlabel('time in days from release')
#i = 0
#for fate in fates:
#    z = nc_fid.variables['depth'][i,:]
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
plt.show()
#
nc_fid.close()
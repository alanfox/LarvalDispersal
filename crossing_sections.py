# -*- coding: utf-8 -*-
"""


"""

from netCDF4 import Dataset
#import numpy as np
import matplotlib.pyplot as plt
import shapefile
#from matplotlib.colors import ListedColormap, BoundaryNorm
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
    
MPA_SOURCE = ([
'North-east Faroe-Shetland Channel',
'East Rockall Bank',
'Faroe-Shetland Sponge Belt',
'Anton Dohrn Seamount',
'Wyville Thomson Ridge',
'Darwin Mounds',
'East Mingulay',
'The Barra Fan and Hebrides Terrace Seamount',
'North West Rockall Bank',
'Geikie Slide and Hebridean Slope',
'Hatton Bank',
'Rosemary Bank Seamount',
'Hatton-Rockall Basin'
])

nc_infile = ('C:/Users/af26/GEBCO/GEBCO_2014_2D_-20.0_40.0_13.0_65.0.nc')

nc_gebco = Dataset(nc_infile,'r')

elevation = nc_gebco.variables['elevation'][:]
lat = nc_gebco.variables['lat'][:]    
lon = nc_gebco.variables['lon'][:] 

mpa_group = read_mpas()
   
# save 200 m contour line
cs = plt.contour(lon,lat,elevation,[-200],colors = 'k',linestyles = 'solid')
# colours land in black

# extract points along the line
# index [0] in get_paths()[0] is a single continuous path if
# contour is broken.
# Need to look at result to check the right section is selected

p = cs.collections[0].get_paths()[6]
v = p.vertices
x = v[:,0]
y = v[:,1]
points = zip(x,y)

# cut this into sections to test
# crossings in each section.

s1 = [z for z in points if z[0] < -6.2 
                                    and z[1] > 54.0 and z[1] < 57.94]
s2  = [z for z in points if z[0] < -6.2 
                                    and z[1] > 57.94]
s3  = [z for z in points if z[0] >= -6.2 and z[0] < 1.4 
                                    and z[1] > 54.0]
s4  = [z for z in points if z[0] >= 1.4 
                                    and z[1] > 54.0 and z[1] < 62.0]
                                                                        
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
#s1 = zip(x1,y1)
#s2 = zip(x2,y2)
#s3 = zip(x3,y3)
#s4 = zip(x4,y4)

nc_gebco.close()

plt.figure(figsize=(6,6.5))

plt.contourf(lon,lat,elevation,[0,5000], colors = 'lightgreen')
plt.contour(lon,lat,elevation,[0], colors = 'k')

cs1 = plt.contour(lon,lat,elevation,[-200,-300,-400,-500,-600,-700,-800,-1000,-2000,-4000],
                  colors = 'grey',linestyles = 'solid')
plt.clabel(cs1, inline=1, inline_spacing=-5, fontsize=6,fmt='%4.0f')                 
plt.xlim([-20, 3])
plt.ylim([50, 65])
#plt.contour(lon,lat,elevation,[-600],colors = 'k')
#plt.contour(lon,lat,elevation,[-1000],colors = 'k')

#plt.plot(x,y, color = 'darkgreen')

plt.plot(zip(*a1)[0],zip(*a1)[1], color = 'darkgreen')
plt.text((zip(*a1)[0][0]+zip(*a1)[0][-1])/2.0,
          (zip(*a1)[1][0]+zip(*a1)[1][-1])/2.0,
            'A1', color = 'darkgreen')
plt.plot(zip(*a2)[0],zip(*a2)[1], color = 'darkgreen')
plt.text((zip(*a2)[0][0]+zip(*a2)[0][-1])/2.0,
          (zip(*a2)[1][0]+zip(*a2)[1][-1])/2.0,
            'A2', color = 'darkgreen')
plt.plot(zip(*a3)[0],zip(*a3)[1], color = 'darkgreen')
plt.text((zip(*a3)[0][0]+zip(*a3)[0][-1])/2.0,
          (zip(*a3)[1][0]+zip(*a3)[1][-1])/2.0,
            'A3', color = 'darkgreen')
plt.plot(zip(*a4)[0],zip(*a4)[1], color = 'darkgreen')
plt.text((zip(*a4)[0][0]+zip(*a4)[0][-1])/2.0,
          (zip(*a4)[1][0]+zip(*a4)[1][-1])/2.0,
            'A4', color = 'darkgreen')
plt.plot(zip(*a5)[0],zip(*a5)[1], color = 'darkgreen')
plt.text((zip(*a5)[0][0]+zip(*a5)[0][-1])/2.0,
          (zip(*a5)[1][0]+zip(*a5)[1][-1])/2.0,
            'A5', color = 'darkgreen')
plt.plot(zip(*a6)[0],zip(*a6)[1], color = 'darkgreen')
plt.text((zip(*a6)[0][0]+zip(*a6)[0][-1])/2.0,
          (zip(*a6)[1][0]+zip(*a6)[1][-1])/2.0,
            'A6', color = 'darkgreen')
plt.plot(zip(*b1)[0],zip(*b1)[1], color = 'darkblue')
plt.text((zip(*b1)[0][0]+zip(*b1)[0][-1])/2.0,
          (zip(*b1)[1][0]+zip(*b1)[1][-1])/2.0,
            'B1', color = 'darkblue')
plt.plot(zip(*b2)[0],zip(*b2)[1], color = 'darkblue', linewidth = 3.0)
plt.text((zip(*b2)[0][0]+zip(*b2)[0][-1])/2.0 + 1.0,
          (zip(*b2)[1][0]+zip(*b2)[1][-1])/2.0 + 0.5,
            'B2', color = 'darkblue')
plt.plot(zip(*b3)[0],zip(*b3)[1], color = 'darkblue')
plt.text((zip(*b3)[0][0]+zip(*b3)[0][-1])/2.0,
          (zip(*b3)[1][0]+zip(*b3)[1][-1])/2.0,
            'B3', color = 'darkblue')
plt.plot(zip(*b4)[0],zip(*b4)[1], color = 'darkblue')
plt.text((zip(*b4)[0][0]+zip(*b4)[0][-1])/2.0,
          (zip(*b4)[1][0]+zip(*b4)[1][-1])/2.0,
            'B4', color = 'darkblue')
plt.plot(zip(*b5)[0],zip(*b5)[1], color = 'darkblue')
#plt.text((zip(*b5)[0][0]+zip(*b5)[0][-1])/2.0 -4.0,
#          (zip(*b5)[1][0]+zip(*b5)[1][-1])/2.0,
#            'B5', color = 'darkblue')
plt.plot(zip(*b6)[0],zip(*b6)[1], color = 'darkblue')
#plt.text((zip(*b6)[0][0]+zip(*b6)[0][-1])/2.0,
#          (zip(*b6)[1][0]+zip(*b6)[1][-1])/2.0,
#            'B6', color = 'darkblue')
plt.plot(zip(*c1)[0],zip(*c1)[1], color = 'coral', linewidth = 3.0)
plt.text((zip(*c1)[0][0]+zip(*c1)[0][-1])/2.0,
          (zip(*c1)[1][0]+zip(*c1)[1][-1])/2.0,
            'C1', color = 'coral')
plt.plot(zip(*s1)[0],zip(*s1)[1], color = 'r', linewidth = 3.0)
plt.text((zip(*s1)[0][0]+zip(*s1)[0][-1])/2.0,
          (zip(*s1)[1][0]+zip(*s1)[1][-1])/2.0,
            'S1', color = 'r')
plt.plot(zip(*s2)[0],zip(*s2)[1], color = 'r', linewidth = 3.0)
plt.text((zip(*s2)[0][0]+zip(*s2)[0][-1])/2.0,
          (zip(*s2)[1][0]+zip(*s2)[1][-1])/2.0,
            'S2', color = 'r')
plt.plot(zip(*s3)[0],zip(*s3)[1], color = 'r', linewidth = 3.0)
plt.text((zip(*s3)[0][0]+zip(*s3)[0][-1])/2.0,
          (zip(*s3)[1][0]+zip(*s3)[1][-1])/2.0,
            'S3', color = 'r')
plt.plot(zip(*s4)[0],zip(*s4)[1], color = 'r')
#plt.text((zip(*s4)[0][0]+zip(*s4)[0][-1])/2.0,
#          (zip(*s4)[1][0]+zip(*s4)[1][-1])/2.0-4.0,
#            'S4', color = 'r')


plt.ylabel('latitude')
plt.xlabel('longitude')

for mpa in mpa_group:
#    print mpa.get_sitename()
    if mpa.get_sitename() in MPA_SOURCE:
        colour = 'k'
        mpa.plot_shape2(colour)                        
#    elif mpa.get_sitename() in MPA_LOPHELIA:
#        colour = 'blue'
#        mpa.plot_shape2(colour)  


plt.show()

#filenameroot = ('G:/Documents/Papers/LopheliaConnectivity/' + 
#                 'crossings_sections_West')
#
#plt.savefig(filenameroot + '.pdf')
#plt.savefig(filenameroot + '.png', dpi = 1000)
                                   
                                    
                                    
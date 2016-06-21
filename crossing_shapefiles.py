# -*- coding: utf-8 -*-
"""

"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import shapefile

# for reading the MPA shapefiles. For drawing the mpa outlines.    

def read_shapefile(filename):
    sf = shapefile.Reader(filename)
    shapes = sf.shapes()
    records = sf.records()
    return shapes, records
            
# read in the topography
            
nc_infile = ('C:/Users/af26/GEBCO/GEBCO_2014_2D_-20.0_40.0_13.0_65.0.nc')

nc_gebco = Dataset(nc_infile,'r')

elevation = nc_gebco.variables['elevation'][:]
lat = nc_gebco.variables['lat'][:]    
lon = nc_gebco.variables['lon'][:] 
   
# using matplotlib routines to pull the contours out of the
# topography file. There must be a more elegant way.

# save 300 m contour line
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

l1 = [z for z in points if z[0] < -6.2 
                                    and z[1] > 54.0 and z[1] < 57.94]
l2  = [z for z in points if z[0] < -6.2 
                                    and z[1] > 57.94]
l3  = [z for z in points if z[0] >= -6.2 and z[0] < 1.4 
                                    and z[1] > 54.0]
l4  = [z for z in points if z[0] >= 1.4 
                                    and z[1] > 54.0 and z[1] < 62.0]


a1 = [[(-20,54.0),(-11.0,54.0),(-9.5,54.0)]]
a2 = [[(-20,58.9),(-13.9,58.9),(-9.21,57.94),(-8.57,57.81),(-7.29,57.60),
              (-6.67,57.48)]]
a3 = [[(-5.0,58.6),(-6.2,59.55),(-8.9,60.9),(-17.0,65.0)]]
a4 = [[(1.4,61.72),(-2.5,62.6),(-10.0,65.0)]]
a5 = [[(0.0,65.0),(4.07,62.0),(5.0,62.0)]]
b1 = [[(-3.3,58.63),(-3.0,59.0)]]
b2 = [[(-3.0,59.0),(-1.2,60.4)]]
b3 = [[(-1.2,60.4),(1.4,61.72)]]
b4 = [[(1.4,61.72),(4.07,62.0)]]
c1 = [[(-12.0,54.0),(-12,65.0)]]
s1 = [l3 + l2 + l1]
s2 = [l4]

nc_gebco.close()

w = shapefile.Writer(shapeType = 3)
w.field('NAME')
#w.autoBalance = 1
w.line(parts=a1)
w.record('A1')
w.line(parts=a2)
w.record('A2')
w.line(parts=a3)
w.record('A3')
w.line(parts=a4)
w.record('A4')
w.line(parts=a5)
w.record('A5')
w.line(parts=b1)
w.record('B1')
w.line(parts=b2)
w.record('B2')
w.line(parts=b3)
w.record('B3')
w.line(parts=b4)
w.record('B4')
w.line(parts=c1)
w.record('C1')
w.line(parts=s1)
w.record('S1')
w.line(parts=s2)
w.record('S2')
w.save('C:/Users/af26/Shapefiles/TransportSections/test')    
                             
prj = open("C:/Users/af26/Shapefiles/TransportSections/test.prj", "w") 
epsg = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]' 
prj.write(epsg) 
prj.close()                                   
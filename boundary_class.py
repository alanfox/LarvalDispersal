# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 14:41:58 2015

@author: af26

Mpa class definitions required for the polcoms particle tracking
model.

Needs to be imported as:

from mpa_class import Mpa



"""
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
import numpy as np

class Boundary:
    
    def __init__(self, name, points, area):
# points is a list of lon,lat pairs (tuples)        
        
        self.name = name        
        self.points = points
        self.area = area
        self.path = mplPath.Path(self.points)
        self.area_path = mplPath.Path(self.area)
        self.ncrossed = 0
        self.nstayed = 0
        
    def get_boundary_name(self):
        return self.name

    def get_points(self):
        return self.points
        
    def get_area(self):
        return self.area
               
    def get_crossed(self):
        return self.ncrossed
        
    def get_stayed(self):
        return self.nstayed
                
    def crosses(self,larva):
        # tests if an object of class Larva_tracks settles in the mpa
        lon = larva.get_lon()
        lat = larva.get_lat()
        life = len(lon)
# cumprod as once dead stays dead
        # first test that larva remains in viable temperature range
        temp = np.array(larva.get_temp())
        t_lower, t_upper = larva.get_temp_range()
        warm = temp >= t_lower
        cool = temp <= t_upper
        alive = np.cumprod(warm * cool)
# select points where larva are ready to settle, at the bed and alive
        larva_points = [(lon[i],lat[i]) for i in range(life) 
                                            if (alive[i])]

        if len(larva_points) == 0:
            return False   
        else:
            x = self.path.intersects_path(mplPath.Path(larva_points))
            if x:
                self.ncrossed = self.ncrossed + 1
                if self.area_path.contains_point(larva_points[-1]):
                    self.nstayed = self.nstayed + 1
                return True
        return False
                
    def plot_shape(self, m, colour):
        x = []
        y = []
#        print self.shape_path.contains_point((-14.0,58.0)), self.record[1]
        for point in self.points:
            x.append(point[0])
            y.append(point[1])
#        pat.Polygon(zip(x,y))          
        m.plot(x,y, latlon = True, color = colour)   

    def plot_shape2(self, colour):
        x = []
        y = []
#        print self.shape_path.contains_point((-14.0,58.0)), self.record[1]
        for point in self.points:
            x.append(point[0])
            y.append(point[1])
#        pat.Polygon(zip(x,y))          
        plt.plot(x,y, color = colour)   

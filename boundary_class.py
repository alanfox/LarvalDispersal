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
#import matplotlib.path as mplPath
#import numpy as np
from shapely.geometry import LineString


class Boundary:
    
    def __init__(self, name, points):
# points is a list of lon,lat pairs (tuples)    
        
        self.name = name        
        self.points = points
        self.path = LineString(points)
        self.ncrossed = 0
        self.nstayed = 0
        self.lon, self.lat = zip(*self.points)
        self.bbox = self.path.bounds
        
    def get_boundary_name(self):
        return self.name

    def get_points(self):
        return self.points
                       
    def get_crossed(self):
        return self.ncrossed
        
    def get_stayed(self):
        return self.nstayed
                
    def crosses(self,larva):
        # tests if an object of class Larva_tracks crosses the boundary
        # tests for returns by looking for odd or even numbers of crossings
        # odd number of crossings means stayed crossed
        # does not account for touching as unlikely 

        larva_live_points = larva.get_live_points()
        larva_live_points_bbox = larva.get_bbox_live_points()

        if len(larva_live_points) == 0:
            return False   
        elif self.bbox_overlap(self.bbox,larva_live_points_bbox):
            path_larva = LineString(larva_live_points)
            x = self.path.intersects(path_larva)
            if x:
                self.ncrossed = self.ncrossed + 1
                crossings = self.path.intersection(path_larva)
                try:
                    self.nstayed = self.nstayed + len(crossings)%2                    
                except TypeError:
                    self.nstayed = self.nstayed + 1
                return True
        return False

    def bbox_overlap(self,bbox1,bbox2):
        hoverlaps = (bbox1[0] <= bbox2[2]) and (bbox1[2] >= bbox2[0])
        voverlaps = (bbox1[3] >= bbox2[1]) and (bbox1[1] <= bbox2[3])
        return hoverlaps and voverlaps        

                
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

# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 14:41:58 2015

@author: af26

Mpa class definitions required for the polcoms particle tracking
model.

Needs to be imported as:

from mpa_class import Mpa



"""
import matplotlib.path as mplPath
from bngtolatlon import OSGB36toWGS84

class Mpa:
    
    def __init__(self, shape, record, area_type):
        
        
        self.shape = shape
        self.record = record
        self.area_type = area_type
        self.bbox = shape.bbox
        self.bbox_points = [[self.bbox[0],self.bbox[1]],
                            [self.bbox[2],self.bbox[1]],
                            [self.bbox[2],self.bbox[3]],
                            [self.bbox[0],self.bbox[3]]]
        self.points = shape.points
# convert from bng to latlon
        if (self.area_type == 'MPA' or self.area_type == 'MAR_SAC' 
                                    or self.area_type == 'SPA'):
            self.bbox = self.bng2lonlat_bbox(self.bbox)
            self.bbox_points = self.bng2lonlat(self.bbox_points)
            self.points = self.bng2lonlat(self.points)
            
        self.bbox_path = mplPath.Path(self.bbox_points)
        self.shape_path = mplPath.Path(self.points)
        self.nsettled = 0
        
    def get_bbox(self):
        return self.bbox
        
    def get_points(self):
        return self.points
        
    def get_sitename(self):
        if self.area_type == 'OFF_SAC':
            return self.record[1]
        if self.area_type == 'MPA':
            return self.record[0]
        if self.area_type == 'MAR_SAC':
            return self.record[0]
        if self.area_type == 'SPA':
            return self.record[0]
        if self.area_type == 'IRISH':
            return self.record[1]

        
    def get_settled(self):
        return self.nsettled
        
    def settles(self,larva):
        # tests if an object of class Larva settles in the mpa
        if larva.ready_to_settle():
            pos = larva.get_position()
            x = pos[0]
            y = pos[1]
            if self.bbox_path.contains_point((x,y)):
                if self.shape_path.contains_point((x,y)):
                    self.nsettled = self.nsettled + 1
                    return True
        return False
        
    def bng2lonlat(self,bng):
        lonlat = []
        for i in range(len(bng)):
            x = OSGB36toWGS84(bng[i][0],bng[i][1])
            y = [x[1],x[0]]
            lonlat.append(y)
        return lonlat
        
    def bng2lonlat_bbox(self,bng):
        ll = OSGB36toWGS84(bng[0],bng[1])
        ur = OSGB36toWGS84(bng[2],bng[3])
        bbox = [ll[1],ll[0],ur[1],ur[0]]
        return bbox
        
    def plot_shape(self, m, colour):
        x = []
        y = []
#        print self.shape_path.contains_point((-14.0,58.0)), self.record[1]
        for point in self.points:
            x.append(point[0])
            y.append(point[1])
#        pat.Polygon(zip(x,y))          
        m.plot(x,y, latlon = True, color = colour)   
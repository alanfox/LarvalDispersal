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
from bngtolatlon import OSGB36toWGS84
import numpy as np

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
                                    or self.area_type == 'SPA'
#                                    or self.area_type == 'Dahl'
                                    ):
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
        if self.area_type == 'Dahl':
            return self.record[1]
       
    def get_settled(self):
        return self.nsettled
        
    def settles(self,larva):
        # tests if an object of class Larva_tracks settles in the mpa
        settleage = larva.get_settleage()
        lon = larva.get_lon()
        lat = larva.get_lat()
        bed = larva.get_bed()
        life = len(lon)
        # first test that larva remains in viable temperature range
        # until settleage
        temp = larva.get_temp()
        t_lower, t_upper = larva.get_temp_range()
        for i in range(settleage):
            if ((temp[i] < t_lower) or (temp[i] > t_upper)):
                return False
        # then check for settling       
        for i in range(settleage,life):
            if ((temp[i] < t_lower) or (temp[i] > t_upper)):
                return False
            elif bed[i] == 1:
                if self.bbox_path.contains_point((lon[i],lat[i])):
#                    self.nsettled = self.nsettled + 1
#                    return True
                    if self.shape_path.contains_point((lon[i],lat[i])):
                        self.nsettled = self.nsettled + 1
                        return True
        return False
        
    def settles_2(self,larva):
#        # tests if an object of class Larva_tracks settles in the mpa
#        lon = larva.get_lon()
#        lat = larva.get_lat()
#        bed = larva.get_bed()
#        life = len(lon)
#        settleage = larva.get_settleage()
#        ofsettleage = np.array(range(life)) >= settleage
#        # first test that larva remains in viable temperature range
#        # until settleage
#        temp = np.array(larva.get_temp())
#        t_lower, t_upper = larva.get_temp_range()
#        warm = temp >= t_lower
#        cool = temp <= t_upper
## cumprod as once dead stays dead
#        alive = np.cumprod(warm * cool)
#        
## select points where larva are ready to settle, at the bed and alive
#        larva_bed_points = [(lon[i],lat[i]) for i in range(life) 
#                                            if (ofsettleage[i] and bed[i]==1
#                                            and alive[i])]

        larva_bed_points = larva.get_bed_points()
        larva_bed_points_bbox = larva.get_bbox_bed_points()

        if len(larva_bed_points) == 0:
            return False   
        elif self.bbox_overlap(self.bbox,larva_bed_points_bbox):
            x = self.shape_path.contains_points(larva_bed_points)
            if any(x):
                self.nsettled = self.nsettled + 1
                return True
        return False

    def bbox_overlap(self,bbox1,bbox2):
        hoverlaps = (bbox1[0] <= bbox2[2]) and (bbox1[2] >= bbox2[0])
        voverlaps = (bbox1[3] >= bbox2[1]) and (bbox1[1] <= bbox2[3])
        return hoverlaps and voverlaps        

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

    def plot_shape2(self, colour):
        x = []
        y = []
#        print self.shape_path.contains_point((-14.0,58.0)), self.record[1]
        for point in self.points:
            x.append(point[0])
            y.append(point[1])
#        pat.Polygon(zip(x,y))          
        plt.plot(x,y, color = colour)   

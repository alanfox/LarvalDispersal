# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 13:55:59 2015

@author: af26

Grid class definitions required for the polcoms particle tracking
model.

Needs to be imported as:

from Grid_class import Grid



"""
import numpy as np
import scipy.optimize as syop

    
class Grid:

    def __init__(self, nc_fid):
        self.nc_fid = nc_fid
        self.latitude = self.nc_fid.variables['latitude'][:]
        self.longitude = self.nc_fid.variables['longitude'][:]
        self.depth = self.nc_fid.variables['depth'][:]
        self.dz = self.nc_fid.variables['dz'][:]
        
# convert depths to more useful format of depth to layer boundary, rather than
# depth of mid-layer point.
        
        self.depth = self.depth + self.dz / 2.0        
        
        
    def get_latitude(self):
        return self.latitude
    
    def get_longitude(self):
        return self.longitude
    
    def get_depth(self):
        return self.depth
    
    def is_on_land(self, i, j):
        
        return np.ma.getmask(self.depth)[0, j, i]
        
    def is_on_bed(self, lon, lat, dep):
        
        d = self.get_total_depth_at_point(lon, lat)
        
        return dep >= d
        
    def get_depths_at_point(self, lon, lat):
        
        i, j = self.get_index_ne(lon, lat)
      
        # (bi)linear interpolation - not pretty but who cares
        #  well, i might because it uses a lot of time!
        
#        # depths at corners
#        d00 = self.depth[:, j-1, i-1]
#        d10 = self.depth[:, j-1, i  ]
#        d01 = self.depth[:, j  , i-1]
#        d11 = self.depth[:, j  , i  ]
#        
#        # interpolation coefficients        
#        dx1 = self.longitude[i] - lon
#        dx0 = lon - self.longitude[i-1]     
#        dy1 = self.latitude[j] - lat
#        dy0 = lat - self.latitude[j-1]
#        dxdy = ((self.longitude[i] - self.longitude[i-1])
#                * (self.latitude[j] - self.latitude[j-1]))
#                
#        # interpolate
#        return  (((d00 * dy1) + (d01 * dy0)) * dx1
#               + ((d10 * dy1) + (d11 * dy0)) * dx0) / dxdy
        
        N = 40
        B = 0.05
        theta = 8.0
        hc = 150.0
        
        hij = self.get_total_depth_at_point(lon,lat)
        if hij <= hc:
            z = np.linspace(hij/float(N),hij,N)
            return z
        else:
            Sk = np.linspace(-1.0,-1.0/float(N),N)
            h = (hij - hc)/hij
            Ck1 = (1.0 - B) * np.sinh(theta * Sk) / np.sinh(theta)
            Ck2 = (B * (np.tanh(theta * (Sk + 0.5)) - np.tanh(theta * 0.5)) / 
                       (2 * np.tanh(0.5* theta)))
                       
            Ck = Ck1 + Ck2
            z = Sk + h * (Ck - Sk)
       
            return - 1.0 * z * hij
                
    def get_total_depth_at_point(self, lon, lat):
        
        i, j = self.get_index_ne(lon, lat)
      
        # (bi)linear interpolation - not pretty but who cares
        
        # bottom depths at corners
        d00 = self.depth[0, j-1, i-1]
        d10 = self.depth[0, j-1, i  ]
        d01 = self.depth[0, j  , i-1]
        d11 = self.depth[0, j  , i  ]
        
        # interpolation coefficients        
        dx1 = self.longitude[i] - lon
        dx0 = lon - self.longitude[i-1]        
        dy1 = self.latitude[j] - lat
        dy0 = lat - self.latitude[j-1]
        dxdy = ((self.longitude[i] - self.longitude[i-1])
                * (self.latitude[j] - self.latitude[j-1]))
                
        # interpolate
        return  (((d00 * dy1) + (d01 * dy0)) * dx1
               + ((d10 * dy1) + (d11 * dy0)) * dx0) / dxdy
                
                
    def get_minlat(self):
        return self.latitude[0]
    
    def get_maxlat(self):
        return self.latitude[-1]
        
    def get_minlon(self):
        return self.longitude[0]
    
    def get_maxlon(self):
        return self.longitude[-1]
        
    def get_nx(self):
        return len(self.longitude)
        
    def get_ny(self):
        return len(self.latitude)
        
    def get_dx(self):
        # assumes regular grid
        return self.longitude[1] - self.longitude[0]
        
    def get_dy(self):
        # assumes regular grid
        return self.latitude[1] - self.latitude[0]
        
    def get_index_ne(self, lon, lat):
        # i coordinate point to NE of larva
        i = np.searchsorted(self.longitude, lon)
        # j coordinate point to NE of larva
        j = np.searchsorted(self.latitude, lat)
        
        return i,j
        
    def get_kindex(self,lon,lat,dep):
        
        d = self.get_depths_at_point(lon,lat)
        
        return np.searchsorted(-d, -dep) - 1

    def get_kindex_1(self,lon,lat,dep):
        
        # Modified level finder which works by finding the root of the
        # expression for s in terms of sigma. Slightly more expensive
        # but could be useful for tests of vertical interpolation I 
        # think (using sigma levels).
        
        def f_builder_1(s,hij):
            B = 0.05
            theta = 8.0
            hc = 150.0
            # builds the function finding the root of 0 = f(sigma) - s
            # needs to be built as it involves hij
            def f2(sigma):
                h = (hij - hc)/hij
                Ck1 = (1.0 - B) * np.sinh(theta * sigma) / np.sinh(theta)
                Ck2 = (B * (np.tanh(theta * (sigma + 0.5)) - np.tanh(theta * 0.5)) / 
                           (2 * np.tanh(0.5* theta)))
                           
                Ck = Ck1 + Ck2
                return (sigma + h * (Ck - sigma)) - s
            return f2

        N = 40    
        hc = 150.0
        hij = self.get_total_depth_at_point(lon,lat)
        if hij <= hc:
            return int(N - dep/(hij/N)), N - dep/(hij/N)
        else:
            s = -1.0 *  dep / hij
            f1 = f_builder_1(s, hij) 
            z = syop.brentq(f1,-1.1,0.1)
            
            
            return int(N * (z + 1.0)), N * (z + 1.0)



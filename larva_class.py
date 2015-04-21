# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 14:56:12 2015

@author: af26
Larva class definitions required for the polcoms particle tracking
model.

Needs to be imported as:

from larva_class import Larva

"""

import numpy as np
from scipy.interpolate import interp1d




class Larva:
    
    
    # each larva will be an instance of this class

    def __init__(self, pos, vel, source, release_day, gridt, 
                 RUN_CONST, SWIM_CONST):

        self.age = 0.0
        self.release_day = release_day
        self.at_bed = False
        self.pos = np.array([pos[0], pos[1], pos[2]])
        #if pos[2] is negative, bed release
        
        if self.pos[2] < 0.0:
            self.bed_release(self.pos, gridt)
            
        self.newpos = np.array([pos[0], pos[1], pos[2]])
        self.vel = np.array([vel[0], vel[1], vel[2]])
        self.xlon_history = []
        self.xlat_history = []
        self.depth_history = []
        self.xlon_history.append(self.pos[0])
        self.xlat_history.append(self.pos[1])
        self.depth_history.append(self.pos[2])
        self.turb = np.zeros((3),dtype = float)
        self.source = source
        self.i = 10000
        self.j = 10000
        self.fu2 = []
        self.fv2 = []
        self.fw2 = []
        
        #extract run constants
        
        self.seconds_in_day = RUN_CONST[0]
        self.m_to_degree = RUN_CONST[1]
        self.dt = RUN_CONST[2]
        self.km = RUN_CONST[3]
        self.vertical_interp = RUN_CONST[4]
        
        
        
        #larval behaviour constants - vary slightly for each larva
        self.swimslow = SWIM_CONST[0]
        self.swimfast = SWIM_CONST[1] 
        self.swimstart = SWIM_CONST[2]
        self.swimmax = SWIM_CONST[3] + np.random.normal(0.0,1.0)
        self.descendage = SWIM_CONST[4] + np.random.normal(0.0,1.0)
        self.fulldescendage = SWIM_CONST[5] + 2 * np.random.normal(0.0,1.0)
        self.minsettleage = SWIM_CONST[6]
        self.deadage = SWIM_CONST[7]
        
        
        

# uncomment to activate settling and dying
#        self.minsettleage = MINSETTLEAGE
#        self.deadage = DEADAGE
        
        
    def bed_release(self, position, gridt):
        
        # note, should not be called if the point is on land or outside the 
        # model region. These conditions are tested outside. Interpolates depth
        # the release point (linear interpolation from surrounding four T 
        # points)
        
        depth = gridt.get_total_depth_at_point(position[0], position[1])
                        
        self.pos[2] = depth - 5.0
        self.at_bed = True
        
        
    def get_position(self):
        return self.pos

    def get_source(self):
        return self.source
    
    def get_release_day(self):
        return self.release_day
    
    def get_track(self):
        return self.xlon_history, self.xlat_history
        
    def get_depth_history(self):
        return self.depth_history        
    
    def get_velocity(self):
        return self.vel
    
    def vertical_interpolation(self,i,j):
        # extend arrays above surface and below bed for interpolation
        # only does this when the larva has moved horizontally to a new 
        # box. Otherwise uses stored values. To reduce recalculation.
    
        if ((i != self.i) or (j != self.j)):
            self.i = i
            self.j = j
            nz = len(gridk[2])
            zlevs = np.insert(np.append(gridk[2],6000.0),0,-5.0)
            u1 = u[i,j]
            v1 = v[i,j]
            w1 = w[i,j]
            u2 = np.insert(np.append(u1,u1[nz-1]),0,u1[0])
            v2 = np.insert(np.append(v1,v1[nz-1]),0,v1[0])
            w2 = np.insert(np.append(w1,w1[nz-1]),0,w1[0])
        
            # cubic spline interpolation 
            self.fu2 = interp1d(zlevs, u2, kind = 'linear')
            self.fv2 = interp1d(zlevs, v2, kind = 'linear')
            self.fw2 = interp1d(zlevs, w2, kind = 'linear')

    def advection(self, gridt, u, v, w):
        
        # presently no interpolation in the horizontal direction
        # interpolation in the vertical
        # this was to save me coding time, reduce run-time
        
        # update velocity

        i, j = gridt.get_index_ne(self.pos[0], self.pos[1])
        
        if self.vertical_interp:
            # interpolate velocities in the vertical
#            self.vertical_interpolation(i,j)
#                  
#            self.vel[0] = self.fu2(self.pos[2])
#            self.vel[1] = self.fv2(self.pos[2])
#            self.vel[2] = self.fw2(self.pos[2])
            pass
            
        else:    
            # or just go with box larva is in
            # k box with larva in
            
            k = gridt.get_kindex(self.pos[0], self.pos[1], self.pos[2])
            # interpolate velocities
            # don't bother for now, just use the gridbox the point is in
            
            self.vel[0] = u[k,j,i]
            self.vel[1] = v[k,j,i]
            self.vel[2] = w[k,j,i]            
        

    def diffusion(self):
        # 2-d horizontal diffusion with constant k. Gaussian randon walk.
        # vertical diffusion, Gaussian random walk, different constant.
    
        rnorm = np.random.normal(0.0,1.0,3)
        self.turb = rnorm * np.sqrt(2.0 * self.km * self.dt)

    def vertical_behaviour(self):
        
        # implements swimming (or could be bouyancy) in the vertical.
        # days for each phase are taken from Larrson et al.
        # basically swims up at increasing rates to surface layer (5 m), 
        # stays there for a couple of weeks then heads back down to the bed
        
        # set percentage chance of heading upwards. If it doesn't go up it goes down.
        
        # set random swimming if in surface layer
        swim = np.zeros((3),dtype = float)                          
        if (self.pos[2] < 10.0):
            percent_up = 0.5
        else:
            if (self.age < self.swimstart):
                percent_up = 0.5
            elif (self.age < self.descendage):
                percent_up = 0.75
            elif (self.age > self.fulldescendage):
                percent_up = 0.25
            else:
                percent_up = 0.75 - 0.5 * ((self.age - self.descendage) 
                                      / (self.fulldescendage - self.descendage))
                                      
               
                
        if (self.age < self.swimstart):
            swimspeed = self.swimslow
        elif (self.age > self.swimmax):
            swimspeed = self.swimfast
        else:
            swimspeed = self.swimslow + ((self.swimfast - self.swimslow)
                                          * (self.age - self.swimstart)
                                          / (self.swimmax - self.swimstart))
                
#        print self.rundays, percent_up
        
        # set swim speed based on percent_up and swimspeed with some random noise
        
#        swimspeed = (swimspeed * 
#                    (float(np.random.choice([1,-1], 1, 
#                                        p=[percent_up,1.0 - percent_up ]))
#                     + np.random.normal(0.0,0.2)))
        swimspeed = (swimspeed * 
                    ((percent_up - (1.0 - percent_up))
                     + np.random.normal(0.0,0.25)))
        swim[2] = swimspeed                                
                                        
        return swim              
                       
    def isoutofarea(self, position, gridt):
        
        i, j = gridt.get_index_ne(position[0], position[1])
        nx = gridt.get_nx()
        ny = gridt.get_ny()

        if i == 0 or i == nx:
            return True
        if j == 0 or j == ny:
            return True
            
        return False
        
    def isonland(self, position, gridu, gridt, u):
        
        i, j = gridt.get_index_ne(position[0], position[1])
        
        if gridu.is_on_land(i,j):
            return True
        elif np.ma.getmask(u)[0,j,i]:
            return True
        elif gridt.is_on_bed(position[0],position[1],position[2]):
            return True
        return False

    def ready_to_settle(self):
        return ((self.age > self.minsettleage) and self.at_bed)
        
    def dead(self):
        return (self.age > self.deadage)
                        
                
    def update(self, dt, rundays, gridu, gridt, u, v, w):
        
        
        self.rundays = rundays
        self.age = self.age + dt / self.seconds_in_day
        self.at_bed = False
        
        # updates the larva position, returns a boolean value - True if 
        # the larva has hit the bed or left the area, otherwise False
        
        # update position
        m_to_degree_lon = self.m_to_degree / np.cos(np.radians(self.pos[1]))
        m_to_degree = np.array([m_to_degree_lon, self.m_to_degree, 1.0])
                      
        self.isonland(self.pos, gridu, gridt, u)

        # advection
        self.advection(gridt, u, v, w)
        self.newpos = self.pos + self.vel * m_to_degree * self.dt
        
#        print 'advection', self.newpos
        
        # test if advected out of area, if so end
        
        if self.isoutofarea(self.newpos, gridt):
            self.pos = self.newpos
            self.xlon_history.append(self.pos[0])
            self.xlat_history.append(self.pos[1])
            self.depth_history.append(self.pos[2])
            return True
            
        # test if advected on to land
        # if it is don't advect
            
        if self.isonland(self.newpos, gridu, gridt, u):
            self.newpos = self.pos
            self.at_bed = True
#            print 'on land', self.newpos
        
        # lock in advection update
        self.pos = self.newpos
#        print 'after advection', self.pos

        # vertical swimming
            
        self.newpos = self.pos - self.dt * self.vertical_behaviour()
        
#        print 'swimming', self.newpos, self.pos
        
        # test if swims on to bed
        # if it is don't swim

        if self.isonland(self.newpos, gridu, gridt, u):
            self.newpos = self.pos
            self.at_bed = True
#            print 'on land', self.newpos 
        
        # test if swims out of surface
        # if it is don't swim

        if (self.newpos[2] < 0.0):
            self.newpos = self.pos
            
        # lock in swimming update
        self.pos = self.newpos
#        print 'after swimming', self.pos
            
                
        # diffusion
        self.diffusion()
        self.newpos = self.newpos + self.turb * m_to_degree
        
#        print 'diffusion', self.newpos
        # test if diffused out of area, if so end
        
        if self.isoutofarea(self.newpos, gridt):
            self.pos = self.newpos
            self.xlon_history.append(self.pos[0])
            self.xlat_history.append(self.pos[1])
            self.depth_history.append(self.pos[2])
            return True
            
        # test if diffused on to land
        # if it is don't diffuse
        
        if self.isonland(self.newpos, gridu, gridt, u):
            self.newpos = self.pos
            self.at_bed = True
#            print 'on land', self.newpos 
            
        # test if diffused out of surface
        # if it is, no vertical diffusion

        if (self.newpos[2] < 0.0):
            self.newpos[2] = self.pos[2]
                        
        # lock in diffusion update
        self.pos = self.newpos
                
        # store the new position in track history
        
        self.xlon_history.append(self.pos[0])
        self.xlat_history.append(self.pos[1])
        self.depth_history.append(self.pos[2])
        
        return False
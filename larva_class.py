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

    def __init__(self, pos, vel, source, release_day, gridt, gridu,
                 RUN_CONST, SWIM_CONST, temperature, salinity, T_LOWER, T_UPPER):

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
        self.bed_history = []
        self.temperature_history = []
        self.salinity_history = []
        self.xlon_history.append(self.pos[0])
        self.xlat_history.append(self.pos[1])
        self.depth_history.append(self.pos[2])
        self.bed_history.append(0)        
        self.turb = np.zeros((3),dtype = float)
        self.source = source
        self.i = 10000
        self.j = 10000
        # nearest u gridpoint
        self.ipos = 10000
        self.jpos = 10000
        self.kpos = 10000
        self.kpos_real = 10000.00
        # nearest t gridpoint
        self.ipost = 10000
        self.jpost = 10000
        
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
        # random element means larvae reach milestones at different times
        self.swimmax = SWIM_CONST[3] + np.random.normal(0.0,1.0)
        self.descendagerange = SWIM_CONST[5]
        self.descendage = SWIM_CONST[4] + np.random.normal(0.0,
                                                self.descendagerange)
        self.minsettleage = SWIM_CONST[6]
        self.deadage = SWIM_CONST[7]
        self.targetdepth = SWIM_CONST[8]
        self.behaviour = SWIM_CONST[9]
        self.t_lower = T_LOWER
        self.t_upper = T_UPPER
                
        # monitor health
        self.isdead = False
        
        # find initial position in grid - updates ipos and ipost etc
        self.update_kji(gridt,gridu)
        
        self.temperature_history.append(
                        temperature[self.kpos,self.jpost,self.ipost])
        self.salinity_history.append(
                        salinity[self.kpos,self.jpost,self.ipost])

        
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
        return self.depth_history, self.bed_history        
    
    def get_temperature_history(self):
        return self.temperature_history       
    
    def get_salinity_history(self):
        return self.salinity_history       
    
    def get_velocity(self):
        return self.vel
    
    def vertical_interpolation(self,i,j,u,v,w):
        # extend arrays above surface and below bed for interpolation
        # interpolation done in sigma space, where velocities stored at 
        # 0.5,1.5,2.5.....39.5
    
#        if ((i != self.i) or (j != self.j)):
        self.i = i
        self.j = j
        u1 = u[:,j,i]
        v1 = v[:,j,i]
        w1 = w[:,j,i]
        nz = len(u1)
        zlevs = np.linspace(-0.5,40.5,42)
        u2 = np.insert(np.append(u1,u1[nz-1]),0,u1[0])
        v2 = np.insert(np.append(v1,v1[nz-1]),0,v1[0])
        w2 = np.insert(np.append(w1,w1[nz-1]),0,w1[0])
            
        # linear interpolation 
        self.fu2 = interp1d(zlevs, u2, kind = 'linear')
        self.fv2 = interp1d(zlevs, v2, kind = 'linear')
        self.fw2 = interp1d(zlevs, w2, kind = 'linear')

    def update_kji(self,gridt,gridu):
        
        # find the position on the grid. Moved here (out of 'advection')
        # because it is needed for assessing temperature and for
        # advection. Finding k is expensive
        
        self.ipos, self.jpos = gridt.get_index_ne(self.pos[0], self.pos[1])
        self.ipost, self.jpost = gridu.get_index_ne(self.pos[0], self.pos[1])
        self.ipost = self.ipost - 1
        self.jpost = self.jpost - 1
# two options for get_kindex. _1 version is needed for vertical
# interpolation but is slower than standard version if no vertical
# interpolation
        if self.vertical_interp:
            self.kpos, self.kpos_real = gridt.get_kindex_1(
                                self.pos[0], self.pos[1], self.pos[2])
        else:
            self.kpos = gridt.get_kindex(self.pos[0], self.pos[1], self.pos[2])
        
        return

    def advection(self, gridt, u, v, w):
        
        # presently no interpolation in the horizontal direction
        # interpolation in the vertical not implemented either
        # this was to save me coding time, reduce run-time 
        
        # update velocity

        i = self.ipos
        j = self.jpos
        
        if self.vertical_interp:
            # interpolate velocities in the vertical
            self.vertical_interpolation(i,j,u,v,w)
            
            k = self.kpos
                                          
            self.vel[0] = self.fu2(self.kpos_real)
            self.vel[1] = self.fv2(self.kpos_real)
            self.vel[2] = self.fw2(self.kpos_real)
            
        else:    
            # or just go with box larva is in
            # k box with larva in
            
            k = self.kpos
            
            self.vel[0] = u[k,j,i]
            self.vel[1] = v[k,j,i]
            self.vel[2] = w[k,j,i]

        # check if larva is dead    
        

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
        
        # set random swimming if at target depth
        swim = np.zeros((3),dtype = float)                          
        percent_up = 0.5
        factor = self.at_target()
        
        if (self.age > self.swimstart and self.age < self.descendage):
            percent_up = percent_up + factor * 0.5
        else:
            percent_up = 0.0
                                               
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
                     + np.random.normal(0.0,0.5)))
        swim[2] = swimspeed                                
                                        
        return swim
              
    def vertical_behaviour_2(self):
        
        # implements swimming (or could be bouyancy) in the vertical.
        # days for each phase are taken from Larrson et al.

        # new vertical behaviour function
        
        # considers vertical behaviour to consist of a gaussian random 
        # component dependent on the swimspeed plus a directional component
        # dependent on swimspeed and age and relation to target depth or 
        # density

        swim = np.zeros((3),dtype = float)                          
        
        if (self.age < self.swimstart):
            swimspeed = self.swimslow
        elif (self.age > self.swimmax):
            swimspeed = self.swimfast
        else:
            swimspeed = self.swimslow + ((self.swimfast - self.swimslow)
                                          * (self.age - self.swimstart)
                                          / (self.swimmax - self.swimstart))
        # 'diffusion' coeff proportional to the square of the swim speed
        # 3.0 is determined by simple 1d diffusion model. This is the random
        #  part of the motion
                                          
        kh = (swimspeed / 2.0)**2
        fact = np.sqrt(2.0 * kh * self.dt)
        zdisp = np.random.normal(0.0, 1.0) * fact
        
        # 'flow' set to a proportion of max swimspeed. This represents the
        # directional part of the motion. Let it fall faster than it rose,
        # or it might not get back to the bed
        if (self.age > self.swimstart and self.age < self.descendage):
            flow = self.swimfast / 50.0
        else:
            flow = -5.0 * self.swimfast / 50.0
        disp = flow * self.dt

        zdisp = zdisp + disp
            
        swim[2] = zdisp  
#        print swim                             
                                        
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
        return (self.isdead or (self.age > self.deadage))
        
    def at_target(self):
        # tests depth against target. Returns 1 if below target depth, 
        # -1 if above and 0 if at target.
        # set tolerance to 10 m for tight depth control
        # set to 5000 m for random swimming before descent 
        # should put this in the input file
        
        if self.pos[2] > self.targetdepth + 10.0:
            return 1.0
        elif self.pos[2] < self.targetdepth - 10.0:
            return -1.0
        else:
            return 0.0
                                        
    def update(self, dt, rundays, gridu, gridt, u, v, w, 
               temperature, salinity):
        
        self.rundays = rundays
        self.age = self.age + dt / self.seconds_in_day
        self.at_bed = False
        
        # updates the larva position, returns a boolean value - True if 
        # the larva has hit the bed or left the area, otherwise False
        
        # update position
        m_to_degree_lon = self.m_to_degree / np.cos(np.radians(self.pos[1]))
        m_to_degree = np.array([m_to_degree_lon, self.m_to_degree, 1.0])
                      
#        self.isonland(self.pos, gridu, gridt, u)
                      
        
        # check for heat death
        
#        print self.kpos,self.jpos,self.ipos
#        print temperature[self.kpos,self.jpos,self.ipos]
        
        self.isdead = (
            (temperature[self.kpos,self.jpos,self.ipos] > self.t_upper 
            or temperature[self.kpos,self.jpos,self.ipos] < self.t_lower)
            )
                      
        if self.isdead:
            return False

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
            self.bed_history.append(0)
            # out of grid, just repeat final temperature
            self.temperature_history.append(self.temperature_history[-1])
            self.salinity_history.append(self.salinity_history[-1])
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
        
        if (self.behaviour == 1):        
            self.newpos = self.pos - self.dt * self.vertical_behaviour()
        elif (self.behaviour == 2):
            self.newpos = self.pos - self.vertical_behaviour_2()
        
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
            self.bed_history.append(0)
            # out of grid, just repeat final temperature
            self.temperature_history.append(self.temperature_history[-1])
            self.salinity_history.append(self.salinity_history[-1])
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
        if self.at_bed:
            self.bed_history.append(1)
        else:
            self.bed_history.append(0)
        # find position in grid
        self.update_kji(gridt,gridu)
        self.temperature_history.append(
                    temperature[self.kpos,self.jpos,self.ipos])
        self.salinity_history.append(
                    salinity[self.kpos,self.jpos,self.ipos])
        
        return False
        
class Larval_tracks:
    
    def __init__(self, lon, lat, dep, bed, rt, fate, temp, sal,
                 MPA_SOURCE, RUN_CONST, SWIM_CONST, 
                 T_LOWER, T_UPPER):
                     
        self.lon = lon
        self.lat = lat
        self.dep = dep
        self.bed = bed
        self.rt = rt
        self.fate = fate
        self.temp = temp
        self.sal = sal
        self.source = MPA_SOURCE
        
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
        # random element means larvae reach milestones at different times
        self.swimmax = SWIM_CONST[3] + np.random.normal(0.0,1.0)
        self.descendagerange = SWIM_CONST[5]
        self.descendage = SWIM_CONST[4] + np.random.normal(0.0,
                                                self.descendagerange)
        self.minsettleage = SWIM_CONST[6]
        self.deadage = SWIM_CONST[7]
        self.targetdepth = SWIM_CONST[8]
        self.behaviour = SWIM_CONST[9]
        self.t_lower = T_LOWER
        self.t_upper = T_UPPER
        
        #calculate and store parameters for use later
        
        self.points = zip(lon,lat)
        life = len(self.lon)
        self.settleage = int(self.minsettleage * self.seconds_in_day / self.dt)
        ofsettleage = np.array(range(life)) >= self.settleage
        # first test that larva remains in viable temperature range
        # until settleage
        temp = np.array(self.temp)
        warm = temp >= self.t_lower
        cool = temp <= self.t_upper
# cumprod as once dead stays dead
        alive = np.cumprod(warm * cool)

# select points where larvae are alive
        self.live_points = [(self.lon[i],self.lat[i]) for i in range(life) 
                                            if (alive[i])]
        
# select points where larva are ready to settle, at the bed and alive
        self.bed_points = [(self.lon[i],self.lat[i]) for i in range(life) 
                                            if (ofsettleage[i] and bed[i]==1
                                            and alive[i])]
                                            
# find boundary boxes                                            
        self.bbox_points = [min(lon),min(lat),max(lon),max(lat)]
        if len(self.bed_points) != 0:
            x,y = zip(*self.bed_points)
            self.bbox_bed_points = [min(x),min(y),max(x),max(y)]
        else:
            self.bbox_bed_points = []                 
        if len(self.live_points) != 0:
            x,y = zip(*self.live_points)
            self.bbox_live_points = [min(x),min(y),max(x),max(y)]
        else:
            self.bbox_live_points = []                 
        
        

    def get_lon(self):
        return self.lon
        
    def get_lat(self):
        return self.lat
        
    def get_bed(self):
        return self.bed
        
    def get_temp(self):
        return self.temp
        
    def get_sal(self):
        return self.sal
        
    def get_temp_range(self):
        return self.t_lower, self.t_upper
             
    def get_settleage(self):
        return self.settleage

    def get_points(self):
        return self.points

    def get_bed_points(self):
        return self.bed_points

    def get_live_points(self):
        return self.live_points
                
    def get_bbox_points(self):
        return self.bbox_points
        
    def get_bbox_bed_points(self):
        return self.bbox_bed_points
        
    def get_bbox_live_points(self):
        return self.bbox_live_points
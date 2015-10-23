
"""

Takes netcdf files output by larval_dispersal_polcoms.py and tests whether
they cross boundaries. For example into the North Sea.

"""
from netCDF4 import Dataset
import numpy as np
import shapefile
import networkx as NX
from boundary_class import Boundary
from larva_class import Larval_tracks
import platform

if platform.system() == 'Windows':
    run_dir = ('C:/Users/af26/LarvalDispersalResults/'
            + 'polcoms2004/Run_1000_behaviour2/')
elif platform.system() == 'Linux':
    run_dir = ('/home/af26/LarvalModelResults/Polcoms1990/Run_test/')

graph_output_dir = run_dir + 'Crossingdata/'
track_input_dir = run_dir + 'Trackdata/'
mpa_name_file = open(run_dir + 'MPA_names_lophelia.txt', 'r') 

input_data_file = open(run_dir + 'input.dat', 'r')

input_dict = {}

for line in input_data_file:
    wordlist = line.split()
    input_dict[wordlist[0]] = wordlist[-1]
        
# larvae are released from MPA_SOURCE

# NUM_LARVAE are released at the start of each day for RELEASE_WINDOW days
# for a single release at time zero set RELEASE_WINDOW negative

NUM_LARVAE = int(input_dict['NUM_LARVAE'])
RELEASE_WINDOW = int(input_dict['RELEASE_WINDOW'])

STARTDAY = int(input_dict['STARTDAY'])

SECONDS_IN_DAY = 60.0 * 60.0 * 24.0
RADIUS_OF_EARTH = 6378160.0
M_TO_DEGREE = 360.0 / (2.0 * np.pi * RADIUS_OF_EARTH)

DT = float(input_dict['DT'])
DTDAYS = DT / SECONDS_IN_DAY

KMX = float(input_dict['KMX'])
KMY = float(input_dict['KMY'])
KMZ = float(input_dict['KMZ'])

KM = np.array([KMX, KMY, KMZ]) #constant diffusion coefficient m2/s

VERTICAL_INTERP = bool(int(input_dict['VERTICAL_INTERP']))
ANIMATE = bool(int(input_dict['ANIMATE']))
DEATH = bool(int(input_dict['DEATH']))

# constants for larval behaviour
# larvae released at bed, head upwards with increasing swim speeds up to 
# average age SWIMMAX. Then swimming gradually directed more downwards from 
# age DESCENDAGE up to DEADAGE.
BEHAVIOUR = int(input_dict['BEHAVIOUR'])    # behaviour index
# 1 is standard Larrson et al behaviour
# 2 is modified Larsson et al with small swimming speeds
TARGETDEPTH = float(input_dict['TARGETDEPTH'])    # target depth
SWIMSLOW = float(input_dict['SWIMSLOW'])          #initial swimming speed
SWIMFAST = float(input_dict['SWIMFAST'])      # max swimming speed
SWIMSTART = float(input_dict['SWIMSTART'])         #age in days at which start swimming
SWIMMAX = float(input_dict['SWIMMAX'])          #average age in days at which max swimming speed is reached
DESCENDAGE = float(input_dict['DESCENDAGE'])       # average age at which probability of heading down starts
                        # to increase
DESCENDAGERANGE = float(input_dict['DESCENDAGERANGE'])    # now fully heading down
MINSETTLEAGE = float(input_dict['MINSETTLEAGE'])     # minimum age at which can settle given suitable 
                        # habitat
DEADAGE = float(input_dict['DEADAGE'])          # Average age at which dead
# just set DEADAGE to a large value if larvae are not dying
if not DEATH:
    DEADAGE = 1000.0

if DEATH:
    if RELEASE_WINDOW > 0:
        NRUNDAYS = DEADAGE + RELEASE_WINDOW + 1
    else:
        NRUNDAYS = DEADAGE + 1
else:
    NRUNDAYS = float(input_dict['NRUNDAYS'])
    
# set temperature bounds for life
# just set wide here, all larvae survive but temperature
# along track is recorded for possible use later

T_LOWER = -10.0
T_UPPER = 100.0     
    
    
# bring constants together for passing to larva class

RUN_CONST = [SECONDS_IN_DAY, M_TO_DEGREE, DT, KM, VERTICAL_INTERP]
SWIM_CONST = [SWIMSLOW,SWIMFAST,SWIMSTART,SWIMMAX,
              DESCENDAGE,DESCENDAGERANGE,
              MINSETTLEAGE,DEADAGE,TARGETDEPTH, BEHAVIOUR]
    
def read_shapefile(filename):
    sf = shapefile.Reader(filename)
    shapes = sf.shapes()
    records = sf.records()
    return shapes, records

# helper functions
        
def group_cross(boundary_sprite_group, larva_object):
    crossed = False
    for boundary in set(boundary_sprite_group):
        if boundary.crosses(larva_object):
            crossed = True
    return crossed

def group_group_cross(larval_sprite_group, boundary_sprite_group):
    for larva in set(larval_sprite_group):
        group_cross(boundary_sprite_group, larva)
              
G = NX.DiGraph()
   
# loop over larval track files
   
for line in mpa_name_file:
    if platform.system() == 'Windows':    
        shapefile_root =   'C:/Users/af26/Shapefiles/' #windows
    elif platform.system() == 'Linux':
        shapefile_root =   '/home/af26/Shapefiles/' #linux
    
    MPA_SOURCE = line.rstrip()
    print MPA_SOURCE

# area checking entry to    
    NorthSea = [(-3.4,58.5),(-3.0,59.0),
                           (-1.2,60.4),(1.8,61.5),
                           (5.8,62.0),(13.0,62.0),
                           (13.0,51.0),(0.0,51.0),
                           (-4.9,57.8),(-3.4,58.5)]

    # set up dictionary of boundaries to test crossing
    
    boundary_group = set([])
    
    boundary_group.add(Boundary('Pentland',[(-3.4,58.5),(-3.0,59.0)],NorthSea))
    boundary_group.add(Boundary('OrktoShet',[(-3.0,59.0),(-1.2,60.4)],NorthSea))
    boundary_group.add(Boundary('ShettoNorChan',[(-1.2,60.4),(1.8,61.5)],NorthSea))
    boundary_group.add(Boundary('NorChantoNor',[(1.8,61.5),(5.8,62.0)],NorthSea))
    
    # initialise larvae. 
    nc_file = (track_input_dir + MPA_SOURCE + '.nc')
    nc_fid = Dataset(nc_file, 'r')
    
    larvae_group = set([])

    # use fates variable to get the nlarvae dimension
    fate = nc_fid.variables['fate'][:]  
    nlarvae = len(fate)
#
#      thei gives a netCDF Dimension object, not sure how to extract the 'size'    
#    nlarvae = nc_fid.dimensions['nlarvae']
    
    for i in range(nlarvae):
        lon = nc_fid.variables['longitude'][i,:]
        lat = nc_fid.variables['latitude'][i,:]
        dep = nc_fid.variables['depth'][i,:]
        bed = nc_fid.variables['at bed'][i,:]
        rt = nc_fid.variables['release day'][i]
        fate = nc_fid.variables['fate'][i]
        temp = nc_fid.variables['temperature'][i,:]
# trick to cope with files with no salinity data (salinity not used)
#        sal = nc_fid.variables['temperature'][i,:]
        sal = nc_fid.variables['salinity'][i,:]

        larvae_group.add(Larval_tracks(lon, lat, dep, bed, rt, fate, 
                                       temp, sal, MPA_SOURCE,
                                       RUN_CONST, SWIM_CONST,
                                       T_LOWER, T_UPPER))  
                                           
    nc_fid.close()
    
    group_group_cross(larvae_group, boundary_group)
#    group_group_enter(larvae_group, shapes)
            
    # output the connectivity graph
        
    # build graph
    G.clear()
    
    G.add_node(MPA_SOURCE)
    
    for boundary in boundary_group:
        ncrossed = boundary.get_crossed()
        nstayed = boundary.get_stayed()
        if ncrossed != 0:
            boundary_name = boundary.get_boundary_name()
            ncrossed = ncrossed
            nstayed = nstayed
            d = {}
            d['ncrossed'] = ncrossed
            d['nstayed'] = nstayed
            G.add_edges_from([(MPA_SOURCE,boundary_name,d)])
                       
    # output graph to file
#
    outfile = open(graph_output_dir + MPA_SOURCE + '_crosses.graphml', 'w')
    NX.write_graphml(G,outfile)
    outfile.close()
    

    


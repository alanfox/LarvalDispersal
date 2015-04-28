
"""

@author: af26

24 February 2015 --- Grid, Mpa and Larva classes separated out into 
separate files.

Uses velocity fields from the POLCOMS 1/6 x 1/9 degree model
to advect and diffuse larvae.

Written in Python 2.7 (though I don't think there is anything here which isn't
Python 3) because I couldn't find the basemap package for Python 3 on Windows.

Only tested on Windows using Anaconda python installation.

The data file path will need to be modified. Data file netcdf4.

Includes larval behaviour. Starting at bed, rising through water 
column, then sinking. 

POLCOMS netcdf arrays are indexed [k,j,i] ie [depth, lat, lon]


"""
from netCDF4 import Dataset
import numpy as np
import shapefile
import networkx as NX
from mpa_class import Mpa
from larva_class import Larval_tracks

run_dir = 'polcoms1993/Run_20150422/'

graph_output_dir = run_dir + 'Networkdata/'
track_input_dir = run_dir + 'Trackdata/'
mpa_name_file = open(run_dir + 'MPA_names.txt', 'r') 

# NUM_LARVAE are released at the start of each day for RELEASE_WINDOW days
# for a single release at time zero set RELEASE_WINDOW negative

NUM_LARVAE = 30
RELEASE_WINDOW = 30

STARTDAY = 32

SECONDS_IN_DAY = 60.0 * 60.0 * 24.0
RADIUS_OF_EARTH = 6378160.0
M_TO_DEGREE = 360.0 / (2.0 * np.pi * RADIUS_OF_EARTH)

DT = 3600.0
DTDAYS = DT / SECONDS_IN_DAY

KM = np.array([1.0, 1.0, 0.0002]) #constant diffusion coefficient m2/s

VERTICAL_INTERP = False
ANIMATE = False
SETTLING = True
DEATH = True

# constants for larval behaviour
# larvae released at bed, head upwards with increasing swim speeds up to 
# average age SWIMMAX. Then swimming gradually directed more downwards from 
# age DESCENDAGE up to DEADAGE.

SWIMSLOW = 0.0          #initial swimming speed
SWIMFAST = 0.003      # max swimming speed
SWIMSTART = 0.0         #age in days at which start swimming
SWIMMAX = 14.0          #average age in days at which max swimming speed is reached
DESCENDAGE = 21.0       # average age at which probability of heading down starts
                        # to increase
FULLDESCENDAGE = 42.0    # now fully heading down
MINSETTLEAGE = 30.0     # minimum age at which can settle given suitable 
                        # habitat
DEADAGE = 63.0          # Average age at which dead
# just set DEADAGE to a large value if larvae are not dying
if not DEATH:
    DEADAGE = 1000.0

if DEATH:
    if RELEASE_WINDOW > 0:
        NRUNDAYS = DEADAGE + RELEASE_WINDOW + 1
    else:
        NRUNDAYS = DEADAGE + 1
else:
    NRUNDAYS = 65
    
# bring constants together for passing to larva class

RUN_CONST = [SECONDS_IN_DAY, M_TO_DEGREE, DT, KM, VERTICAL_INTERP]
SWIM_CONST = [SWIMSLOW,SWIMFAST,SWIMSTART,SWIMMAX,DESCENDAGE,FULLDESCENDAGE,
              MINSETTLEAGE,DEADAGE]
    
def read_shapefile(filename):
    sf = shapefile.Reader(filename)
    shapes = sf.shapes()
    records = sf.records()
    return shapes, records

# helper functions
        
def group_settle(mpa_sprite_group, larva_object):
    settled = False
    for mpa in set(mpa_sprite_group):
        if mpa.settles(larva_object):
            settled = True
    return settled

def group_group_settle(larval_sprite_group, mpa_sprite_group):
    for larva in set(larval_sprite_group):
        group_settle(mpa_sprite_group, larva)
              
G = NX.DiGraph()
   
# loop over larval track files
   
for line in mpa_name_file:
    MPA_SOURCE = line[0:-1]
        
    print MPA_SOURCE
    
    # set up group of mpas to test for settling
    
    mpa_group = set([])
    
    # offshore SAC
    shapes, records = read_shapefile('C:/Users/af26/Shapefiles/UK_SAC_MAR_GIS_20130821b/UK_SAC_MAR_GIS_20130821b/SCOTLAND_SAC_OFFSHORE_20121029_SIMPLE3')
    for i in range(len(shapes)):
        mpa_group.add(Mpa(shapes[i], records[i],'OFF_SAC'))
        
    # SAC with marine components
    shapes, records = read_shapefile('C:/Users/af26/Shapefiles/UK_SAC_MAR_GIS_20130821b/UK_SAC_MAR_GIS_20130821b/SCOTLAND_SACs_withMarineComponents_20130821_SIMPLE3')
    for i in range(len(shapes)):
        mpa_group.add(Mpa(shapes[i], records[i],'MAR_SAC'))
        
    # Nature conservation MPA
    shapes, records = read_shapefile('C:/Users/af26/Shapefiles/MPA_SCOTLAND_ESRI/MPA_SCOTLAND_SIMPLE3')
    for i in range(len(shapes)):
        mpa_group.add(Mpa(shapes[i], records[i],'MPA'))
        
    # Irish SACs
    shapes, records = read_shapefile('C:/Users/af26/Shapefiles/SAC_ITM_WGS84_2015_01/SAC_Offshore_WGS84_2015_01')
    for i in range(len(shapes)):
        mpa_group.add(Mpa(shapes[i], records[i],'IRISH'))
        
    # initialise larvae. 
    nc_file = (track_input_dir + MPA_SOURCE + '.nc')
    nc_fid = Dataset(nc_file, 'r')
    
    # Using grids of larvae at the same depth around a central point.
    
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

        larvae_group.add(Larval_tracks(lon,lat,dep,bed,rt,fate, MPA_SOURCE,
                                           RUN_CONST, SWIM_CONST))  
                                           
    nc_fid.close()

    group_group_settle(larvae_group, mpa_group)
            
    # output the connectivity graph
    
    
    # build graph
    G.clear()
    
    G.add_node(MPA_SOURCE)
    
    for mpa in mpa_group:
        nsettled = mpa.get_settled()
        if nsettled != 0:
            mpa_name = mpa.get_sitename()
            weight = float(nsettled) / float(nlarvae)
            G.add_weighted_edges_from([(MPA_SOURCE,mpa_name,weight)])
                       
    # output graph to file
#
    outfile = open(graph_output_dir + MPA_SOURCE + '.graphml', 'w')
    NX.write_graphml(G,outfile)
    outfile.close()
    

    



"""

Takes netcdf files output by larval_dispersal_polcoms.py and builds 
a network of connections between MPAs.

Mostly testing whether larvae have entered another MPA while at the 
bed and in a 'settling' phase.

Also checks that larvae remain within the temperature range in which
they are viable.


"""
from netCDF4 import Dataset
import numpy as np
import shapefile
import networkx as NX
from mpa_class import Mpa
from larva_class import Larval_tracks
import platform

if platform.system() == 'Windows':
    run_dir = ('C:/Users/af26/LarvalDispersalResults/'
            + 'polcoms1990/Run_1000_baseline/')
elif platform.system() == 'Linux':
    run_dir = ('/home/af26/LarvalModelResults/Polcoms1990/Run_test/')

graph_output_dir = run_dir + 'Networkdata/'
track_input_dir = run_dir + 'Trackdata/'
mpa_name_file = open(run_dir + 'MPA_names.txt', 'r') 

input_data_file = open(run_dir + 'input.dat', 'r')

input_dict = {}

for line in input_data_file:
    wordlist = line.split()
    input_dict[wordlist[0]] = wordlist[-1]
        
mpa_name_file = open(run_dir + 'MPA_names.txt', 'r') 

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
              MINSETTLEAGE,DEADAGE,TARGETDEPTH]
    
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
    if platform.system() == 'Windows':    
        MPA_SOURCE = line[0:-1]
        shapefile_root =   'C:/Users/af26/Shapefiles/' #windows
    elif platform.system() == 'Linux':
        MPA_SOURCE = line[0:-2]
        shapefile_root =   '/home/af26/Shapefiles/' #linux
    
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

    # Mikael Dahl's lophelia sites
    shapes, records = read_shapefile('C:/Users/af26/Shapefiles/MikaelDahl/MikaelDahl_1')
    for i in range(len(shapes)):
        mpa_group.add(Mpa(shapes[i], records[i],'Dahl'))
        
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
        temp = nc_fid.variables['temperature'][i,:]
# trick to cope with files with no salinity data (salinity not used)
        sal = nc_fid.variables['temperature'][i,:]
#        sal = nc_fid.variables['salinity'][i,:]

        larvae_group.add(Larval_tracks(lon, lat, dep, bed, rt, fate, 
                                       temp, sal, MPA_SOURCE,
                                       RUN_CONST, SWIM_CONST,
                                       T_LOWER, T_UPPER))  
                                           
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
    

    


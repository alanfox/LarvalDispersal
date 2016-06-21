
"""
Created on Wed 24 June 2015

@author: af26

Reads track data from netCDF4 format files and produces an animation

"""
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
import shapefile
#from matplotlib.colors import ListedColormap, BoundaryNorm
from mpa_class import Mpa
from grid_class import Grid
import platform
import datetime as datetime
import matplotlib as matplotlib

matplotlib.rcParams['font.sans-serif'] = 'Arial'

# for reading the MPA shapefiles. For drawing the mpa outlines.    

def read_shapefile(filename):
    sf = shapefile.Reader(filename)
    shapes = sf.shapes()
    records = sf.records()
    return shapes, records
            
    
# set up group of mpas
def read_mpas():
    mpa_group = set([])
    
    if platform.system() == 'Windows':    
        shapefile_root =   'C:/Users/af26/Shapefiles/' #windows
    elif platform.system() == 'Linux':
        shapefile_root =   '/home/af26/Shapefiles/' #linux
    
    shapes, records = read_shapefile(shapefile_root + 
                    'UK_SAC_MAR_GIS_20130821b/UK_SAC_MAR_GIS_20130821b/' + 
                    'SCOTLAND_SAC_OFFSHORE_20121029_SIMPLE3')
    for i in range(len(shapes)):
        mpa_group.add(Mpa(shapes[i], records[i],'OFF_SAC'))
        
    # SAC with marine components
    shapes, records = read_shapefile(shapefile_root + 
                    'UK_SAC_MAR_GIS_20130821b/UK_SAC_MAR_GIS_20130821b/' + 
                    'SCOTLAND_SACs_withMarineComponents_20130821_SIMPLE3')
    for i in range(len(shapes)):
        mpa_group.add(Mpa(shapes[i], records[i],'MAR_SAC'))
        
    # Nature conservation MPA
    shapes, records = read_shapefile(shapefile_root + 
                    'MPA_SCOTLAND_ESRI/MPA_SCOTLAND_SIMPLE3')
    for i in range(len(shapes)):
        mpa_group.add(Mpa(shapes[i], records[i],'MPA'))
        
    # Irish SACs
    shapes, records = read_shapefile(shapefile_root + 
                    'SAC_ITM_WGS84_2015_01a/SAC_Offshore_WGS84_2015_01')
    for i in range(len(shapes)):
        mpa_group.add(Mpa(shapes[i], records[i],'IRISH'))
        
    # Mikael Dahl's lophelia sites
    shapes, records = read_shapefile(shapefile_root + 
                    'MikaelDahl/MikaelDahl_1')
    for i in range(len(shapes)):
        mpa_group.add(Mpa(shapes[i], records[i],'Dahl'))
        
    return mpa_group
    

# INITIALISE STUFF
# need t grids for the 2d histogram

#nc_filedzu = 'C:/Users/af26/PolcommModelData/depth_dz_S12run420_UV.nc'
nc_filedzt = 'C:/Users/af26/PolcommModelData/depth_dz_S12run420_TS.nc'
#nc_fiddzu = Dataset(nc_filedzu, 'r')
nc_fiddzt = Dataset(nc_filedzt, 'r')

    # read in and calculate the model grid variables
    
gridt = Grid(nc_fiddzt)
#gridu = Grid(nc_fiddzu)

nc_fiddzt.close()
#nc_fiddzu.close()

latitude = gridt.get_latitude()
longitude = gridt.get_longitude()

# read in GEBCO bathymetry

nc_infile = ('C:/Users/af26/GEBCO/GEBCO_2014_2D_-20.0_40.0_13.0_65.0.nc')

nc_gebco = Dataset(nc_infile,'r')

bathymetry = nc_gebco.variables['elevation'][:]
lat_bath = nc_gebco.variables['lat'][:]    
lon_bath = nc_gebco.variables['lon'][:]    

nc_gebco.close()

# lophelia mpa subset for plotting
MPA_LOPHELIA = (['Anton Dohrn Seamount',
               'Darwin Mounds',
               'East Mingulay',
               'East Rockall Bank',
               'Faroe-Shetland Sponge Belt',
               'Hatton Bank',
               'North West Rockall Bank',
               'Rosemary Bank Seamount',
               'South-East Rockall Bank SAC',
               'Wyville Thomson Ridge',
               'Norwegian Boundary Sediment Plain',
               'Central Fladen',
               'South-West Porcupine Bank SAC',
               'Belgica Mound Province SAC','Hovland Mound Province SAC',
               'North-West Porcupine Bank SAC','Porcupine Bank Canyon SAC'])

MPA_SCOT =    (['Anton Dohrn Seamount',
               'Darwin Mounds',
               'East Mingulay',
               'East Rockall Bank',
               'Faroe-Shetland Sponge Belt',
               'Hatton Bank',
               'North West Rockall Bank',
               'Rosemary Bank Seamount',
               'South-East Rockall Bank SAC',
               'Wyville Thomson Ridge',
               'Norwegian Boundary Sediment Plain',
               'Central Fladen'])

MPA_RB =      (['Anton Dohrn Seamount',
               'East Rockall Bank',
               'Hatton Bank',
               'North West Rockall Bank',
               'South-East Rockall Bank SAC'])

MPA_NW =      (['Darwin Mounds',
               'Faroe-Shetland Sponge Belt',
               'Rosemary Bank Seamount',
               'Wyville Thomson Ridge'])

MPA_EM =      (['East Mingulay'])

MPA_NS =      (['Norwegian Boundary Sediment Plain',
               'Central Fladen'])
MPA_PORC =    (['South-West Porcupine Bank SAC',
               'Belgica Mound Province SAC','Hovland Mound Province SAC',
               'North-West Porcupine Bank SAC','Porcupine Bank Canyon SAC'])
MPA_SOUTH =   (['Aviles Canyon',
                'Galicia Mounds'
                ])

#MPA_GROUPS = ([MPA_SCOT,MPA_RB,MPA_NW,MPA_EM])
MPA_GROUPS = ([MPA_SOUTH,MPA_SOUTH,MPA_PORC,MPA_PORC])
#MPA_GROUPS = ([MPA_PORC,MPA_PORC,MPA_PORC])

RUN_NAME = ['baseline','behaviour2','baseline','behaviour2']
#RUN_NAME = ['baseline','behaviour2','doublelife']

def add_annual_data_to_heat_map(MPA_SOURCE, run_dir, run_name):
    x = []
    y = []
    nc_fid = {}
    for mpa in MPA_SOURCE:
        filename = run_dir + 'Trackdata/' + mpa + '.nc'
        nc_fid = Dataset(filename, 'r')
        
    
    # dead or alive
        temperature = nc_fid.variables['temperature'][:,:]
        not_dead = temperature > min_temp
        not_dead = np.cumprod(not_dead,axis = 1)
            
        releasedays = nc_fid.variables['release day'][:]
        
# loop over all larvae
        ii = 0
        for releaseday in releasedays:
            lat = nc_fid.variables['latitude'][ii,:]
            lon = nc_fid.variables['longitude'][ii,:]
            temperature = nc_fid.variables['temperature'][ii,:]
            if run_name == 'doublelife':
                days = range(70,126)
            else:
                days = range(35,63)

            for day in days:
                i = day * daysteps
                if (not_dead[ii,i]):
                    x.append(lon[i])
                    y.append(lat[i])
            ii += 1
        nc_fid.close()
        
    return x,y

def plot_heatmap(x,y,MPA_SOURCE,axarr):
    
    for mpa in mpa_group:
#        print mpa.get_sitename()
        if mpa.get_sitename() in MPA_SOURCE:
            colour = 'limegreen'
            mpa.plot_shape3(axarr,colour)                        
        elif mpa.get_sitename() in MPA_LOPHELIA:
            colour = 'red'
            mpa.plot_shape3(axarr,colour)  
                                  
    counts, xedges, yedges, Image = axarr.hist2d(x,y,[longitude,latitude],
                                                 cmin = 0.0001,
                                                 vmin = 0.0, vmax = 50000)
    axarr.tick_params(length = 1.0)
    for axis in ['top','bottom','left','right']:
        axarr.spines[axis].set_linewidth(0.5)
    axarr.set_xlim([-20, 13])
    axarr.set_ylim([40, 55])
#    axarr.set_ylim([50, 61])
    for tl in axarr.get_yticklabels():
        tl.set_size(6)
    for tl in axarr.get_xticklabels():
        tl.set_size(6)
# save 300 m contour line
    cs = axarr.contour(lon_bath,lat_bath,bathymetry,[-1000,-200],
                     colors = 'grey',linestyles = 'solid', linewidths = 0.5)
#    axarr.clabel(cs,fontsize = 4, inline_spacing = -2, fmt='%4.0f')                 
    # colours land in black
    axarr.contourf(lon_bath,lat_bath,bathymetry,[0,5000], colors = 'lightgreen')
    axarr.contour(lon_bath,lat_bath,bathymetry,[0,5000], 
                  colors = 'black', linewidths = 0.5)

    return Image
    
daysteps = 24.0
    
# minimum survivable temperature
    
min_temp = 3.0
# read the mpa shapefiles
  
mpa_group = read_mpas()    

#4 plots
#fig, axarr = plt.subplots(2,2,figsize = (7,4.96),sharex='col', sharey='row')
fig, axarr = plt.subplots(2,2,figsize = (4.48,3.2),sharex='col', sharey='row')
#3 plots
#fig, axarr = plt.subplots(3,figsize = (3.5,7.0),sharex='col')

icount = 0
                
for MPA_SOURCE in MPA_GROUPS:
    x = []
    y = []
    z = []
    
         
    for year in range(1965,2005):
        year_str = str(year)
        print year
        
        # timesteps per day
        
        #startday = 32.0
        
        starttime = datetime.datetime(year,2,1)
        
        if platform.system() == 'Windows':
            run_dir = ('E:/af26/LarvalDispersalResults/'
                    + 'polcoms'+year_str+'/Run_1000_' + RUN_NAME[icount] + '/')
        elif platform.system() == 'Linux':
            run_dir = ('/home/af26/LarvalModelResults/Polcoms1990/Run_test/')
            
        # open the larval track datasets. Put them in a dictionary.
        x1,y1 = add_annual_data_to_heat_map(MPA_SOURCE, run_dir, RUN_NAME[icount])    
        
        x = x + x1
        y = y + y1
#4 plots
    Image = plot_heatmap(x,y,MPA_SOURCE,axarr[icount//2,icount%2])
#3 plots
#    Image = plot_heatmap(x,y,MPA_SOURCE,axarr[icount])
    icount = icount + 1
    
#4 plots
axarr[0,0].set_ylabel('latitude', fontsize = 6)
axarr[1,0].set_ylabel('latitude', fontsize = 6)
axarr[1,0].set_xlabel('longitude', fontsize = 6)
axarr[1,1].set_xlabel('longitude', fontsize = 6)
textstr = 'A'
axarr[0,0].text(-19, 64,'E', fontsize = 8)
axarr[0,1].text(-19, 64,'F', fontsize = 8)
axarr[1,0].text(-19, 64,'G', fontsize = 8)
axarr[1,1].text(-19, 64,'H', fontsize = 8)
# 3 plots
#axarr[0].set_ylabel('latitude', fontsize = 6)
#axarr[1].set_ylabel('latitude', fontsize = 6)
#axarr[2].set_ylabel('latitude', fontsize = 6)
#axarr[2].set_xlabel('longitude', fontsize = 6)
#axarr[0].text(-19, 60,'A', fontsize = 8)
#axarr[1].text(-19, 60,'B', fontsize = 8)
#axarr[2].text(-19, 60,'C', fontsize = 8)

# Now adding the colorbar
#cbaxes = fig.add_axes([0.05, 0.03, 0.92, 0.03]) 
cbaxes = fig.add_axes([0.07, 0.04, 0.925, 0.03]) 
#cb = plt.colorbar(Image, ax = axarr[0,0], orientation='horizontal', cax = cbaxes)
cb = plt.colorbar(Image, ax = axarr[0], orientation='horizontal', 
                  cax = cbaxes)
cb.ax.tick_params(labelsize=6)        
plt.setp(cb.ax.get_xticklabels()[-1], visible=False)

cb.outline.set_linewidth(0.5)

#plt.tight_layout(pad=1)
# 4 plots
#plt.subplots_adjust(left=0.05, bottom=0.13, right=0.97, top=0.98, 
#                    wspace=0.025, hspace=0.025)        
plt.subplots_adjust(left=0.07, bottom=0.17, right=0.995, top=0.995, 
                    wspace=0.025, hspace=0.025)        
# 3 plots
#plt.subplots_adjust(left=0.1, bottom=0.11, right=0.97, top=0.98, 
#                    wspace=0.025, hspace=0.025)        

plt.show()

#filenameroot = ('C:/Users/af26/Documents/Papers/LopheliaConnectivity/' + 
#                 'heat_maps_doublelife')
#
#plt.savefig(filenameroot + '.pdf')
#plt.savefig(filenameroot + '.eps', format='eps', dpi=1200)
#plt.savefig(filenameroot + '.png', dpi = 1000)
#            

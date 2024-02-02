""" FireConsts
This is the module containing all constants used in this project
Also used for run control: which region is being used?
"""

# region
region = ''
stralpha = 'alpha1000'   # string of alpha parameter

# project directories
dirhome = ''               # home directory
dirpjcode = dirhome + ''        # project code directory
dirdata = dirhome + ''         # input data directory
dirextdata = ''
dirNOAA20 = dirextdata + ''
lakedir = dirextdata + ''
lcdir = dirextdata + ''

# set (output) project data directory
dirpjdata = dirhome + ''
if region:
    dirpjdata = dirpjdata + region + '/'
if stralpha:
    dirpjdata = dirpjdata + stralpha + '/'
    
dir_temp = dirpjdata + 'temp/'          # directory for temporary file


# shapefiles for regional clipping
shp_all = {'AK':dirdata+'shapefiles/Alaska.gpkg', 'NT':dirdata+'shapefiles/nwt.gpkg',
           'CA':dirdata+'shapefiles/canada.gpkg', 'YAK':dirdata+'shapefiles/yakutia.gpkg',
           'RU':dirdata+'shapefiles/russia.gpkg', 'full':dirdata+'shapefiles/wwf_borealarctic_diss_05deg.gpkg'}
if region:
    if not region == 'custom':
        shp_fnm = shp_all[region]
    else:
        shp_fnm = None
        regext = [-118,60.4,-114,62.4]
else:
    shp_fnm = shp_all['full']

# which peatmap to use
peatmap = 'hugelius' # 'peatml' or 'hugelius' or 'yu_hug' (default is a merged version)

# which satellite to use
sat = 'SNPP' # 'SNPP' | 'NOAA20' | 'VIIRS'

# parameters used for fire pixel clustering
EARTH_RADIUS_KM = 6371.0  # earth radius, km
SPATIAL_THRESHOLD_KM = 1  # threshold of spatial distance (between pixels, km) used to do fire initial clustering

# fire type and visualization parameters
lc_dict = { # dictionary for simplifying ESA CCI land cover types
           0:0, 20:0, 21:0, 22:0, # no data, water, ice, bare
           1:1, 2:1, 3:1, 4:1, # cropland
           5:2, 6:3, 7:4, 8:5, 9:6, # forest
           10:7, 11:7, 12:7, # shrubs & mosaic landscape
           13:8, 14:9, 15:10, # grassland, lichen, sparse vegetation
           16:11, 17:11, 18:12, # flooded forests and shrubs
           19:13} # urban
FTYP = {0:'Other', 1:'Crop', 2:'forest-peat', 3:'forest-no peat', 4:'shrubs/mosaic/open-peat', 
        5:'shrubs/mosaic/open-no peat', 6:'grassland:peat', 7:'grassland-no peat', 8:'urban'} 

# temporal and spatial distances for fire object definition
maxoffdays = {0:4, 1:5, 2:10, 3:10, 4:10, 5:10, 6:10, 7:3, 8:5}   # fire becomes inactive after this number of consecutive days without active fire detection
sleepdays = 30 # this is the additional threshold that involves checking for active fireline
CONNECTIVITY_THRESHOLD_KM = {0:1, 1:1, 2:2.5, 3:2.5, 4:2.5, 5:2.5, 6:1, 7:1, 8:1} # km, corresponding to fire type (forest/shrub: 5km; other: 1km)
sleeperthresh = 1 #km, distance a sleeper is allowes to spread since its last active fire detection

# alpha parameter
valpha = 1000      # 1000 m (default)
lalpha = 100


# lengths and areas related to VIIRS pixl data
VIIRSbuf = 187.5   # fire perimeter buffer (deg), corresponding to 375m/2 at lat=30
fpbuffer = 100     # buffer use to determine fire line pixels (deg), is added to VIIRS buffer!
flbuffer = 100     # buffer for determination of fire line length, ~500m
extbuffer = 1000   # buffer to define interior/exterior region, 1000 m
area_VI = 0.141    # km2, area of each 375m VIIRS pixel

# mcd bffer (not used)
MCD64buf = 231.7   # fire perimeter buffer (deg), corresponding to 463.31271653 m/2

# parameter used for determine large fires
falim = 0.1    # large fire area threshold (in km2)

# diagnostic data name and types saved in geojson files (in addition to geometries)
dd = {'fireid':'int',                  # id
      'mergid':'int',               # this is the id in the large fire database
      'ftype':'str',                # fire type
      'isactive':'int',             # active status
      't_inactive':'int',           # how long has it been inactive
      'isignition':'int',           # is this a new ignition?
      'invalid':'int',              # invalid status
      'n_pixels':'int',             # number of total pixels
      'n_newpixels':'int',          # number of new pixels
      'farea':'float',              # fire size
      'fperim':'float',             # fire perimeter length
      'flinelen':'float',           # active fire front line length
      'duration':'float',           # fire duration
      'pixden':'float',             # fire pixel density
      'meanFRP':'float',            # mean FRP of the new fire pixels
      'FRP95':'float',              # 95th percentile of FRP of the new fire pixels
      'tst_year':'int',             # t_st[0]
      'tst_month':'int',
      'tst_day':'int',
      'tst_ampm':'str',
      'ted_year':'int',             # t_ed[0]
      'ted_month':'int',
      'ted_day':'int',
      'ted_ampm':'str',
      'ted_doy':'int',
      'lake_area':'float',          # area within fire perimeter covered by lakes
      'lake_no':'int',
      'lake_border':'float',        # fire perimeter length that ends at lake
      'lake_border_tot':'float'
      }

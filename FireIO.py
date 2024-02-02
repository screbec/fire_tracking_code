""" FireIO
This module include functions used to read and save data
"""
# ------------------------------------------------------------------------------
#%% Read active fire data
# ------------------------------------------------------------------------------



def viirs_pixel_size(sample,band='i',rtSCAN_ANGLE=False):
    '''calculate approximate size of i-band (375-m) or m-band (750-m) VIIRS pixel
        Adapted from L. Giolio's IDL code
    Usage: DS, DT = viirs_pixel_size(200,band='m',rtSCAN_ANGLE=False)

    Parameters
    ----------
    sample : int
        sample number
    band : str, 'i'|'m'
        i (default) or m band
    rtSCAN_ANGLE : bool, default = False
        flag to return scan angle

    Returns
    -------
    DT : float
        length in along-track dimension [km]
    DS : float
        length in along-scan dimension [km]
    '''
    import numpy as np

    # set constants
    earth_radius = 6371.   # earth radius [km]
    h = 824.               # SUOMI-NPP orbit altitude [km]
    pt = 0.361             # nadir pixel resolution [km]
    scan_rate = 56.06/6304 # scan rate (deg)
    st1,st2,st3 = 1184, 1920, 3200  # sample tier steps
    sb1,sb2,sb3 = 0, 3552, 5024  # sample base for each tier
    ps_const = 0.128  # constant to convert zone number to along-scan deg for 1 pixel

    # adjust constants for m band
    if band == 'm':
        pt *= 2
        scan_rate *= 2
        st1 /= 2
        st2 /= 2
        st3 /= 2
        sb1 /= 2
        sb2 /= 2
        sb3 /= 2
        ps_const *= 2

    # derive more constants
    st = pt/h              # along-track deg for 1 pixel[rad]
    r = earth_radius + h   # r for satellite[km]

    # calculate along-scan degrees
    abs_sample = (sample>=st3)*(sample+1-st3) + (sample<st3)*(st3-sample)
    zone1 = abs_sample <= st1
    zone2 = np.logical_and(abs_sample > st1,abs_sample <= st2)
    zone3 = abs_sample > st2
    ps = ps_const*(3*zone1 + 2*zone2 + zone3)
    ss = ps/h              # along-scan deg for 1 pixel[rad]
    abs_scan = zone1*(abs_sample*3 * scan_rate) + \
               zone2*((sb2+((abs_sample-st1)*2))*scan_rate) + \
               zone3*((sb3+(abs_sample-st2))*scan_rate)      # scan angle, in degree
    theta = np.deg2rad(abs_scan)  # scan angle, in radian
    cos_theta = np.cos(theta)

    # calculate pixel size DS and DT
    temp = (earth_radius/r)**2.-np.sin(theta)**2.
    sqrt_temp = np.sqrt(temp)
    DS = earth_radius * ss * (cos_theta/sqrt_temp-1.)
    DT = r * st * (cos_theta-sqrt_temp)

    if rtSCAN_ANGLE == True:
        scan_angle = np.rad2deg(theta)
        return DS, DT, scan_angle
    else:
        return DS, DT

def save_af(data,d,sensor='SNPP'):
    ''' Save a daily allfires object to a pickle file

    Parameters
    ----------
    data : pd with all fires of a month
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    sensor: Suomi NPP ('SNPP') or NOAA-20 ('NOAA20')
    '''

    import pickle
    import os
    from FireConsts import dirpjdata
    
    # create head from sensor
    if sensor == 'NOAA20':
        head = 'VJ14IMGML.'
    else:
        head = 'VNP14IMGML.'
    
    # get output file name
    dir_temp = os.path.join(dirpjdata,'temp',str(d.year)) + '/'
    fnm = os.path.join(dir_temp,head+d.strftime('%Y%m')+'.pkl')

    # check folder
    check_filefolder(fnm)

    # save
    with open(fnm,'wb') as f:
        pickle.dump(data, f)

def load_af_pkl(d,sensor='SNPP'):
    ''' Function used to load regional filtered data to a temporary directory
    '''
    from FireConsts import dirpjdata
    import os
    import pickle
    
    # create head from sensor
    if sensor == 'NOAA20':
        head = 'VJ14IMGML.'
    else:
        head = 'VNP14IMGML.'

    # temporary file name
    dir_temp = os.path.join(dirpjdata,'temp',str(d.year)) + '/'
    fnm_tmp = os.path.join(dir_temp,head+d.strftime('%Y%m')+'.pkl')

    # load and return
    if os.path.exists(fnm_tmp):
        with open(fnm_tmp,'rb') as f:
            df = pickle.load(f)
            return df
    else:  # if no presaved file, return None
        return None

def read_VNP14IMGML(d):
    ''' read monthly VNP14IMGML data

    Parameters
    ----------
    d : date from which year and month will be axtracted

    Returns
    -------
    df : dataframe
        the monthly dataframe of VIIRS active fires
    '''
    from FireConsts import dirextdata
    import os
    import pandas as pd

    # set monthly file name
    dirFC = os.path.join(dirextdata,'VIIRS','VNP14IMGML05') + '/'
    fnmFC = os.path.join(dirFC,'VNP14IMGML.'+d.strftime('%Y%m')+'.C1.05.txt')

    # read
    if os.path.exists(fnmFC):
        try:
            usecols = ['YYYYMMDD','HHMM','Lat','Lon','Line','Sample','FRP','Confidence','Type','DNFlag']
            df = pd.read_csv(fnmFC,parse_dates=[['YYYYMMDD','HHMM']],usecols=usecols,skipinitialspace=True)
        except ValueError:
            usecols = ['YYYYMMDD','HHMM','Lat','Lon','DT','DS','FRP','Confidence','Type','DNFlag']
            df = pd.read_csv(fnmFC,parse_dates=[['YYYYMMDD','HHMM']],usecols=usecols,skipinitialspace=True)
            #### TIME PARSING DOES NOT WORK WITH NEW FILES!
        
        # in Collection 04, some FRP values are '*******', causing problems
        if df.dtypes.FRP.name == 'object':
            df.FRP = df.FRP.replace('*******',0).astype('float')
        
        # filter for type and quality (use types 0 (vf) and 3 (offshore))
        df = df.loc[((df['Type'] == 0) | (df['Type'] == 3)) &
                    ((df['Confidence'] == 'nominal') | (df['Confidence'] == 'high'))]
        
        # compute pixel dimensions form sample and line
        if 'Sample' in df:
            df['DT'],df['DS'] = viirs_pixel_size(df['Sample'].values)
            df = df.drop(columns=['Sample','Line'])
        
        # compute normalised frp
        df.FRP = df.FRP/(df.DT*df.DS)
        
        # add the satellite information
        df = df.assign(Sat = 'SNPP')
        return df
    else:
        return None

def read_VJ114IMGML(yr,mo):
    ''' read monthly VNP14IMGML data

    Parameters
    ----------
    yr : int
        year
    mo : int
        month

    Returns
    -------
    df : dataframe
        the monthly dataframe of VIIRS active fires
    '''
    from FireConsts import dirNOAA20
    import os
    import pandas as pd

    # set monthly file name
    fnmFC = os.path.join(dirNOAA20,str(yr),'VJ114IMGML_'+str(yr)+str(mo).zfill(2)+'.txt')

    # read
    usecols = ['year','month','day','hh','mm','lon','lat','mask',
               'line','sample','frp']
    def parser(yr,mo,dy, h, m):
        return pd.to_datetime(yr + '-' + mo + '-' + dy + ' ' + h + ':' + m,
                              format='%Y-%m-%d %H:%M')
    if os.path.exists(fnmFC):
        # df = pd.read_csv(fnmFC,usecols=usecols,skipinitialspace=True)
        df = pd.read_csv(fnmFC, parse_dates={'YYYYMMDD_HHMM': ['year','month','day','hh','mm']},
                         date_parser=parser,usecols=usecols,skipinitialspace=True)
        # harmonise column names
        df = df.rename(columns={'lat':'Lat','lon':'Lon','frp':'FRP',
                                'line':'Line','sample':'Sample'})
        # compute pixel size from sample
        df['DT'],df['DS'] = viirs_pixel_size(df['Sample'].values)
        # filter for confidence
        df = df.loc[df['mask'] > 7] # 7:low conf., 8:nominal, 9:high
        # add the satellite information
        df = df.assign(Sat = 'NOAA20')
        return df
    else:
        return None

def AFP_setDN(gdf):
    ''' Function to set daynight column (using local hour) and update the df
    also removes geometry column

    Parameters
    ----------
    gdf : pandas DataFrame
        fire pixel data, with 'Lon' and 'YYYYMMDD_HHMM' column

    Returns
    -------
    df : pandas DataFrame
        the DataFrame with 'DN' column
    '''
    import pandas as pd
    import numpy as np
    
    df = pd.DataFrame(gdf.drop(columns='geometry'))
    
    # convert utc time to local solar time (rough estimate, can be 15 mins off)
    df['diff_lst'] = pd.to_timedelta(df.Lon /15, unit = 'hours')
    df['datetime_lst'] = df['YYYYMMDD_HHMM']+df['diff_lst']
    
    # extract separate date and hour columns
    df['date_lst'] = df.datetime_lst.dt.date
    df['hour_lst'] = df.datetime_lst.dt.hour
            
    # set am/pm flag based on local hour
    df = df.assign(DN = np.where(((df.hour_lst > 6) & (df.hour_lst < 18)), 'day','night'))

    return df

def AFP_regfilter(df,strlat='Lat',strlon='Lon'):
    ''' filter fire pixels using a given shp_Reg

    Parameters
    ----------
    df : pandas DataFrame
        fire pixel data with lat and lon information
    shp_Reg : geometry
        the geometry of the region shape
    strlat : str
        the column name for the latitude
    strlon : str
        the column name for the longitude

    Returns
    -------
    df_filtered : pandas DataFrame
        the filtered fire pixels
    '''
    from shapely.geometry import Point
    import geopandas as gpd
    from FireConsts import shp_fnm
    
    # set extent from region shapefile
    if shp_fnm:
        shp = get_any_shp(shp_fnm)
        regext = shp.bounds
    else:
        from FireConsts import regext
    
    # preliminary spatial filter and quality filter
    newfirepixels = df.loc[(df[strlat] >= regext[1]) & (df[strlat] <= regext[3]) &
                          (df[strlon] >= regext[0]) & (df[strlon] <= regext[2])]
    point_data = [Point(xy) for xy in zip(newfirepixels[strlon], newfirepixels[strlat])]
    df_filtered = gpd.GeoDataFrame(newfirepixels, geometry=point_data, crs=4326)

    # if shp is not a rectangle, do detailed filtering (within shp)
    if shp_fnm:
        import math
        if not math.isclose(shp.minimum_rotated_rectangle.area, shp.area):
            df_filtered = df_filtered[df_filtered['geometry'].within(shp)]

    return df_filtered

def AFP_spurfilter(df, t):
    '''filter for known spurious fire pixels'''
    # there is spurious data on June 24, 2020, 11:06
    if (t[0] == 2020) and (t[1] == 6):
        df = df.loc[~(df['YYYYMMDD_HHMM'] == '2020-06-24 11:06:00')]
    return df

def AFP_gasflarefilter(gdf):
    """ Function used to filter out known areas of gas flaring
    Parameters
    ----------
    gdf : geopandas GeoDataFrame
        fire pixel data
    Returns
    -------
    gdf_filtered : geopandas GeoDataFrame
        the filtered fire pixels
    """

    # filter non-veg fires using a pre-derived mask
    gdf_gf = read_gasflare()
    buf_gf = 2000  # buffer for each gas flare (now set to 2000m)
    gdf_filtered = gdf[(gdf.geometry.within(gdf_gf.iloc[0].geometry.buffer(buf_gf)) == False)]

    return gdf_filtered

def AFP_toprj(gdf):
    '''Transforms lat/lon coordinates to projected coordinate system for computations in m
    for global studies probably proj:cea would be a good choice
    for the boreals we use North Pole LAEA (epsg:3571)'''
    
    gdf = gdf.to_crs(3571)
    gdf['x'] = gdf.geometry.x  
    gdf['y'] = gdf.geometry.y
    
    return gdf

def read_AFPVIIRS(t, sat='SNPP',cols=['Lat','Lon','FRP','Sat','DT','DS','YYYYMMDD_HHMM','DN','x','y']):
    ''' Read half-daily fire location for a region from monthly NOAA20 VIIRS fire data.
        For the monthly VIIRS data, in order to reduce repeat calculation, we
        do the spatial and quality filtering once and save the filtered data to
        files in a temporary directory. For most time steps, we read the pre-saved
        VIIRS data.

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the intialization
    region : 2-element tuple [regnm, regshp].
        Regnm is the name of the region; Regshp can be
         - a geometry
         - a four-element list showing the extent of the region [lonmin,latmin,lonmax,latmax]
         - a country name
    cols : list
        the column names of the output dataframe
    Returns
    -------
    df_AFP : pandas DataFrame
        list of active fire pixels (filtered) from NOAA20
    '''

    from datetime import date

    # read from pre-saved file
    d = date(*t[:-1])

    # load from pre-saved file
    df = load_af_pkl(d,sensor=sat)

    # if no pre-saved file, read from original VNP14IMGML file and do/save spatial filtering
    if df is None:
        
        # read VJ114IMGML monthly file
        if sat=='SNPP':
            df = read_VNP14IMGML(d)
        elif sat=='NOAA20':
            df = read_VJ114IMGML(t[0],t[1])
        else:
            print('please set SNPP or NOAA20 for sat')
        
        if df is not None: # NOAA data can be empty
            # filter for known spurious pixels
            df = AFP_spurfilter(df, t)
            
            # do regional filtering
            gdf = AFP_regfilter(df)
            
            # transform all coordinates to projected
            gdf = AFP_toprj(gdf)
            
            # filter our potential gas flares
            gdf = AFP_gasflarefilter(gdf)
            
            # set ampm
            df = AFP_setDN(gdf)
            
            # save to temporary file
            save_af(df,d,sensor=sat)
    
    if df is not None:
        # extract active pixels at current time step  (day and ampm filter)
        df = df.loc[(df['date_lst']==d)]
        if t[-1] == 'PM':   # pm corresponds to daytime (~1.30pm) detections
            df = df.loc[df['DN'] == 'day']
        else:
            df = df.loc[df['DN'] == 'night']
    
        # return selected columns; need to change column names first
        df = df[cols]

    return df

def read_AFP(t,src='SNPP'):
    ''' The wrapper function used to read and extract half-daily active fire
    pixels in a region.

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' during the intialization
    src : str
        'SNPP' | 'NOAA20' | 'VIIRS'
    region : 2-element tuple [regnm, regshp].
        Regnm is the name of the region; Regshp can be
         - a geometry
         - a four-element list showing the extent of the region [lonmin,latmin,lonmax,latmax]
         - a country name
    Returns
    -------
    vlist : pandas DataFrame
        extracted active fire pixel data used for fire tracking
    '''
    
    if src == 'VIIRS':
        import pandas as pd
        vlist_SNPP = read_AFPVIIRS(t,sat='SNPP')
        if t[0] > 2017: # NOAA data is only availbale form 2018 onwards
            vlist_NOAA20 = read_AFPVIIRS(t,sat='NOAA20')
        else:
            vlist_NOAA20 = None
        if vlist_NOAA20 is None: # NOAA20 data is missing for some years and months
            vlist = vlist_SNPP
        else:
            vlist = pd.concat([vlist_NOAA20,vlist_SNPP],ignore_index=True)
    elif src == 'SNPP':
        vlist = read_AFPVIIRS(t,sat=src)
    elif src == 'NOAA20':
        vlist = read_AFPVIIRS(t,sat=src)
    else:
        print('Please set src to SNPP, NOAA20, or VIIRS')
        return None
    
    # if there is not fire data for this month available we return an empty list
    if vlist is None: # can be empty when only using NOAA20
        vlist = [] # (prevents the len(afp) check in FireForward rom throwing error)

    return vlist.reset_index(drop=True)

# ------------------------------------------------------------------------------
#%% Read other datasets
# ------------------------------------------------------------------------------
def read_gasflare():
    """ read NOAA global gas flaring dataset
    Returns
    -------
    gdf : GeoDataFrame
        the points of non vegetation fire location
    """
    import geopandas as gpd
    from FireConsts import dirdata

    fnm = dirdata + "gas_flare/gas_flare_all_3571.gpkg"
    gdf = gpd.read_file(fnm)

    return gdf

def load_burnraster(year,region,sensor,xoff=0,yoff=0,xsize=None,ysize=None):
    '''get annual circumpolar modis or viirs burned/unburned tif
    optional: clip using UL pixel and pixel dimensions
    
    Parameters
    ----------
    year: int
    sensor: str, mcd64 for burned area only, viirs for burned area + active fires
    xoff, yoff: int, UL pixel coordinates for clipping
    xsize,ysize: int, pixels in x and y direction for clipping
        
    Returns
    -------
    arr: np.array, clipped image'''
    
    if sensor == 'mcd64':
        import gdal
        from FireConsts import mcd64dir
        
        fnm = mcd64dir + 'mcd64_' + str(year) + '.tif'
        ds = gdal.Open(fnm)
        #arr = ds.ReadAsArray(xsize=xsize, ysize=ysize)
        arr = ds.ReadAsArray(xoff=xoff, yoff=yoff, xsize=xsize, ysize=ysize)
    else:
        path = get_path(str(year),region)
        fnm = path + '/Summary/viirs_' + str(year) + '.tif'
        arr = fnm
        
    return arr

def get_any_shp(filename):
    ''' get shapefile of any region
    '''
    from FireConsts import dirdata
    import geopandas as gpd
    import os
    
    # find the california shapefile
    dirshape = os.path.join(dirdata,'shapefiles')
    statefnm = os.path.join(dirshape,filename)
    
    # read the geometry
    shp = gpd.read_file(statefnm).iloc[0].geometry
    
    return shp

def get_LCT_fnm(year):
    '''get land cover filename from year'''
    
    from FireConsts import lcdir
    
    # read ESA CCI LC 300m data
    if year > 2020: year = 2020 # no land cover data available after 2020
    if year in range(2016,2021): # years 2016-2020
        fnmLCT = lcdir+'C3S-LC-L4-LCCS-Map-300m-P1Y-'+str(year)+'-v2.1.1.nc'
    else: # years 2012-2015
        fnmLCT = lcdir+'ESACCI-LC-L4-LCCS-Map-300m-P1Y-'+str(year)+'-v2.0.7cds.nc'
    
    return fnmLCT

def get_LCT(locs,year):
    ''' Get land cover type for active fires
    
    Parameters
    ----------
    locs : list of lists (nx2)
        lat and lon values for each active fire detection
        
    Returns
    -------
    vLCT : list of ints
        land cover types for all input active fires
    '''
    from FireConsts import lc_dict
    import rasterio
    # from osgeo import gdal
    
    # get filename for lc dataset of the year
    fnmLCT = get_LCT_fnm(year)
    
    # read dataset
    sds = 'NETCDF:"'+fnmLCT+'":lccs_class'
    dataset = rasterio.open(sds) # rasterio version
    vLCT = dataset.sample(locs, indexes = 1)
    vLCT = [lc_dict[int(lc[0]/10)] for lc in vLCT] # simplify classes
    
    return vLCT

def get_peatstatus(locs):
    ''' Check if fire burns in peatland
    
    Parameters
    ----------
    locs : list of lists (nx2)
        lat and lon values for each active fire detection
        
    Returns
    -------
    vPEAT : list of ints
        peatland y/n for all fire pixels
    '''
    from FireConsts import dirdata, peatmap
    import rasterio
    
    # read dataset
    if peatmap == 'hugelius':
        fnm = dirdata + 'peat/peat_greater_10.tif'
    elif peatmap == 'peatml':
        fnm = dirdata + 'peat/peatml_greater_10.tif'
    elif peatmap == 'yu_hug':
        fnm = dirdata + 'peat/peat_yu_hug.tif'
    else: #use merged version
        fnm = dirdata + 'peat/peat_merge.tif'
    dataset = rasterio.open(fnm) # rasterio version
    vPEAT = dataset.sample(locs, indexes = 1)
    vPEAT = [lc[0] for lc in vPEAT] # list values
    
    return vPEAT

# ------------------------------------------------------------------------------
#%% read and load object, gdf and summary related files
# ------------------------------------------------------------------------------
def check_filefolder(fnm):
    ''' if the folder containing a file does not exist, make it

    Parameters
    ----------
    fnm : str
        file name
    '''
    import os

    # folder name
    dirnm = os.path.dirname(fnm)

    # create the folder if needed
    if not os.path.exists(dirnm):
        os.makedirs(dirnm)

def check_folder(dirnm):
    ''' if the folder does not exist, make it

    Parameters
    ----------
    dirfnm : str
        folder name
    '''
    import os
    
    # create the folder if needed
    if not os.path.exists(dirnm):
        os.makedirs(dirnm)

def get_path(year):
    ''' returns the base output directory for the year and region '''
    from FireConsts import dirpjdata
    path = dirpjdata+year
    return path

def correct_dtype(gdf,op='',geom=True):
    ''' correct the datatype for gdfs loaded from geojson files
    
    Parameters
    ----------
    gdf : gpd DataFrame
        the gdf directly read from the geojson file
    op : str
        the option for the geojson data -
        '': fire perimeter and attributes
        'FL': active fire line
        'NFP': new fire pixels
    
    Returns
    -------
    gdf : gpd DataFrame
        the gdf with correct datatype as given in FireConsts
    '''
    from FireConsts import dd
    
    # explicitly set the attributes data types
    if op == '':
        for v,tp in dd.items():
            if v == 'geometry' and not geom:
                continue
            gdf[v] = gdf[v].astype(tp)
    else:
        gdf['fireid'] = gdf['fireid'].astype('int')
        gdf['mergid'] = gdf['mergid'].astype('int')
    
    return gdf

def get_fobj_fnm(t):
    ''' Return the fire object pickle file name at a time step
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    
    Returns
    ----------
    fnm : str
        pickle file name
    '''
    from datetime import date
    
    d = date(*t[:-1])
    path = get_path(d.strftime('%Y'))
    fnm = path+'/Serialization/'+d.strftime('%Y%m%d')+t[-1]+'.pkl'
    return fnm

def check_fobj(t):
    ''' Check if the pickle file storing a daily allfires object exists

    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    '''
    import os

    fnm = get_fobj_fnm(t)
    return os.path.exists(fnm)

def save_fobj(allfires, t, activeonly=False):
    """ Save a daily allfires object to a pickle file

    Parameters
    ----------
    allfires : obj of Allfires class
        daily allfires
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    activeonly : bool
        the flag to save activeonly or full pickle file
    """

    import pickle
    from FireObj import Allfires

    if activeonly:
        # we only want to save active and sleeper fires
        allfires_out = Allfires(t)  # create a new allfires object
        # allfires_out.update_t(t)  # correct the time (previous time step was used in the intialization)
        fids_out = allfires.fids_active + allfires.fids_sleeper  # active fire or sleeper
        id_dict = []
        allfires_index = 0
        for fid in fids_out:
            allfires_out.fires.append(allfires.fires[fid])
            # allfires_out.fires[fid] = allfires.fires[fid]
            id_dict.append((allfires_index, fid))
            allfires_index += 1

        # copy the heritage from the original allfires object
        allfires_out.heritages = allfires.heritages
        allfires_out.mergedates = allfires.mergedates
        allfires_out.id_dict = id_dict
    else:
        allfires_out = allfires

    # get output file name
    fnm = get_fobj_fnm(t)

    # check folder
    check_filefolder(fnm)

    # save
    with open(fnm, "wb") as f:
        pickle.dump(allfires_out, f)

def load_fobj(t):
    ''' Load a daily allfires object from a pickle file

    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization

    Returns
    ----------
    data : obj of Allfires class
        daily allfires
    '''
    import pickle

    # get file name
    fnm = get_fobj_fnm(t)

    # load data
    with open(fnm,'rb') as f:
        data = pickle.load(f)
    return data

def get_gdfobj_fnm(t,op=''):
    ''' Return the fire object pickle file name at a time step
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    op : str
        the option for the geojson data -
        '': fire perimeter and attributes
        'FL': active fire line
        'NFP': new fire pixels
    Returns
    ----------
    fnm : str
        gdf file name
    '''
    from datetime import date

    d = date(*t[:-1])
    path = get_path(d.strftime('%Y'))
    if op == '':
        fnm = path+'/Snapshot/'+d.strftime('%Y%m%d')+t[-1]+'.gpkg'
    else:
        fnm = path+'/Snapshot/'+d.strftime('%Y%m%d')+t[-1]+'_'+op+'.gpkg'

    return fnm

def check_gdfobj(t,op=''):
    ''' Check if the gpkg file storing a daily allfires attributes exists

    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    op : str
        the option for the geojson data -
        '': fire perimeter and attributes
        'FL': active fire line
        'NFP': new fire pixels
    '''
    import os
    
    fnm = get_gdfobj_fnm(t,op=op)

    return os.path.exists(fnm)

def save_gdfobj(gdf,t,param='',fid='',op=''):
    ''' Save geopandas to a gpgk file

    Parameters
    ----------
    gdf : geopandas DataFrame
        daily allfires diagnostic parameters
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    param : string
        if empty: save daily allfires diagnostic dataframe
        if 'large': save daily single fire 
        if 'ign': save ignition layer
        if 'final': save final perimeters layer
    fid : int,
        fire id (needed only with param = 'large')
    '''

    # get file name
    if param == '':
        fnm = get_gdfobj_fnm(t,op=op)
    elif param == 'large':
        fnm = get_gdfobj_sf_fnm(t,fid,op=op)
    else:
        from datetime import date
        # get path
        d = date(*t[:-1])
        path = get_path(d.strftime('%Y'))
        # get file name
        fnm = path+'/Summary/'+param+d.strftime('%Y')+'.gpkg'
        
    # check folder
    check_filefolder(fnm)
    
    # reset the index column so a new fid column is created
    gdf = gdf.reset_index()

    # save file
    gdf.to_file(fnm, driver='GPKG')

def load_gdfobj(t='',op='',geom=True):
    ''' Load daily allfires diagnostic dataframe as geopandas gdf

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    Returns
    ----------
    gdf : geopandas DataFrame
        daily allfires diagnostic parameters
    '''
    import geopandas as gpd

    # get file name
    fnm = get_gdfobj_fnm(t,op=op)

    # read data as gpd DataFrame
    gdf = gpd.read_file(fnm)
    
    # correct the datatype
    gdf = correct_dtype(gdf,op=op,geom=geom)
    
    if not geom:
        gdf.drop(columns='geometry', inplace=True)

    # set fireid as the index
    gdf = gdf.set_index('fireid')

    return gdf

def save_FP_txt(df,t):
    
    # get filename of new fire pixels product
    fnm = get_gdfobj_fnm(t,op='NFP')
    fnm = fnm[:-4]+'txt' # change ending to txt
    
    # check folder
    check_filefolder(fnm)
    
    # write out
    if len(df) > 0:
        df.to_csv(fnm)

def load_lake_geoms(t,fid):
    ''' Load final perimeters as geopandas gdf

    Parameters
    ----------
    t: time tuple,
        needed for the year
    fid : int, 
        the fire fid
    Returns
    ----------
    geoms : tuple (geom, geom)
        geometry of all lake perimeters within the fire
    '''
    import geopandas as gpd
    
    year = str(t[0])
    path = get_path(year)
    
    # get file name
    fnm_lakes = path+'/Summary/lakes'+year+'.gpkg'
    
    # read data and extract target geometry
    gdf = gpd.read_file(fnm_lakes)
    gdf = gdf.set_index('mergid')
    if fid in gdf.index:
        gdf = gdf.loc[fid]
        geom_lakes = gdf['geometry']
    else:
        geom_lakes = None
        
    return geom_lakes

def save_merge_dict(mergedict,t):
    
    import pandas as pd
    
    year = str(t[0])
    path = get_path(year)
    
    # get file name
    fnm = path+'/Summary/merge_ted_dict.csv'
    
    # check folder
    check_filefolder(fnm)
    
    # write out
    df = pd.DataFrame.from_dict(mergedict, orient = 'index', 
                                columns = ['ted_year','ted_month','ted_day','ted_ampm'])
    df.index.names = ['fireid']
    df.to_csv(fnm)

def save_id_ted_dict(id_ted_dict,t):
    
    import pandas as pd
    
    year = str(t[0])
    path = get_path(year)
    
    # get file name
    fnm = path+'/Summary/id_ted_dict.csv'
    
    # check folder
    check_filefolder(fnm)
    
    # write out
    df = pd.DataFrame.from_dict(id_ted_dict, orient = 'index', 
                                columns = ['ted_year','ted_month','ted_day','ted_ampm'])
    df.index.names = ['fireid']
    df.to_csv(fnm)
    
def load_id_ted_dict(t):
    
    import pandas as pd
    
    year = str(t[0])
    path = get_path(year)
    
    # get file name
    fnm = path+'/Summary/id_ted_dict.csv'
    
    # read csv
    df = pd.read_csv(fnm)
    
    # combine the date and time columns to one attribute
    df['ted'] = df['ted_year'].astype(str)+df['ted_month'].astype(str).str.zfill(2)+df['ted_day'].astype(str).str.zfill(2)+df['ted_ampm']
    
    # turn back into dictionary
    id_ted_dict = dict(zip(df.fireid,df.ted))
    
    return id_ted_dict

def get_gdfobj_sf_fnm(t,fid,op=''):
    ''' Return the single fire fire object pickle file name at a time step
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    fid : int
        fire id
    Returns
    ----------
    fnm : str
        gdf file name
    '''
    from datetime import date
    
    d = date(*t[:-1])
    path = get_path(d.strftime('%Y'))
    
    if op == '':
        fnm = path+'/Largefire/F'+str(int(fid))+'_'+d.strftime('%Y%m%d')+t[-1]+'.gpkg'
    else:
        fnm = path+'/Largefire/F'+str(int(fid))+'_'+d.strftime('%Y%m%d')+t[-1]+'_'+op+'.gpkg'

    return fnm

def get_gdfobj_sf_fnms_year(year,fid,op=''):
    ''' Return the single fire fire object pickle file name at a time step
    Parameters
    ----------
    year : int
        year
    fid : int
        fire id
    Returns
    ----------
    fnms : list of str
        geojson file names
    '''
    from glob import glob
    
    path = get_path(year)
    if op == '':
        fnms = glob(path+'/Largefire/F'+str(int(fid))+'_??????????.gpkg')
    else:
        fnms = glob(path+'/Largefire/F'+str(int(fid))+'_??????????_'+op+'.gpkg')

    return fnms


def load_gdfobj_sf(t,fid,op=''):
    ''' Load single fire daily allfires diagnostic dataframe to a geojson file

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    fid : int
        fire id
    Returns
    ----------
    gdf : geopandas DataFrame
        time series of daily single fire diagnostic parameters
    '''
    import geopandas as gpd
    import pandas as pd

    # get file name
    fnm = get_gdfobj_sf_fnm(t,fid,op=op)

    # read data as gpd DataFrame
    gdf = gpd.read_file(fnm)

    # parse time step to Datetime and set as index
    gdf['index'] = pd.to_datetime(gdf['index'])
    gdf = gdf.set_index('index')

    return gdf

def get_summary_sf_fnm(t):
    ''' Return the fire summary file name at year end
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    Returns
    ----------
    fnm : str
        summary netcdf file name
    '''
    import os

    path = get_path(str(t[0]))
    fnm = os.path.join(path,'Summary','fsummary_sf.csv')

    # check folder
    check_filefolder(fnm)

    return fnm

def save_summary_sf(ds,t):
    ''' Save summary info as of t a netcdf file

    Parameters
    ----------
    ds : xarray dataset
        year end summary dataset
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    '''
    # get file name
    fnm = get_summary_sf_fnm(t)

    # check folder
    check_filefolder(fnm)

    # save netcdf file
    ds.to_csv(fnm)

def get_summary_fnm(t):
    ''' Return the fire summary file name at year end
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    Returns
    ----------
    fnm : str
        summary netcdf file name
    '''
    from datetime import date
    import os

    d = date(*t[:-1])
    path = get_path(str(t[0]))
    fnm = os.path.join(path,'Summary','fsummary_'+d.strftime('%Y%m%d')+t[-1]+'.nc')

    # check folder
    check_filefolder(fnm)

    return fnm

def get_summary_fnm_lt(t):
    ''' Return the latest time step before current time when summary file exists
    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization

    Returns
    ----------
    pt : tuple, (year, month, day, pmpm)
        the lastest time step with summary file
    '''
    import os
    from glob import glob
    import FireTime
    
    path = get_path(str(t[0]))
    
    # if there's no summary file for this year, return the first time step of the year
    fnms = glob(os.path.join(path,'Summary','fsummary_*.nc'))
    if len(fnms) == 0:
        return None

    # if there is, find the nearest one
    endloop = False
    pt = FireTime.t_nb(t,nb='previous')
    while endloop == False:
        if os.path.exists(get_summary_fnm(pt)):
            return pt
        else:
            pt = FireTime.t_nb(pt,nb='previous')
            if pt[0] != t[0]:
                return None

def check_summary(t):
    ''' Check if the summary file storing a daily allfires object exists

    Parameters
    ----------
    t : tuple, (year,month,day,str)
        the day and 'AM'|'PM' during the intialization
    '''
    import os

    # get file name
    fnm = get_summary_fnm(t)

    # return if it's present
    return os.path.exists(fnm)

def save_summary(ds,t):
    ''' Save summary info as of t a netcdf file

    Parameters
    ----------
    ds : xarray dataset
        year end summary dataset
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'
    '''
    # get file name
    fnm = get_summary_fnm(t)

    # check folder
    check_filefolder(fnm)

    # save netcdf file
    ds.to_netcdf(fnm)

def load_summary(t):
    ''' Load summary info from a netcdf file at t

    Parameters
    ----------
    t : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM'

    Returns
    -------
    ds : xarray dataset
        summary dataset
    '''
    import xarray as xr

    # get file name
    fnm = get_summary_fnm(t)

    # read data as xarray Dataset
    ds = xr.open_dataset(fnm)

    return ds

def save_summarycsv(df,year,op='heritage'):
    ''' save summary csv files

    Parameters
    ----------
    df : pandas DataFrame
        the data
    year : int
        year
    op : str
        option, 'heritage'|'large'
    '''
    import os
    
    path = get_path(str(year))
    
    fnm = os.path.join(path,'Summary','Flist_'+op+'_'+str(year)+'.csv')
    check_filefolder(fnm)

    df.to_csv(fnm)

def read_summarycsv(year,op='heritage'):
    ''' read summary csv files

    Parameters
    ----------
    year : int
        year
    op : str
        option, 'heritage'|'large'

    Returns
    -------
    df : pandas DataFrame
        the data
    '''
    import pandas as pd
    import os
    
    path = get_path(str(year))
    
    fnm = os.path.join(path,'Summary','Flist_'+op+'_'+str(year)+'.csv')
    df = pd.read_csv(fnm,index_col=0)
    return df

def get_lts_VNP14IMGTDL(year=None):
    from FireConsts import dirextdata
    from datetime import date
    from glob import glob
    import os

    if year == None:
        year = date.today().year

    dirFC = os.path.join(dirextdata,'VIIRS/VNP14IMGTDL') + '/'
    fnms = glob(os.path.join(dirFC,'SUOMI_VIIRS_C2_Global_VNP14IMGTDL_NRT_'+str(year)+'*.txt'))
    fnms.sort()
    DOY_lts = int(os.path.splitext(os.path.basename(fnms[-1]))[0][-3:])

    return DOY_lts

def get_lts_serialization(year=None):
    ''' get the time with lastest pkl data
    '''
    from datetime import date
    from glob import glob
    import os

    if year == None:
        year = date.today().year
    
    path = get_path(str(year))

    fnms = glob(path+'/Serialization/*.pkl')

    if len(fnms) > 0:
        fnms.sort()
        fnm_lts = os.path.basename(fnms[-1])

        lts = [int(fnm_lts[0:4]),int(fnm_lts[4:6]),int(fnm_lts[6:8]),fnm_lts[8:10]]
    else:
        lts = None

    return lts

#%% other functions related to read/write

def save2gtif(arr, outfile, cols, rows, geotrans, proj):
    '''write out a geotiff'''
    
    import gdal
    
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(outfile, cols, rows, 1, gdal.GDT_Byte)
    ds.SetGeoTransform(geotrans)
    ds.SetProjection(proj)
    band = ds.GetRasterBand(1)
    band.WriteArray(arr)

def geo_to_polar(lon_arr, lat_arr):
    '''transform lists of geographic lat lon coordinates to polar LAEA grid (projected)'''
    
    import numpy as np
    import pyproj
    
    proj4str = ("epsg:3571")
    p_modis_grid = pyproj.Proj(proj4str)
        
    x_arr, y_arr = p_modis_grid(lon_arr, lat_arr)
    x_arr = np.array(x_arr)
    y_arr = np.array(y_arr)
    
    return x_arr, y_arr

def polar_to_geo(x_arr, y_arr):
    '''transform lists of geographic lat lon coordinates to polar LAEA grid (projected)'''
    
    import numpy as np
    import pyproj
    
    proj4str = ("epsg:3571")
    p_modis_grid = pyproj.Proj(proj4str)
        
    lon_arr, lat_arr = p_modis_grid(x_arr, y_arr,inverse=True)
    lon_arr = np.array(lon_arr)
    lat_arr = np.array(lat_arr)
    
    return lon_arr, lat_arr

def world2Pixel(gt, Xgeo, Ygeo):
    ''' Uses a geomatrix (gdal.GetGeoTransform()) to calculate the pixel 
    location of a geospatial coordinate'''
    
    import numpy as np
    
    gt = list(gt)
    Xpx = np.rint((Xgeo - gt[0]) / gt[1]).astype(int)
    Ypx = np.rint((Ygeo - gt[3]) / gt[5]).astype(int)
    
    return (Xpx, Ypx)

def pixel2World(gt, Xpixel, Ypixel):
    '''Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate the
    geospatial coordinate of a pixel location'''
    
    Xgeo = gt[0] + Xpixel*gt[1] + Ypixel*gt[2]
    Ygeo = gt[3] + Xpixel*gt[4] + Ypixel*gt[5]
    
    return (Xgeo, Ygeo)
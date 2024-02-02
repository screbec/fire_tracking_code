# -*- coding: utf-8 -*-
"""
Writes a geodataframe of final fier perimeters for a calendar year
Created on Fri Jan 28 10:52:41 2022

@author: Rebecca C. Scholten
"""

def correct_final_ftype(gdf):
    '''correct the final fire type using the complete perimeter
    and supplementing the ESA-CCI land cover with CAVM for tundra regions'''
    
    import FireIO
    from FireConsts import dirdata
    import numpy as np
    import rasterio
    import rasterio.mask
    
    # grab year from the gdf
    year = gdf.iloc[0].tst_year
    
    # read CCI LC dataset
    fnmLCT = FireIO.get_LCT_fnm(year)
    ds = rasterio.open('NETCDF:"'+fnmLCT+'":lccs_class')
    
    # read tundra dataset
    ds_tundra = rasterio.open(dirdata+'cavm/cavm_coarse.tif')
    proj_tundra = ds_tundra.crs
    
    # read peat dataset
    fnm_peat = dirdata + 'peat/peat_greater_10.tif'
    ds_peat = rasterio.open(fnm_peat)
    
    # dictionaries for reclassifying
    tundra_dict = { # dictionary for naming of tundra classes
               1:'barren tundra', 
               2:'graminoid tundra', 
               3:'shrub tundra', 
               4:'wetland tundra'}
    boreal_dict = { # dictionary for naming of boreal classes
               0:'other',
               1:'cropland', 
               2:'forest', 
               3:'shrub/mosaic', 
               4:'grassland',
               5:'sparse',
               6:'urban'}
    lc_dict = { # dictionary for simplifying ESA CCI land cover types
               0:0, 20:0, 21:0, 22:0, # no data, water, ice, bare
               1:1, 2:1, 3:1, 4:1, # cropland
               5:2, 6:2, 7:2, 8:2, 9:2, # forest
               10:3, 11:3, 12:3, # shrubs & mosaic landscape
               13:4,  # grassland
               14:5, 15:5, # lichen, sparse vegetation
               16:2, 17:2, 18:3, # flooded forests and shrubs
               19:6} # urban

    
    # reproject gdf to cavm for clipping
    gdf_reproj_tundra = gdf.to_crs(proj_tundra)
    gdf_reproj_boreal = gdf.to_crs(4326)
    
    # add a new column for writeput to the dataframe
    gdf['lcc_final'] = None
    
    for fire in gdf.index:
        
        # get geom in wgs84
        geom = gdf_reproj_boreal.geometry[fire]
        
        # check for peat percentage
        peat_arr, peat_trans = rasterio.mask.mask(ds_peat, [geom],
                                                  all_touched = True, crop=True, nodata=255)
        peat_arr = peat_arr.astype(float)
        peat_arr[peat_arr == 255] = np.nan
        gdf.loc[fire,'peat'] = np.nanmean(peat_arr)
        
        boreal = True
        perim_tundra = gdf_reproj_tundra.geometry[fire]      
        # check if the fire is a tundra fire
        try:
            tundra_arr, trans_tundra = rasterio.mask.mask(ds_tundra, [perim_tundra], 
                                                          all_touched = True, crop=True, nodata=0)
            if np.sum(tundra_arr) > 0:
                # we have a tundra fire
                values, counts = np.unique(tundra_arr, return_counts=True)
                lcc = values[np.argmax(counts)]
                if lcc > 0:
                    boreal = False
                    lcc = tundra_dict[lcc]
        except:
            pass
        # if the fire is not in tundra, we take the ESA CCI lc class
        if boreal:
            arr, trans = rasterio.mask.mask(ds, [geom], 
                                            all_touched = True, crop=True, nodata=255)
            arr = (arr/10).astype(int)
            values, counts = np.unique(arr[~(arr==25)], return_counts=True)
            lcc = values[np.argmax(counts)]
            lcc = lc_dict[lcc]
            lcc = boreal_dict[lcc]
            
        # assign land cover class to datbase entry
        gdf.loc[fire,'lcc_final'] = lcc
    
    return gdf

def find_all_end(tst,ted):
    ''' find all final perimeters points in the half-daily snapshots and
    save them to one gdf
    
    Parameters
    ----------
    tst : Start date of fire season
    ted: End date of fire season
    
    Returns
    -------
    gdf_all : geopandas DataFrame
        the gdf containing half daily fire basic attributes and final perimeter
    '''
    import FireTime,FireIO
    import geopandas as gpd
    
    # initialise list of already checked fire ids
    checked_ids = []
    id_ted_dict = {}
    
    endloop = False  # flag to control the ending of the loop
    creategdf = True # needed since ted cannot be used anymore
    t = list(ted)    # t is the time (year,month,day,ampm) for each step
    while endloop == False:
        #print(t)
        # read daily gdf
        if FireIO.check_gdfobj(t,op=''):
            gdf = FireIO.load_gdfobj(t)
            gdf_active = gdf[gdf.isactive == 1]
            
            # append daily row to gdf_all
            if creategdf:
                gdf_all = gdf_active
                creategdf = False
            else:
                # exclude ids that have an older active perimeter
                gdf_active = gdf_active[~gdf_active['mergid'].isin(checked_ids)]
                gdf_all = gdf_all.append(gdf_active)
            
            # add ids that have been written out to checked_ids, these will be skipped next time
            if len(gdf_active) > 0:
                for id in gdf_active['mergid'].tolist():
                    id_ted_dict[id] = tuple(t)
                checked_ids = list(set(id_ted_dict))
        
        #  - if t reaches tst, set endloop to True to stop the loop
        if FireTime.t_dif(t,tst)==0:
            endloop = True
        
        #  - update t with the previous time stamp
        t = FireTime.t_nb(t,nb='previous')
    
    # clip inner lakes
    for fid in gdf_all.index:
        if gdf_all.loc[fid].lake_area > 0: # check if there are lakes within before loading lake file
            mergid = gdf_all.loc[fid].mergid
            lakediss = FireIO.load_lake_geoms(tst,mergid)
            geom = gdf_all.loc[[fid]].iloc[0].geometry
            geom = geom.difference(lakediss)
            gdf_all.loc[[fid],'geometry'] = gpd.GeoSeries([geom],crs=3571).values
    
    # lastly we also save the id_dict for the summary file later
    FireIO.save_id_ted_dict(id_ted_dict,ted)
    
    return id_ted_dict, gdf_all



def make_gdf_fperim_simple(allfires, active_only=True):
    ''' Create geopandas DataFrame for fire basic attributes and fire perimeter (hull).
            Use saved gdf files for previous time step and update active fires only.

    Parameters
    ----------
    allfires : Allfires object
        the Allfires object to be used to modify gdf

    Returns
    -------
    gdf : geopandas DataFrame
        the gdf containing half daily fire basic attributes and fire perimeter
    '''
    import geopandas as gpd

    # initialize the GeoDataFrame with fire attribute names (in dd) as columns
    gdf = gpd.GeoDataFrame(columns=['fireid','mergid'],crs='epsg:3571', geometry=[])
    heritage = dict(allfires.heritages)
    if active_only:
        firecounter = range(allfires.number_of_activefires)
    else:
        firecounter = allfires.fids_valid

    # update the hull of each active fire as the geometry column
    for fire_no in firecounter:
        fid = allfires.fires[fire_no].id
        # print(fid)
        gdf.loc[fid,'fireid'] = fid
        if fid in heritage.keys():
            gdf.loc[fid,'mergid'] = heritage[fid]
        else:
            gdf.loc[fid,'mergid'] = fid
        fhull = allfires.fires[fire_no].hull
        if fhull is not None:
            if fhull.geom_type == 'MultiPolygon':
                gdf.loc[fid,'geometry'] = gpd.GeoDataFrame(geometry=[fhull]).geometry.values
            else:
                gdf.loc[fid,'geometry'] = fhull
    gdf = gdf.dissolve(by = 'mergid')
    # make sure the data types are correct
    gdf['fireid'] = gdf.index

    return gdf

  
def save_gdf_trng(tst,ted):
    ''' Wrapper to create and save all ignitions as gdf
        1) final perimeters from VIIRS only
        2) final perimeters combining VIIRS and MCD64
        3) final perimeters with clipped lakes

    Parameters
    ----------
    tst : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at start time
    ted : tuple, (int,int,int,str)
        the year, month, day and 'AM'|'PM' at end time
    '''
    import FireIO
    
    # find all final perimeters and write out to gdf
    id_ted_dict,gdf = find_all_end(tst,ted)
    gdf = correct_final_ftype(gdf)
    FireIO.save_gdfobj(gdf,tst,param='final_viirs')
    


if __name__ == "__main__":
    ''' The main code to record time series of geojson data for a fire
    '''

    import time
    t1 = time.time()
    # set the start and end time
    tst=(2020,6,1,'AM')
    ted=(2020,8,31,'PM')

    # for each day during the period,
    # create and save geojson files for temporal evolution of large fires
    save_gdf_trng(tst=tst,ted=ted)

    t2 = time.time()
    print(f'{(t2-t1)/60.} minutes used to run code')

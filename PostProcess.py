# -*- coding: utf-8 -*-
"""
Processing steps to remove surface water bodies from fire perimeters
Created on Wed Jan 19 14:14:52 2022

@author: Rebecca Scholten
"""

def retrieve_water_tile(ext, year):
    '''Retrieve the Global Surface Water tile number
    from the bounding box of a geometry (in lat lon!)
    
    Parameters
    ----------
    ext : tuple
        the extent of the geometry of the fire to process
    year : year
        year from which GSW data should be taken
    
    Returns
    -------
    tiles : tuple of lists
        lists containing all overlapping x and y tiles
    '''
    import math
    import numpy as np
    
    # starting lats and lons for tiles
    tilesx = np.arange(-180, 180, 10) 
    tilesy = np.flip(np.arange(40, 80, 10))
    
    # extract bounding box of fire
    minx, miny, maxx, maxy = ext
    
    # extract nearest start lat and lon form bounding box
    minx = math.floor(minx/10)*10 
    maxx = math.floor(maxx/10)*10 
    miny = math.floor(miny/10)*10 
    maxy = math.floor(maxy/10)*10 
    
    # find indices
    tilex = np.where(tilesx == minx)[0].tolist()
    tiley = np.where(tilesy == miny)[0].tolist()
    if maxx - minx > 0:
        tilex.append(np.where(tilesx == maxx)[0].tolist()[0])
    if maxy - miny > 0:
        tiley.append(np.where(tilesy == maxy)[0].tolist()[0])
    
    if len(tilex) > 0 and len(tiley) > 0:
        # some fires extrend outside of the water dataset (north of 80N in Greenland)
        tiles = (tilex, tiley)
    
        lake_files = tile_2_fnm(tiles, year)
    else:
        lake_files = None
    
    return lake_files

def tile_2_fnm(tiles, year):
    '''Retrieve filenames of all GSW tiles returned by 
    retrieval_water_tile
    
    Parameters
    ----------
    tiles : tuple of 2 lists
        2 lists containing the tile numbers in x and y direction
    year : year
        year from which GSW data should be taken
    
    Returns
    -------
    fnms : list of strings
        lists containing paths of all overlapping tile files
    '''
    from FireConsts import lakedir
    import itertools, os
    
    # there is no surface water informaion for 2021, so we replace that with 2020
    if year > 2020:
        year = 2020
    
    tilex, tiley = tiles
    tilex = [str(tile*4).zfill(3) for tile in tilex]
    tiley = [str(tile*4).zfill(2) for tile in tiley]
    all_combinations = [list(zip(each_permutation, tiley)) for each_permutation in itertools.permutations(tilex, len(tiley))]
    
    fnms = []
    basepath = lakedir+str(year)+'/GSW_300m_'
    for tile in all_combinations:
        
        fnm = basepath+str(tile[0][1])+'_'+str(tile[0][0])+'.gpkg'
        # check if path exists
        if os.path.exists(fnm):
            fnms.append(fnm)
        # missing file means the tile does not contain any lakes
        # (or has not been processed! make sure to run poygonise_water_test on the year before)
    
    return fnms

def dissolve_lake_geoms(geom, ext, year):
    '''Read all lake geometries within the geometry boundaries and 
    if necessary dissolve them to one MultiPolygon
    
    Parameters
    ----------
    geom : geometry
        the geometry of the fire to process
    year : year
        year from which GSW data should be taken
    
    Returns
    -------
    lakediss : geopandas dataframe
        can be empty (no overlap), or contains one entry (index = 1) with a 
        Polygon or MultiPolygon
    '''
    import geopandas as gpd
    import pandas as pd
    from shapely.geometry import MultiPolygon, Polygon
    
    # retrieve GSW tile paths
    fnm_lakes = retrieve_water_tile(ext, year) # retrieve water tile paths
    
    # load and dissolve water features
    if fnm_lakes:
        if len(fnm_lakes) > 1:
            temp = []
            for lake in fnm_lakes:
                temp.append(gpd.read_file(lake, mask=geom))
            lakes_test = pd.concat(temp, ignore_index = True)
        else:
            lakes_test = gpd.read_file(fnm_lakes[0], mask=geom)
        
        lakes_test = lakes_test.assign(diss = 1)
        lakediss = lakes_test.dissolve(by = 'diss')
        
        # at the edges of files we can have artificial holes in the lakes due to mosaicking
        if len(fnm_lakes) > 1:
            if len(lakediss) > 0:
                geom_orig = lakediss.loc[1].geometry
                # cnt = 0
                # for p in geom_orig:
                #     if len(p.interiors)>0:
                #         print(cnt)
                #         print(p)
                #     cnt+=1
                if type(geom_orig) == Polygon:
                    geom_fix = Polygon(geom_orig.exterior)
                else:
                    geom_fix = MultiPolygon(Polygon(p.exterior) for p in geom_orig)
                lakediss.set_geometry([geom_fix], inplace = True)
            
    else:
        lakediss = None
    
    return lakediss
   
def final_lakes(gdf, t):
    '''Wrapper to creade a geodataframe of final surface water elements
    
    Parameters
    ----------
    gdf : pandas geodataframe
        geodataframe of final fire perimeters
    t : tuple
        time tuple
    
    Returns
    -------
    gdf_lakes: geopandas geodataframe
        geodataframe containing only the lakes fore each fire perimeter id
        (this is used for clipping daily large fire perimeters)
    '''
    
    import geopandas as gpd
    from shapely.geometry import box
    import FireIO
    
    # extract year
    year = t[0]
    lake_list = []
    
    # loop over all final perimeters
    for row in gdf.itertuples():
        fireid = row.mergid
        # print(fid)
        geom = row.geometry
        # fetch bounding box in lat lon
        ext = gpd.GeoSeries([geom], crs=3571).to_crs(4326).geometry[0].bounds
        
        # read and dissolve lakes for the geom
        lakediss = dissolve_lake_geoms(geom, ext, year) 
        
        if not isinstance(lakediss, type(None)):
            if len(lakediss) > 0:
                lake_geom = box(*geom.bounds).buffer(2000).intersection(lakediss.loc[1].geometry)
                lake_list.append((fireid, lake_geom))
    
    # turn lake geom and fire id into geopandas dataframe
    fids, geoms = zip(*lake_list)
    d = {'fireid': fids, 'geometry': geoms, 'mergid': fids}
    gdf_lakes = gpd.GeoDataFrame(d, crs="EPSG:3571")
    
    # save to file
    FireIO.save_gdfobj(gdf_lakes,t,param='lakes')
    
    return gdf_lakes

def clip_lakes_1fire_outer(geom, lakediss, lake_idx=None):
    '''Clip surface water areas from a fire geometry
    (only corrects the outer perimete)
    
    Parameters
    ----------
    geom : geometry
        the geometry of the fire to process
    lakediss : geopandas dataframe
        can be empty (no overlap), or contains one entry (index = 1) with a 
        Polygon or MultiPolygon
    
    Returns
    -------
    geom: clipped geometry (only outer lakes clipped)
    lake_area: area within fire covered by lakes
    lake_no: number of lakes within fire
    intsct_length: length of fire perimeter that ends at lake
    intsct_length: total length of fire perimeter bordering lake areas
    '''
    import shapely
    import FireClustering
    
    # initialise
    intsct_length = 0       # length of perimeter bordering lake for this time step
    intsct_length_tot = 0   # total length of intersection incl inner lakes
    lake_area = 0           # total lake area within the fire
    lake_no = 0             # number of lakes within perimeter
    
    # if its a polygon and not multipolygon we put it in a list
    if not isinstance(lakediss, shapely.geometry.multipolygon.MultiPolygon):
        lakediss = shapely.geometry.multipolygon.MultiPolygon([lakediss])
    
    # check if an rtree for fast search of intersecting lakes is given
    if not lake_idx:
        lake_idx = FireClustering.build_rtree(lakediss) # else create one
    
    # extract potential lakes in geometry bounding box
    lake_fids = FireClustering.idx_intersection(lake_idx, geom.bounds)
    
    # loop over lakes
    for i in lake_fids:
        lake = lakediss[i]
        if lake.intersects(geom):
            # if lakes are within fire scar we just compute lake statistics
            if lake.within(geom):
                lake_no += 1
                lake_area += lake.area
                intsct_length_tot += lake.length
            # if lakes are intersecting fire scar we clip
            else:
                geom = geom.difference(lake)
                outer_perim = lake.intersection(geom).length
                intsct_length += outer_perim
                intsct_length_tot += outer_perim
    
    return geom, lake_area, lake_no, intsct_length, intsct_length_tot

def clip_lakes_1fire(geom, lakediss):
    '''Clip all surface water elements from a fire geometry
    (including those completely within the perimeter)
    
    Parameters
    ----------
    geom : geometry
        the geometry of the fire to process
    lakediss : geopandas dataframe
        can be empty (no overlap), or contains one entry (index = 1) with a 
        Polygon or MultiPolygon
    lakes_only: bool,
        return only lake r
    
    Returns
    -------
    tiles : tuple of lists
        lists containing all overlapping x and y tiles
    '''
    # import pyproj
    import copy
    
    # clip lakes
    geom = geom.difference(lakediss)
    
    geom_nolakes = copy.deepcopy(geom) # geometry without lakes
    #geod = pyproj.Geod(ellps='WGS84') # geoid used for computing distances in m (if EPSG:4326)
    intsct_length = 0 # length of perimeter bordering lake
    lake_area = 0 # total lake area within the fire
    
    # additional attributes
    intsct = geom_nolakes.difference(geom) # linestring of intersection
    
    # if not intsct.is_empty: # if there are lakes in the area
    if intsct.geom_type == 'Polygon': # if only one lake is present we put it in a list for the loop
        intsct = [intsct]
    # check which lakes are inside firescar and which are outside
    for i in range(len(intsct)):
        lake = intsct[i]
        if lake.within(geom_nolakes): # if lakes are within fire scar
            #temp = geod.geometry_area_perimeter(lake)
            #lake_area += abs(temp[0])
            # intsct_length += temp[1]
            lake_area = lake.area
            intsct_length += lake.length
        else: # if lakes are bordering fire scar
            true_intsct = lake.intersection(geom)
            intsct_length += true_intsct.length
            # intsct_length += geod.geometry_area_perimeter(true_intsct)[1]
    
    out = (geom, lake_area, intsct_length)
    
    return out

if __name__ == "__main__":
    ''' The main code to record daily geojson data for a time period
    '''
    pass
    

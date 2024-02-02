""" FireObj
This is the module containing the object definitions used in the system

This module include:
1. MODULES AND UTILITY FUNCTIONS
    a. Other project modules used in this module
    b. Object related utility functions (related to time step)
2. FOUR LAYERS OF OBJECTS
    a. Allfires:  the class of all fire events
    b. Fire:      the class of a fire event
    c. Cluster:   the class of active fire pixel cluster (only for supporting)
    d. FirePixel: the class of an active fire pixel
"""

# ------------------------------------------------------------------------------
# 1. MODULES AND UTILITY FUNCTIONS
# ------------------------------------------------------------------------------

# a. Project modules used
import numpy as np
# from pyproj import Geod
import FireIO
import FireVector
from FireTime import t_dif
from FireConsts import maxoffdays,sleepdays,fpbuffer,flbuffer

# ------------------------------------------------------------------------------
# 2. FOUR LAYERS OF OBJECTS
# ------------------------------------------------------------------------------

# a. Object - Allfires
class Allfires:
    """ Class of allfire events at a particular time step
    """

    # initilization
    def __init__(self,t):
        ''' Initiate the object with current time
        Parameters
        ----------
        t : tuple, (int,int,int,str)
            the year, month, day and 'AM'|'PM'
        '''
        # time
        self.t = t # itialize the object with previous time step

        # a list of Fire objects
        self.fires = []

        # Initialise lists of fire ids which change in each time step
        self.fids_expanded = []   # a list of ids for fires with expansion at current time step
        self.fids_new = []        # a list of ids for all new fires formed at current time step
        self.fids_merged = []     # a list of ids for fires with merging at current time step
        self.fids_invalid = []    # a list of ids for fires invalidated at current time step

        # cumulative recordings
        self.heritages = []       # a list of fire heritage relationships (source, target)
        self.mergedates = []      # a list of merge dates
        self.id_dict = []         # this list relates the list position of the fire in the allfires object to the fire id
                                  # (list_index, fireid) (list position can be variable when writing out only active fires)

    # properties
    @property
    def cday(self):
        ''' Datetime date of current time step
        '''
        from datetime import date
        return date(*self.t[:-1])

    @property
    def ampm(self):
        ''' Ampm indicator of current time step
        Parameters
        ----------
        ampm : str, 'AM'|'PM'
           the ampm option calculated from t
        '''
        return self.t[-1]

    @property
    def number_of_fires(self):
        ''' Total number of fires (active and inactive) at this time step
        '''
        return len(self.fires)

    @property
    def fids_active(self):
        ''' List of active fire ids
        '''
        return [f.id for f in self.fires if f.isactive is True]
    
    @property
    def fids_sleeper(self):
        ''' List of fire ids that may reactivate
        '''
        return [f.id for f in self.fires if f.mayreactivate is True]

    @property
    def number_of_activefires(self):
        ''' Total number of active fires at this time step
        '''
        return len(self.fids_active)

    @property
    def activefires(self):
        ''' List of active fires
        '''
        return [self.fires[fid] for fid in self.fids_active]

    @property
    def fids_valid(self):
        ''' List of valid (non-invalid) fire ids
        '''
        return [f.id for f in self.fires if f.invalid is False]

    @property
    def number_of_validfires(self):
        ''' Total number of valid fires at this time step
        '''
        return len(self.fids_valid)

    @property
    def validfires(self):
        ''' List of valid fires
        '''
        return [self.fires[fid] for fid in self.fids_valid]

    @property
    def fids_updated(self):
        ''' List of fire id which is updated at this time step
            (expanded, new, merged, invalid)
        '''
        fids_updated = list(set(self.fids_expanded+self.fids_new+
                                self.fids_merged+self.fids_invalid))
        return fids_updated

    @property
    def fids_ne(self):
        ''' List of fire id which is newly formed or expanded
               at this time step
        '''
        fids_ne = sorted(set(self.fids_expanded+self.fids_new))
        return fids_ne

    # functions to be run before tracking VIIRS active fire pixels at each time step
    def update_t(self,t):
        ''' Update the time (cday and ampm) for the Allfire object.
        Parameters
        ----------
        t : tuple, (int,int,int,str)
            the year, month, day and 'AM'|'PM'
        '''
        self.t = list(t)  # current date and ampm

    def update_t_allfires(self,t):
        ''' Update the time (t) for each Fire object in the Allfire object.
        Parameters
        ----------
        t : tuple, (int,int,int,str)
            the year, month, day and 'AM'|'PM'
        '''
        for f in self.fires:
            f.t = list(t)

    def reset_newpixels(self):
        ''' Reset newpixels to [] for each fire.
        '''
        for f in self.fires:
            f.newpixels = []

    def cleanup(self,t):
        ''' Clean up Allfires obj at each time step
        - update t (for allfires and all fires)
        - clean up lists to record fire changes
        - clean up newpixels for each fire
        '''
        # time updated to t
        self.update_t(t)           # update t for allfires
        self.update_t_allfires(t)  # update t

        # reset the fids used to record changes
        self.fids_expanded = []   # a list of ids for fires with expansion at current time step
        self.fids_new = []        # a list of ids for all new fires formed at current time step
        self.fids_merged = []     # a list of ids for fires with merging at current time step
        self.fids_invalid = []    # a list of ids for fires invalidated at current time step
        
        # reset the new pixels
        self.reset_newpixels()
        


    # functions to be run after tracking VIIRS active fire pixels at each time step
    def record_fids_change(self,
         fids_expanded=None, fids_new=None, fids_merged=None, fids_invalid=None):
        ''' Update the list of ids for fires with changes at current time step.
        Parameters
        ----------
        fids_expanded : list
            ids of expanded fires
        fids_new : list
            ids of new formed fires
        fids_merged : list
            ids of fires with other fires merging to them
        fids_invalid : list
            ids of fires invalidated (due to merging with other fires)
        '''
        if fids_expanded:
            self.fids_expanded = fids_expanded     # expanded fires
        if fids_new:
            self.fids_new = fids_new               # new formed fires
        if fids_merged:
            self.fids_merged = fids_merged         # fires with merging with other fires
        if fids_invalid:
            self.fids_invalid = fids_invalid       # fires invalidated due to merging with other fires

    def invalidate_statfires(self):
        ''' Original: If pixel density of an active fire is too large, assume it's static
                fires and invalidate it.
            Now: 
                (1) invalidate fires that only burn one pixel at a single time step
                (2) invalidate fires if they have been inactive longer than sleepertime
        '''
        for f in self.validfires:
            if (not f.isactive) & (not f.mayreactivate) & (f.duration < 1) & (f.n_pixels == 1):
            # if (f.pixden > 20) & (f.farea < 20):
                # invalidate the fire
                f.invalid = True

                # add the fire id into the fids_invalid list
                self.fids_invalid.append(f.id)
    
    def record_flines(self):
        '''adds the current fire line to the fline records
        '''
        for f in self.validfires:
            if f.ftype in [0,1,5,7,8]: # fires in non-peat and non-forest fire fronts can't sleep
                fline = f.fline_prior_list + [None]
            else:
                fline = f.fline_prior_list + [f.fline]
            if len(fline) > sleepdays*2:
                fline = fline[1:]
            f.fline_prior_list = fline
                
    def record_spreadrates(self):
        '''records average spread rate of each fire for each time step
        ----- this is not currently in use -----
        '''
        for f in self.activefires:
            f.spreadrate_ts.append(f.spreadavg)

# b. Object - Fire
class Fire:
    """ Class of a single fire event at a particular time step
    """

    # initilization
    def __init__(self, id, t, pixels):
        ''' Initialize Fire class with active fire pixels locations and t. This
                is only called when fire clusters forming a new Fire object.
        Parameters
        ----------
        id : int
            fire id
        t : tuple, (int,int,int,str)
            the year, month, day and 'AM'|'PM'
        pixels : list (nx5)
            x coord, y coord, line, sample, and FRP values of active fire pixels
        '''
        # initialize fire id 
        self.id = id
        self.mergeid = id

        # initialize current time, fire start time, and fire final time
        tlist = list(t)
        self.t = tlist  # current time
        self.t_st = tlist
        self.t_ed = tlist

        # initialize pixels
        fpixels = pixels
        self.pixels = fpixels      # all pixels
        self.newpixels = fpixels   # new detected pixels
        self.ignpixels = fpixels   # pixels at ignition

        # initialize hull using the pixels
        locs = [p.loc for p in fpixels]  # list of [x,y] coordinates
        locs_geo = [p.loc_geo for p in fpixels]  # list of [lon,lat] coordinates
        hull = FireVector.cal_hull(locs) # the hull from all locs
        self.hull = hull   # record hull
        self.prior_hull = None # hull of prior time step, used for spread rate computation
        self.spreadrate_ts = [] # time series of average spread rates for the fire
        self.fline_prior_list = []
        
        # initialize the exterior pixels (pixels within the inward extbuffer of
        #    the hull; used for saving time for hull calculation of large fires)
        self.extpixels = FireVector.cal_extpixels(fpixels,hull)

        # always set validate at initialization
        self.invalid = False
        
        # get and record fm1000 value at ignition (at the time of initilization)
        # self.stFM1000 = FireIO.get_stFM1000(hull,locs,tlist)

        # LCTmax
        vLCT = FireIO.get_LCT(locs_geo, self.t[0])
        self.LCTmax = max(set(vLCT), key = vLCT.count)
        
        # peat status
        vPEAT = FireIO.get_peatstatus(locs_geo)
        self.peat = 1 if 1 in vPEAT else 0 # if fire is intersecting at least 1 10x10km2 peat tile

    # properties
    @property
    def cday(self):
        ''' Current day (datetime date)
        '''
        from datetime import date
        return date(*self.t[:-1])
    
    @property
    def cdoy(self):
        ''' Current day of year (Julian day, int)
        '''
        return self.cday.timetuple().tm_yday

    @property
    def ampm(self):
        ''' Current ampm flag, 'AM'|'PM'
        '''
        return self.t[-1]

    @property
    def duration(self):
        ''' Time difference between first and last active fire detection
        '''
        duration = t_dif(self.t_st,self.t_ed) + 0.5
        return duration

    @property
    def t_inactive(self):
        ''' Time difference between current time and the last active fire detection
        '''
        t_inactive = t_dif(self.t_ed,self.t)
        return t_inactive

    @property
    def isactive(self):
        ''' Fire active status
        '''
        # invalidated fires are always inactive
        if self.invalid:
            return False
        # otherwise, set to True if no new pixels detected for 5 consecutive days
        return (self.t_inactive <= maxoffdays[self.ftype])
    
    @property
    def mayreactivate(self):
        ''' Fire sleeper status
        '''
        cond1 = self.invalid # invalidated fires are always inactive
        cond2 = self.isactive # not active anymore
        cond3 = self.t_inactive > sleepdays # incative time longer than sleeper days
        
        # flines = [fline for fline in self.fline_prior_list if fline is not None]
        # cond4 = len(flines) == 0 # at least one valid fire line in correct fuel type in sleeperdays
        
        if cond1 or cond2 or cond3:
            return False
        else:
            return True
    
    @property
    def isignition(self):
        ''' Is the current timestep the ignition?
        '''
        if len(self.newpixels) == 0: # in this case t_st = t_ed because fire is inactive
            ign = 0
        else:
            ign = (t_dif(self.t_st,self.t_ed) == 0)*1
        
        return ign
        
    @property
    def locs(self):
        ''' List of fire pixel locations (x,y)
        '''
        return [p.loc for p in self.pixels]

    @property
    def n_pixels(self):
        ''' Total number of fire pixels'''
        return len(self.pixels)

    @property
    def newlocs(self):
        ''' List of new fire pixels locations (x,y)
        '''
        return [p.loc for p in self.newpixels]
    
    @property
    def newlocs_geo(self):
        ''' List of new fire pixels locations (lon,lat)
        '''
        return [p.loc_geo for p in self.newpixels]
    
    @property
    def newpixelatts(self):
        ''' List of new fire pixels attriputes (frp, (DS))
        '''
        return [(p.frp, p.DS, p.DT, p.datetime, p.sat) for p in self.newpixels]

    @property
    def n_newpixels(self):
        ''' Total number of new fire pixels
        '''
        return len(self.newpixels)

    @property
    def extlocs(self):
        ''' List of exterior fire pixel locations (x,y)
        '''
        return [p.loc for p in self.extpixels]

    @property
    def n_extpixels(self):
        ''' Total number of exterior fire pixels
        '''
        return len(self.extpixels)

    @property
    def ignlocs(self):
        ''' List of fire pixel locations (x,y) at ignition time step
        '''
        return [p.loc for p in self.ignpixels]

    @property
    def n_ignpixels(self):
        ''' Total number of ignition fire pixels
        '''
        return len(self.ignpixels)
    
    @property
    def farea(self):
        ''' Fire spatail size of the fire event (km2)
        '''
        import FireVector
        from FireConsts import area_VI

        # get hull
        fhull = self.hull

        # If no hull, return area calculated from number of pixels
        if fhull is None:
            return self.n_pixels * area_VI
        # otherwise, use calConcHarea to calculate area,
        #   but no smaller than area_VI (sometimes calculated hull area is very mall)
        else:
            return max(FireVector.calConcHarea(fhull),area_VI)
    
    @property
    def fperim(self):
        ''' Perimeter length of fire hull
        '''
        # get hull
        fhull = self.hull

        if fhull is None:  # if no hull, return zero
            perim = 0
        else:  # otherwise, use the hull length
            perim = fhull.length
        return perim
    
    @property
    def pixden(self):
        ''' Fire pixel density (number of pixels per km2 fire area)
        '''
        farea = self.farea
        if farea > 0:
            return self.n_pixels/farea
        else:
            return 0

    @property
    def meanFRP(self):
        ''' Mean FRP of the new fire pixels
        '''
        frps = [p.frp for p in self.newpixels]
        if len(frps) > 0:
            m = sum(frps)/len(frps)
        else:
            m = 0
        return m
    
    @property
    def FRP95(self):
        ''' Mean FRP of the new fire pixels
        '''
        frps = [p.frp for p in self.newpixels]
        if len(frps) > 0:
            m = np.percentile(frps,95)
        else:
            m = 0
        return m

    @property
    def ftype(self):
        ''' Fire type (as defined in FireConsts) derived using LCTmax and stFM1000
        '''
        # get the dominant land cover type
        LCTmax = self.LCTmax
        
        # determine fire type using land cover and peatland status
        peatstatus = self.peat
        
        if peatstatus == 1:
            peat_dict = {0:0,1:1,2:2,3:2,4:2,5:2,6:2,7:4,8:6,9:4,10:4,11:2,12:4,13:8}
        else:
            peat_dict = {0:0,1:1,2:3,3:3,4:3,5:3,6:3,7:5,8:7,9:5,10:5,11:3,12:5,13:8}
        
        return peat_dict[LCTmax]

    @property
    def ftypename(self):
        ''' Fire type name
        '''
        from FireConsts import FTYP
        return FTYP[self.ftype]

    @property
    def flinepixels(self):
        ''' List of all fire pixels near the fire perimeter (fine line pixels)
        '''
        from shapely.geometry import Point, MultiLineString
        from FireConsts import VIIRSbuf
        
        # get pixels of last active fire detection
        nps = self.newpixels
        
        # get hull
        fhull = self.hull
        
        if fhull is None: # if no hull, return empty list
            return []
        else:  # otherwise, extract the pixels nearl the hull
            nps_buf = [FireVector.addbuffer(Point(p.x,p.y), VIIRSbuf+fpbuffer) for p in nps]
            # if hull is a polygon, return new pixels near the hull
            if fhull.type == 'Polygon':
                lr = fhull.exterior
                return [p for i,p in enumerate(nps) if lr.intersects(nps_buf[i])]
            
            # if hull is a multipolygon, return new pixels near the hull
            elif fhull.type == 'MultiPolygon':
                lr = MultiLineString([x.exterior for x in fhull])
                # mlr = FireVector.addbuffer(mlr,fpbuffer)
                return [p for i,p in enumerate(nps) if lr.intersects(nps_buf[i])]
            
            # instead of buffer: compute all distances
            # dist = [lr.distance(Point(p.x,p.y)) for p in nps]
            # return [p for i,p in enumerate(nps) if dist[i] <= fpbuffer+VIIRSbuf]
        
        

    @property
    def fline(self):
        ''' Active fire line MultiLineString shape (segment of fire perimeter with active fires nearby)
        '''
        from shapely.geometry import MultiLineString#,Polygon,Point,MultiPoint
        from FireConsts import VIIRSbuf
        
        if len(self.flinepixels)==0: # if for whatever reason we don't have an active fire line
            return None
        
        # get fireline pixel locations
        flinelocs = [p.loc for p in self.flinepixels]
        flinelocsMP = FireVector.doMultP(flinelocs,VIIRSbuf+flbuffer)
        
        # get the hull
        fhull = self.hull
        
        # calculate the fire line
        if fhull is None: # if no hull, return None
            return None
        else:  # otherwise, create shape of the active fire line
            if fhull.type == 'MultiPolygon':
                # extract exterior of fire perimeter
                mls = MultiLineString([plg.exterior for plg in fhull])
                # return the part which intersects with  bufferred flinelocsMP
                # return mls.intersection(flinelocsMP.buffer(flbuffer))
                # flinelocsMP_buf = FireVector.addbuffer(flinelocsMP,flbuffer)
                fline = mls.intersection(flinelocsMP)
            
            elif fhull.type == 'Polygon':
                mls = fhull.exterior
                # return mls.intersection(flinelocsMP.buffer(flbuffer))
                # flinelocsMP_buf = FireVector.addbuffer(flinelocsMP,flbuffer)
                fline = mls.intersection(flinelocsMP)
            else:  # if fhull type is not 'MultiPolygon' or 'Polygon', return flinelocsMP
                fline = flinelocsMP
            
            return fline

    @property
    def flinelen(self):
        ''' The length of active fire line
        '''
        try:
            flinelen = self.fline.length/1000 # in km
        except:
            flinelen = 0
        
        return flinelen
    
    @property
    def fline_prior(self):
        '''assembles the fire lines of prior time steps to one multilinestring'''
        from shapely.ops import linemerge
        from shapely.geometry import MultiLineString
        
        flines = [fline for fline in self.fline_prior_list if fline is not None]
        x = []
        for fline in flines:
            if isinstance(fline, MultiLineString):
                x += [line for line in fline]
            else:
                x += [fline]
        if len(x) == 0:
            return None
        else:
            return linemerge(x)
    
    @property
    def spreadrate(self):
        ''' spread rate, computed as shortest distances of all points on the fire line to the prior hull
        output: a list of distances [m]
        '''
        from shapely.geometry import Point
        
        if self.prior_hull:
            # get new fire line pixels
            nps = self.flinepixels
            
            # compute all distances
            dist = [self.prior_hull.distance(Point(p.x,p.y)) for p in nps]
            
            return dist
        else: # this is the ignition time step
            return None
    
    @property
    def spreadavg(self):
        '''Average spread rate
        '''
        spread = self.spreadrate
        if spread is not None:
            if len(spread) > 0:
                spread = sum(spread)/len(spread)
            else:
                spread = 0
            return spread
        else: # this is the ignition time step
            return np.nan
    
    @property
    def spread95(self):
        ''' 95th percentile of spread rate
        '''
        spread = self.spreadrate
        if spread is not None:
            if len(spread) > 0:
                spread = np.percentile(np.array(self.spreadrate), 95)
            else:
                spread = 0
            return spread
        else: # this is the ignition time step
            return np.nan
    
    def updatefhull(self, newlocs):
        """ Update the hull using old hull and new locs
        """
        import FireVector
        from shapely.wkt import loads, dumps
        from shapely.ops import unary_union
        from shapely.geometry import Polygon, MultiPolygon
        
        hull = FireVector.cal_hull(newlocs)
        # use the union to include hull in past time step
        phull = self.hull
        try:
            self.hull = unary_union([phull,hull])
        except:
            if isinstance(hull, Polygon):
                hull = loads(dumps(hull, rounding_precision=2))
            else:
                hull = MultiPolygon([loads(dumps(poly, rounding_precision=2)) for poly in hull])
            #phull = phull.apply(lambda x: shapely.wkt.loads(shapely.wkt.dumps(x, rounding_precision=2)))
            self.hull = unary_union([phull,hull])
    
    def updateextpixels(self, newpixels):
        """ Update the external pixels
        """
        import FireVector
        
        pextpixels = self.extpixels
        self.extpixels = FireVector.cal_extpixels(pextpixels + newpixels, self.hull)
    
    def updateftype(self):
        """ Update fire type
        # do not use ftype as property since it may mess up when t updates (without pixel addition)
        """
        import FireIO
        import random
        year = self.t[0]
        if self.n_newpixels > 0: # only update if there are new pixels
            if (self.n_newpixels < 1000):
                uselocs = self.newlocs_geo
            else:
                # we can do a random sample of 1000 new pixels (it's likely going to be a forest fire anyways)
                uselocs = random.sample(self.newlocs_geo, 1000)
            
            # get all LCT/peat for the fire pixels
            vLCT = FireIO.get_LCT(uselocs,year)
            vPEAT = FireIO.get_peatstatus(uselocs)
            
            # extract the LCT with most pixel counts
            self.LCTmax = max(set(vLCT), key = vLCT.count)
            
            # check if at leat one peat pixel is intersecting the new fire loctions
            self.peat = 1 if 1 in vPEAT else 0
    
    def add_prior_hull(self,src_hull):
        '''Adds a prior_hull of a merged fire 
        for accurate spread rate calculation'''
        
        from shapely.ops import unary_union
        
        tgt_hull = self.prior_hull
        if src_hull is None:
            self.prior_hull = tgt_hull
        elif tgt_hull is None:
            self.prior_hull = src_hull
        else:
            self.prior_hull = unary_union([tgt_hull, src_hull])
        
        

# c. Object - Cluster
class Cluster:
    """ class of active fire pixel cluster at a particular time
    """

    # initilization
    def __init__(self, id, pixels, t):
        ''' initilization

        Parameters
        ----------
        id : int
            cluster id number
        pixels : 3-element list
            (y, x, FRP) of AF pixels
        t : tuple, (int,int,int,str)
            the year, month, day and 'AM'|'PM'
        '''
        from datetime import date
        self.cday = date(*t[:-1])  # current date
        self.ampm = t[-1]          # current ampm
        self.id = id
        self.pixels = pixels       # (x,y, FRP,...)

    # properties
    @property
    def locs(self):
        ''' List of pixel locations (y,x)
        '''
        return [(p.x,p.y) for p in self.pixels]

    # @property
    # def centroid(self):
    #     ''' Centroid of the cluster (y, x)
    #     '''
    #     return FireClustering.cal_centroid(self.locs)

    @property
    def n_pixels(self):
        ''' Number of total pixels
        '''
        return len(self.pixels)

    @property
    def hull(self):
        ''' Fire concave hull (alpha shape)
        '''
        hull = FireVector.cal_hull(self.locs)
        return hull
    
    @property
    def b_box(self):
        ''' Bounding box of concave hull
        '''
        b_box = self.hull.bounds
        return b_box

# d. Object - FirePixel
class FirePixel:
    """ class of an acitve fire pixel, which includes
        loc : location (x, y)
        atts : line & sample or viirs pixel, fire radiative power
        t : time (y,m,d,ampm) of record
        origin : the fire id originally recorded (before merging)
    """
    def __init__(self,x,y,lon,lat,frp,DS,DT,dt,Sat,origin):
        self.x = x
        self.y = y
        self.lon = lon
        self.lat = lat
        self.frp = frp          # frp
        self.DS = DS
        self.DT = DT
        self.sat = Sat          # satellite
        self.datetime = dt      # YYYYMMDD_HHMM
        self.origin = origin    # originate fire id
    
    @property
    def loc(self):
        return (self.x,self.y)
    
    @property
    def loc_geo(self):
        return (self.lon,self.lat)

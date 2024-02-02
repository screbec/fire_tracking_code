""" FireClustering
This module include all functions used for doing fire clustering
Authors: Yang Chen, Rebecca Scholten
"""

def remove_self(inds, dist):
    ''' Remove self from the index and distance arrays

    Parameters
    ----------
    inds : np array of np array
        indices of neighbors for each point (self included)
    dist : np array of np array
        distance of neighbors for each point (self included)

    Returns
    -------
    new_inds : np array of np array
        indices of neighbors (self excluded)
    new_dist : np array of np array
        distance of neighbors (self excluded)
    '''
    import numpy as np

    new_inds = []
    new_dist = []
    for i in range(len(inds)):
       pos = np.where(inds[i] == i)
       new_inds.append(np.delete(inds[i], pos))
       new_dist.append(np.delete(dist[i], pos))

    return np.array(new_inds,dtype=object), np.array(new_dist,dtype=object)

def build_rtree(geoms, fids=[]):
    '''Builds Rtree from a shapely multipolygon shape
    and optionally uses list of fids as identifier'''
    import rtree
    
    idx = rtree.index.Index() # create new index
    for ind, geom in enumerate(geoms):
        if len(fids) > 0:
            idx.insert(ind, geom.bounds, fids[ind])
        else:
            idx.insert(ind, geom.bounds, ind)
    
    return idx

def idx_intersection(idx, bbox):
    ''' 
    Finds all objects in an index whcih bounding boxes intersect with a geometry's bounding box
    '''
    intersections = list(idx.intersection(bbox, objects = True))
    if len(intersections)>0:
        fids, bbox = zip(*[(item.object, item.bbox) for item in intersections])
    else:
        fids = []
    return fids


def compute_all_spatial_distances(data, max_thresh_km):
    ''' Derive neighbors (and distances) for each point (with x/y)

    Parameters
    ----------
    data : pd DataFrame
        point location with 'x' and 'y' columns
    max_thresh_km : float
        maximum distance threshold (km) used for classifying neighbors

    Returns
    -------
    inds : np array of np array
        indices of neighbors (self exluded)
    dist : np array of np array
        distance (km) of neighbors (self excluded)
    '''
    import numpy as np
    from sklearn.neighbors import BallTree

    x = data.x.values
    y = data.y.values
    
    X = np.stack([x, y], axis=-1)

    bt = BallTree(X, leaf_size=20) # we can actually use euclidean now since it is projected

    inds, dist = bt.query_radius(X, r = max_thresh_km*1000, return_distance=True)
    inds, dist = remove_self(inds, dist)
    return inds, dist * 1000

def do_clustering(data, max_thresh_km):
    ''' Do initial clustering for fire pixels

    Parameters
    ----------
    data : list (nx2)
        latitude and longitude values of all fire pixels
    max_thresh_km : float
        maximum distance threshold (km) used for classifying neighbors

    Returns
    -------
    point_to_cluster_id : list
        cluster id for each fire point
    '''
    import numpy as np
    import pandas as pd

    # value to fill in pixels without clustering
    NO_CLUSTER_VAL = -1

    # if number of points is 1 or 2, each point is one cluster
    num_points = len(data)
    if num_points < 3:
        cluster_id = list(range(num_points))
        return cluster_id

    # convert list to pd DataFrame
    dfdata = pd.DataFrame(data,columns=['x','y'])

    # initialization
    cluster_id_counter = 0
    point_to_cluster_id = np.full(num_points, fill_value=NO_CLUSTER_VAL, dtype=np.int64)

    # compute and sort neighbor pixels for each pixel
    neighbor_inds, neighbor_spatial_dists = compute_all_spatial_distances(dfdata, max_thresh_km)
    # neighbor_inds, neighbor_spatial_dists = sort_neighbors(neighbor_inds, neighbor_spatial_dists)

    # include all possible pixels in cluster
    to_check = np.full(num_points, fill_value=1, dtype=np.int8)
    while np.sum(to_check) > 0:
        start_ind = np.argmax(to_check == 1) # catch first index to check
        
        neighbors_to_search = list(neighbor_inds[start_ind])
        all_neighbors = neighbors_to_search
        
        if len(all_neighbors) == 0:  # if no neighbor, record the current pixel as a separate cluster
            point_to_cluster_id[start_ind] = cluster_id_counter
            cluster_id_counter += 1
            to_check[start_ind] = 0
            
        else:  # if with neighbor, record all neighbors
            # find all neighbours of neighbours:
            searched_neighbours = [start_ind]
            while len(neighbors_to_search) > 0:
                # take the first of these
                px = neighbors_to_search[0]
                searched_neighbours.append(px)
                px_neighbors = list(neighbor_inds[px])
                all_neighbors = list(set(all_neighbors + px_neighbors))
                neighbors_to_search = list(set(all_neighbors).difference(searched_neighbours))
            # now we have all pixels in this cluster in all_neighbors
            point_to_cluster_id[all_neighbors] = cluster_id_counter
            cluster_id_counter += 1
            to_check[all_neighbors] = 0
    
    return point_to_cluster_id.tolist()

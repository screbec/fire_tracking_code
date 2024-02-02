a============= Arctic-boreal fire object tracking system using VIIRS active fires =============

[NAME]: Fire object tracking system using VIIRS active fires, optimised for Arctic-boreal regions

[VERSION]: 2.0

[AUTHORS]: Rebecca Scholten, Yang Chen, Casey Graff, Shane Coffield

[LANGUAGE]: Python 3.0 with external packages of
  * Numpy (1.17.5)
  * Pandas (1.0.1)
  * Geopandas (0.7.0)
  * Xarray (0.15.0)
  * Scipy (1.4.1)
  * Shapely (1.7.1)
  * Gdal (3.0.4)
  * rasterio (1.1.0)
  * Pyproj (2.5.0)
  * sklearn

[LIST OF FILES]:
* FireRun.py:        module to control different runs

- Tracking algorithm
* FireObj.py:        module containing the object definitions
* FireMain.py:       main module for running the fire object tracking along time
* FireConsts.py:     constants and options used for the project
* FireIO.py:         functions used to read and save data
* FireTime.py:       functions used to handle times within the tracking algorithm
* FireVector.py:     functions used for vector related calculations
* FireClustering.py: functions used for doing fire clustering
* FireLog.py:        module containing all logging info
* PostProcess.py:    module containing functions for post-processing using surface water data

- Outputs
* FireGdf.py:        module for creating geopackage files of all perimeters at each time
                       step
* FireGdf_ign.py:    module for creating annual geopackage of ignitions
* FireGdf_final:     module for creating annual geopackage of final perimeters
* FireGdf_sfs.py:    module for creating geopackage summary for temporal evolution
                       of each single fire
* FireSummary.py:    module for creating a summary of all valid fires up to a
                       time step (usually at year end)


[COPYRIGHT]: This code is open-source and can be freely used for research purposes

[CONTACTS]: rebecca.scholten@posteo.de

=================================== End ========================================

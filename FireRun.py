""" FireRun
Module to control different runs
"""
def Yearbatchrun(year,tst=None,ted=None,restart=False,final_only=False):
    ''' Run the code for each single year
    '''
    import FireMain,FireGdf,FireGdf_sfs,FireGdf_ign,FireGdf_final,FireSummary
    import time
    from FireConsts import sat

    t1 = time.time()
    # set the start and end time
    if tst is None: tst = (year,1,1,'AM')
    if ted is None: ted = (year,10,31,'PM')
    # if year == 2012: tst = (year,1,20,'AM')
    
    if not final_only:
        # Run the time forward and record daily fire objects .pkl data and fire attributes .GeoJSON data
        FireMain.Fire_Forward(tst=tst,ted=ted,restart=True,sat=sat)
        t2 = time.time()
        print(f'{(t2-t1)/60.} minutes to run algorithm')
    
        # Run to save geojson files for each time step
        FireGdf.save_gdf_trng(tst=tst,ted=ted,fperim=True,fline=True,NFP=True)
        #FireGdf_merge.save_gdf_trng(tst=tst,ted=ted,fperim=True)
        t3 = time.time()
        print(f'{(t3-t2)/60.} minutes to save gpkg files')
    
        # Run to save ignition point layer for each time step
        FireGdf_ign.save_gdf_trng(tst,ted)
    FireGdf_final.save_gdf_trng(tst,ted)
    t31 = time.time()
    print(f'{(t31-t3)/60.} minutes to save ignitions and final perimeters')
    
    # # Run to save large fire geosjon files for each time step
    # FireGdf_sfs.save_gdf_trng(tst=tst,ted=ted,fperim=True)
    # FireGdf_sfs.save_gdf_trng(ted=ted,fperim=True)
    # t4 = time.time()
    # print(f'{(t4-t31)/60.} minutes to save large fires')
    
    

    # # Run to save year end summary and file lists
    FireSummary.save_sum_1d(tst,ted)
    # FireSummary.add_heritage(ted)
    # FireSummary.add_largefirelist(ted)
    # t5 = time.time()
    # print(f'{(t5-t4)/60.} minutes to generate summary')


    t2 = time.time()
    print(f'{(t2-t1)/60.} minutes to run code')

if __name__ == "__main__":
    ''' The main code to run time forwarding for a time period
    '''
    import sys
    sys.path.insert(1, 'D:/fire_atlas/1_code/fireatlas')
    import os
    if 'GDAL_DATA' not in os.environ:
        os.environ['GDAL_DATA'] = r'C:/Users/rebec/anaconda3/envs/py3work/Library/share/gdal' 
        os.environ['PROJ_LIB'] = r'C:/Users/rebec/anaconda3/envs/fireatlas/Library/share/proj' 
    
    for year in range(2023,2024):
        Yearbatchrun(year)
    
    

# Finding the Data Volume Path
'''
Script to create the full, final product data file path from an observation ID. The input required is EITHER the observation ID or the source coordinates.
This then outputs the file path to the mosaick.
'''
## Imports
from astroquery.esasky import ESASky

from astropy.coordinates import SkyCoord
import astropy.units as u

import glob
import os
from tqdm import tqdm

## Functions

## Main Function
def find_filepath(
        source_coords = None, 
        obsid = None, 
        calib_level = 3, 
        mission = 'JWST-NEAR-IR', 
        filt = 'F115W'
        ):

    # Path to the data volume
    data_volume_path = '/data/user/jwst_jwst01_data/jwst-artifacts-repository'

    # Checking that either a Source ID or an Observation ID is provided.    
    if not source_coords and not obsid:
        Exception('You must provide either an observation ID or a Source Coordinate!')
    
    # Checking that the source coordinates are a SkyCoord object.
    if source_coords:
        if not type(source_coords) == SkyCoord:
            Exception('Please provide the source coordinates as a SkyCord object')

    # If no observation ID was provided.
    if obsid == None:
        
        # Query ESA Sky in the specified mission (default JWST NEAR IR) and get a table of observation files.
        query_result = ESASky.query_region_maps(
            position = source_coords,
            radius = 1 * u.arcsecond,
            missions = [mission]
        )
        
        # Convert the obs table into a pandas dataframe and select the wanted filter and level 3 data.
        obsid = (
            query_result[mission]
            .to_pandas()
            .query('calibrationlevel == @calib_level and filter == @filt')
            .observation_id
            .iloc[0]
        )

    # Extract the proposal ID from the found or given observation ID.
    proposal_id = obsid.split('-')[0]

    # Find the specific file name for the specified filter.
    wild_path = f'{data_volume_path}/{proposal_id}/{obsid}*_i2d.fits.gz'
    file_paths = glob.glob(wild_path)
    file_path = ''
    for i in file_paths:
        if filt.lower() in i:
            file_path = i

    if len(file_path) < 0.5:
        raise Exception('Provided Observation ID does not link to any calibration 3 files.')
    
    return file_path
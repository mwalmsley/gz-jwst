"""
Just a refactor of the below 'large download' script from the STScI GitHub repo.
https://github.com/spacetelescope/notebooks/blob/master/notebooks/MAST/Astroquery/large_downloads/companion_script.py
"""
from astroquery.mast import Observations

import os
import glob

import typing



def get_filenames_for_proposal(proposal_id: str) -> List[str]:
    """
    Get the filenames for a given proposal ID.
    """

    volume_path = '/data/user/jwst_jwst01_data/jwst-artifacts-repository'

    # Get the list of observations for this proposal
    df = (
        Observations.query_criteria(
            dataproduct_type = ['image'],
            obs_collection = 'JWST',
            proposal_id = proposal_id,
            calib_level = 3,
            dataRights = 'Public'
        )
        .to_pandas()
    )

    obs_ids = list(df.obs_ids)
    
    folder_path = f'{volume_path}/jw0{proposal_id}/'

    # May remove assertion as will be annoying if it kills pipeline.
    assert os.path.exists(folder_path)

    filenames = []
    
    for i in obs_ids:

        # A lot of variation in the filenames on DLs. Use glob to avoid the unnecessary variation.
        filename_list = glob.glob(f'{folder_path}/{i}*_i2d.fits.gz')

        # File missing, but continuing.
        if len(filename_list) == 0:
            Warning(f'Files missing between MAST and Datalabs: {i}')
            continue 
        
        # This shouldn't happen, tbh. If it does, I'm worried.
        if len(filename_list) > 1:
            Warning(f'More than one file found for observation ID {i}. Skipping for safety...')
            continue
        
        # Append the full filepath to the initialised list.
        filenames.append(filename_list[0])

    return filenames

def download_filenames(filenames, save_dir, max_files=3):
       manifest = Observations.download_products(
               files
               , curl_flag=False
               , productType=['SCIENCE']
               )
    #    TODO use requests or similar to download a few
    # or use curl_flag=True for script version, TBD
    # output should be level 3 products download to appropriate directory
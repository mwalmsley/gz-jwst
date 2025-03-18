# Imports

import pandas as pd
import os
import warnings
from tqdm import tqdm
from astropy.io import fits
from astropy import coordinates as coords
from astropy.wcs import WCS, FITSFixedWarning
from astropy.nddata import Cutout2D
from astropy.table import Table
import numpy as np




def make_single_band_cutout(catalog_loc, image_loc, output_dir = 'output/', cutout_size = (200, 200), id_col = 'ID', ra_col = 'RA', dec_col = 'DEC', prefix = '', suffix = '', fill_value = 0.0):

    # Part 0: make sure things exist
    assert os.path.exists(catalog_loc), "catalog_loc does not exist"
    assert os.path.exists(image_loc), "image_loc does not exist"

    
    # Part I: read catalog fits
    tbl = Table.read(
        input=catalog_loc,
        hdu=1,
        character_as_bytes=False,
        format='fits'
    )
    names = [name for name in tbl.colnames if len(tbl[name].shape) <= 1]
    df_catalog = (tbl[names].to_pandas())    


    # Part II: load downloaded image
    ff = fits.open(image_loc)

    ff_img = ff[1].data
    ff_header = ff[1].header

    with warnings.catch_warnings():
        # Ignore a warning on using DATE-OBS in place of MJD-OBS
        warnings.filterwarnings('ignore', message="'datfix' made the change", category=FITSFixedWarning)
        img_wcs = WCS(ff[1].header)
    
    ff.close()




    # Part III: iterate through target dataframe and make cutouts
    for idx, data in tqdm(df_catalog.iterrows(), total=df_catalog.shape[0]):
        # generate output filename
        cutout_filename = f'{output_dir}{prefix}{int(data[id_col])}{suffix}.fits'
    
        # make cutout
        c = coords.SkyCoord(data[ra_col], data[dec_col], unit='deg')  # defaults to ICRS frame
        try: #TODO do this better
            img_cutout = Cutout2D(ff_img, c, cutout_size, wcs=img_wcs, mode='partial', fill_value=fill_value)
        except:
            continue
    
        # create new FITS
        hdr = ff_header
        empty_primary = fits.PrimaryHDU(header=hdr)
    
        image_hdu1 = fits.ImageHDU(img_cutout.data)
    
        hdul = fits.HDUList([empty_primary, image_hdu1])
    
        # Update the FITS header with the cutout WCS
        hdul[1].header.update(img_cutout.wcs.to_header())
        hdul[1].header['EXTTYPE'] = 'IMAGE' 
    
        # Write the cutout to a new FITS file
        hdul.writeto(cutout_filename, overwrite=True)

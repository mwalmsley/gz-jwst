# Creating the 1D Cutouts

## Imports
# Astropy Stuff
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.visualization import ZScaleInterval, SqrtStretch
from astropy.nddata import Cutout2D
from astropy.wcs import WCS

# Base Packages
import numpy as np
from matplotlib import cm
from PIL import Image
import os

## Functions

## Main Function
def create_cutout(ras, decs, sizes, source_ids, fits_path = None, save_loc = None):

    # Check that the FITS file has been provided.
    if not fits_path:
        Exception('Path to a FITS file is required!')

    # Check that we have a save location, if not create one.
    if not save_loc:
         save_loc = os.getcwd()
         Warning('No save location provided. Saving cutouts in the current working directory.')

    # Open the FITS file.
    with fits.open(fits_path) as hdul:
            header = hdul[1].header
            data = hdul[1].data

    # Create the WCS and a list of SkyCoords.
    w = WCS(header)
    coords = SkyCoord(ra = ras * u.deg, dec = decs * u.deg, frame = 'icrs')

    save_locs = {}
    for source_id, coord, size in zip(source_ids, coords, sizes):

        # Create the cutout.
        cutout = Cutout2D(data, coord, size, wcs = w, mode = 'partial')

        # Calculate the hyper parameters of the image.
        interval = ZScaleInterval(nsamples=5000, contrast=0.001)
        stretch = SqrtStretch()
        im_norm = stretch(interval(cutout.data))

        # Save the image.
        im = (Image.fromarray(np.uint8(cm.Greys_r(im_norm)*255))).convert('RGB')
        im.save(f'{save_loc}/{source_id}.jpeg')
        save_locs[source_id] = f'{save_loc}/{source_id}.jpeg'
    
    return save_locs
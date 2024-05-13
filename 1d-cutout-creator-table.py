# Creating the 1D Cutouts Assuming with a DataFrame / Astropy Table architecture

## Imports
# Archive requirements
from astroquery.esasky import ESASky

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
def create_cutout(ra, dec, size, source_id, fits_path = None, save_loc = None):

    # Checking that a FITS file has been entered.
    if not fits_path:
        Exception('Path to a FITS file is required!')

    if not save_loc:
         save_loc = os.getcwd()
         Warning('No save location provided. Saving cutouts in the current working directory.')

    # Open the FITS file.
    with fits.open(fits_path) as hdul:
        header = hdul[1].header
        data = hdul[1].data

    # Create the WCS and the SkyCoord at the centre of the cutout.
    w = WCS(header)
    coord = SkyCoord(ra = ra * u.deg, dec = dec * u.deg, frame = 'icrs')
    
    cutout = Cutout2D(data, coord, size, wcs = w, mode = 'partial')

    # Create the image parameters.
    interval = ZScaleInterval(nsamples=5000, contrast=0.001)
    stretch = SqrtStretch()
    im_norm = stretch(interval(cutout.data))

    # Image created and saved as a .jpeg.
    im = (Image.fromarray(np.uint8(cm.Greys_r(im_norm)*255))).convert('RGB')
    im.save(f'{save_loc}/{source_id}.jpeg')

    return f'{save_loc}/{source_id}.jpeg'
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
def create_cutout(ra, dec, sizes, source_ids, fits_path = None):

    # Checking that a FITS file has been entered.
    if not fits_path:
        Exception('Path to a FITS file is required!')

    # Open the FITS file.
    with fits.open(fits_path) as hdul:
        header = hdul[1].header
        data = hdul[1].data

    # Create the WCS and the SkyCoord at the centre of the cutout.
    w = WCS(header)
    coords = SkyCoord(ra = ra * u.deg, dec = dec * u.deg, frame = 'icrs')
    
    # Initialise a dictionary which will link the source ID to the output flux values.
    flux_values = {}

    # Loop through all given parameters and create a dictionary of cutouts.
    for source_id, coord, size in zip(source_ids, coords, sizes):

        # Create the cutout.
        cutout = Cutout2D(data, coord, size, wcs = w, mode = 'partial')

        # Calculate the hyper parameters of the image.
        flux_values[source_id] = cutout.data
    
    return flux_values
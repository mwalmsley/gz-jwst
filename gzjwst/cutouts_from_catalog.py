"""
Each cutout function should expect a catalog row
i.e.: a path to a (in general, multiple) fits files, plus coordinate and sizing metadata
Should also expect parameters controlling the cutout colouring, dynamic range, etc
"""

import pandas as pd
import numpy as np
from PIL import Image


from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
import astropy.units as u


def make_single_band_cutout(galaxy: pd.Series, data: np.array, header: fits.header.Header):
    # Check that the FITS file has been provided.
    if not data or not header:
        Exception('Path to a FITS file is required!')

    # Note, may have to update the entries extracted here depending on name in final file.
    ra = galaxy.RA
    dec = galaxy.DEC
    size = galaxy.size

    # Create the WCS and a list of SkyCoords.
    w = WCS(header)
    coord = SkyCoord(ra = ra * u.deg, dec = dec * u.deg, frame = 'icrs')

    # Create the Cutout
    cutout = Cutout2D(data, coord, size, wcs = w, mode = 'partial')

    flux_values = cutout.data

    cutout_resize = np.asarray((Image.fromarray(flux_values)).resize((300,300)))

    return cutout_resize    


# TODO etc for other band combinations, copy from Hayley's notebook


def select_2d_cutout():
    pass # TODO David or Hayley's slice code to select 2D flux values from a mosaic


"""
Each cutout function should expect a catalog row
i.e.: a path to a (in general, multiple) fits files, plus coordinate and sizing metadata
Should also expect parameters controlling the cutout colouring, dynamic range, etc
"""

import pandas as pd
import numpy as np

from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
import astropy.units as u


def make_single_band_cutout(galaxy: pd.Series, data: np.array, header: fits.header.Header):
    # Note, may have to update the entries extracted here depending on name in final file.
    # Assumed that ra, dec and size are in deg, deg and arcseconds respectively. Change if thought!
    ra = galaxy.RA
    dec = galaxy.DEC
    size = galaxy.size

    # Create the WCS and a list of SkyCoords.
    w = WCS(header)
    coord = SkyCoord(ra = ra * u.deg, dec = dec * u.deg, frame = 'icrs')

    # Create the Cutout
    cutout = Cutout2D(data, coord, size * u.arcsec, wcs = w, mode = 'partial')

    flux_values = cutout.data

    return flux_values 


# TODO etc for other band combinations, copy from Hayley's notebook


def select_2d_cutout():
    pass # TODO David or Hayley's slice code to select 2D flux values from a mosaic


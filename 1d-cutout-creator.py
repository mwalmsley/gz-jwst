# Creating the 1D Cutouts

## Imports
# Astropy Stuff
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.wcs import WCS

## Functions

## Main Function
def create_cutout(ras, decs, size, source_id, fits_path = None):

    # Check that the FITS file has been provided.
    if not fits_path:
        Exception('Path to a FITS file is required!')

    # Open the FITS file.
    with fits.open(fits_path) as hdul:
            header = hdul[1].header
            data = hdul[1].data

    # Create the WCS and a list of SkyCoords.
    w = WCS(header)
    coord = SkyCoord(ra = ras * u.deg, dec = decs * u.deg, frame = 'icrs')

    # Create the Cutout
    cutout = Cutout2D(data, coord, size, wcs = w, mode = 'partial')

    flux_values = cutout.data

    return flux_values


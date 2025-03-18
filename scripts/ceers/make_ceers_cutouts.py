import pandas as pd
import os
import warnings
from tqdm import tqdm
from astropy.io import fits
from astropy import coordinates as coords
from astropy.wcs import WCS, FITSFixedWarning
from astropy.nddata import Cutout2D
from astropy.table import Table


# read catalogue
tbl = Table.read(
    input='CEERS_optap.fits',
    hdu=1,
    character_as_bytes=False,
    format='fits'
)
names = [name for name in tbl.colnames if len(tbl[name].shape) <= 1]
df_ceers = (tbl[names].to_pandas())

# load downloaded Mosaic
fits_file = 'hlsp_ceers_jwst_nircam_nircam2all_f200w_v0.5_i2d.fits.gz'
if os.path.exists(fits_file):
    ff = fits.open(fits_file)

    ff_img = ff[1].data
    ff_header = ff[1].header

    with warnings.catch_warnings():
        # Ignore a warning on using DATE-OBS in place of MJD-OBS
        warnings.filterwarnings('ignore', message="'datfix' made the change", category=FITSFixedWarning)
        img_wcs = WCS(ff[1].header)
    
    ff.close()


# select some objects
df_targets = df_ceers[(df_ceers['RA']>=214.90) & (df_ceers['RA']<=214.91) & (df_ceers['DEC']>=52.90) & (df_ceers['DEC']<=52.91)]

# specify output directory
CUTOUT_FILE_PATH = '../../data/ceers/cutouts/'

# specify size of cutout
cutout_size = (200, 200)

# iterate through target dataframe and make cutouts
for idx, data in tqdm(df_targets.iterrows(), total=df_targets.shape[0]):
    # generate output filename
    cutout_filename = CUTOUT_FILE_PATH + str(int(data['ID'])) + '_CEERS.fits'

    # make cutout
    c = coords.SkyCoord(data['RA'], data['DEC'], unit='deg')  # defaults to ICRS frame
    img_cutout = Cutout2D(ff_img, c, cutout_size, wcs=img_wcs, mode='partial', fill_value=0.0)

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



from gzjwst import catalog_for_survey
import glob

# dev ceers data 
filter_band_kwarg = 'f200w'
list_of_relevant_files = glob.glob('data/dev_data/hlsp_ceers*{}*i2d.fits.gz'.format(filter_band_kwarg))
assert list_of_relevant_files, f"No files found"
filename = list_of_relevant_files[0]
print(filename)
# exit()

# all ceers images
# folder_path = 'data/ceers/hlsp'  # relative path should work for everyone
# filter_band_kwarg = 'f200w'
# list_of_relevant_files = glob.glob(f'{folder_path}/*_{filter_band_kwarg}_*/*.fits.gz')




sep_config = catalog_for_survey.SEPConfig() # defaults
# override here if needed

processed_table_given_bkg, cleaned_segmap_given_bkg = catalog_for_survey.extract_objects_from_mosaic(
    filename, use_existing_bkg=True, fits_ext=2, chop_mosaic=False, sep_config=sep_config
)
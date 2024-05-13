"""
Contents of source_extraction.py should go in here
"""

def get_catalog(mosaic_dir, detection_kwargs):
    # TODO load mosaics from that directory
    objects, segmap = identify_objects(image_data, **detection_kwargs)
    # TODO save segmap to survey/catalogs/segmap.fits or similar
    # TODO convert objects to pandas table and save to survey/catalogs/catalog.csv or similar
    # catalog should include relative paths to the original mosaics

def identify_objects(image_data,nsigma,min_area,deb_n_thresh,deb_cont,filter_kwarg,filter_size):
    pass

# etc from source_extraction.py
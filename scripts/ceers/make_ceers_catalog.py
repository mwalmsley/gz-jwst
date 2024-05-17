from gzjwst import catalog_for_survey as se # have to fix relative importing
import glob


folder_path = '/d1/manth145/jwst/ceers/HLSP'
filter_band_kwarg = 'f200w'
sep_kwarg_dict={'nsigma':3, 
                'min_area': 10, 
                'deb_n_thresh': 32, 
                'deb_cont': 0.0001, 
                'filter_kwarg': 'tophat', 
                'filter_size': 7}

list_of_relevant_files = glob.glob(f'{folder_path}/*_{filter_band_kwarg}_*/*.fits.gz')

processed_table_given_bkg, cleaned_segmap_given_bkg = se.extract_objects_from_mosaic(list_of_relevant_files[0], use_existing_bkg=True, fits_ext=2, chop_mosaic=False, sep_kwarg_dict=sep_kwarg_dict)
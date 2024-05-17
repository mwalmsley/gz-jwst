"""
Contents of source_extraction.py should go in here
"""
import sep
import numpy as np
from astropy.io import fits, ascii
import sep # the source extraction module.
from astropy.convolution import Tophat2DKernel, Gaussian2DKernel, Box2DKernel #kernels used for convolution during source extraction
from astropy.stats import gaussian_fwhm_to_sigma
import matplotlib.pyplot as plt
from copy import deepcopy
from scipy import ndimage
import warnings
from astropy.wcs import WCS
from astropy.table import Table
import os

# sep.set_extract_pixstack(1000000)

# Notes from team below:
# def get_catalog(mosaic_dir, detection_kwargs):
#     # TODO load mosaics from that directory
#     objects, segmap = identify_objects(image_data, **detection_kwargs)
#     # TODO save segmap to survey/catalogs/segmap.fits or similar
#     # TODO convert objects to pandas table and save to survey/catalogs/catalog.csv or similar
#     # catalog should include relative paths to the original mosaics

# def identify_objects(image_data,nsigma,min_area,deb_n_thresh,deb_cont,filter_kwarg,filter_size):
#     pass

# etc from source_extraction.py

def identify_objects(image_data, bkg_data=None, params={}):
    '''
    This function performs source identification using the python-based module named SEP (Barbary et al., 2016).
    :param image_data: provide the image data as an numpy.ndarray.
    :param nsigma: significance above the computed background where sources are identified.
    :param min_area: minimum area of the contiguous regions to be identified as a source.
    :param deb_n_thresh: number of thresholds to be applied during deblending sources.
    :param deb_cont: deblend contrast ratio
    :param param_dict: parameter dictionary which contains user specified information about source extraction parameters.
    :return: object list identified in the image, segmentation map with each object labeled from 1, 2, 3 ... The catalog of objects is
    ordered as per the segmap number.
    '''
    nsigma = params['nsigma']
    min_area = params['min_area']
    deb_n_thresh = params['deb_n_thresh']
    deb_cont = params['deb_cont']
    filter_kwarg = params['filter_kwarg']
    filter_size = params['filter_size']
    
    byte_swaped_data = image_data.byteswap().newbyteorder() # as suggested by the SEP documentation.
    
    if bkg_data is None:
        global_bkg = sep.Background(byte_swaped_data) #calculating the global background
        glob_rms = global_bkg.globalrms
        print(glob_rms)
    else:
        bkg, bkg_mask, rms_data = bkg_data
        global_bkg = deepcopy(bkg.byteswap().newbyteorder())
        bkg_mask[np.where(bkg_mask>1)] = 1
        bkg_mask+=1
        bkg_mask[np.where(bkg_mask>1)] = 0
        
        masked_rms_data = rms_data * bkg_mask
        glob_rms = np.nanpercentile(masked_rms_data[masked_rms_data!=0], 50)
        print(f'The global median RMS = {glob_rms}')

    bkg_subtracted = byte_swaped_data - global_bkg #background subtracted data

    if filter_kwarg.lower() not in ['tophat','gauss','boxcar']:#check if the filter keyword exists in list of supported filters.
        warnings.warn('The filter %s is not supported as of yet, defaulting to tophat of radius 5') #if doesn't exist, then warn the user
        source_kernel = Tophat2DKernel(5) #default to the Tophat kernel of radius 5 pixels
    elif filter_kwarg.lower() == 'tophat':
        source_kernel = Tophat2DKernel(filter_size)
    elif filter_kwarg.lower() == 'gauss':#for gaussian filter, the size provided should be a fwhm and it will be converted to a sigma internally.
        _gauss_sigma = gaussian_fwhm_to_sigma * filter_size
        source_kernel = Gaussian2DKernel(_gauss_sigma)#initiating a gaussian kernel
    elif filter_kwarg.lower() == 'boxcar':#boxcar kernel
        source_kernel = Box2DKernel(filter_size)# initiating a boxcar kernel of side length = filter_size

    ## Extract the objects and generate a segmentation map.
    objects, segmap = sep.extract(bkg_subtracted, nsigma, err=glob_rms, minarea=min_area, deblend_nthresh=deb_n_thresh, deblend_cont=deb_cont,segmentation_map=True, filter_kernel=source_kernel.array)

    return objects, segmap # return the objects and segmentation map.


def extract_objects_from_mosaic(filename, fits_ext, use_existing_bkg=False, sep_kwarg_dict={'nsigma':3, 'min_area': 10, 'deb_n_thresh': 32, 'deb_cont': 0.0001, 'filter_kwarg': 'tophat', 'filter_size': 7}, chop_mosaic=False, chop_loc=4500, write_info=True, overwrite_existing_catalogs=True):
    image_data, image_header = fits.getdata(filename, memmap=True, ext=fits_ext, header=True)
    image_wcs = WCS(image_header)
    if use_existing_bkg:
        background_info = (fits.getdata(filename, memmap=True, ext=9), fits.getdata(filename, memmap=True, ext=10), fits.getdata(filename, memmap=True, ext=3))
    else:
        background_info = None
    
    if chop_mosaic:
        mosaic_part1, mosaic_part2 = chop_larger_mosaic(image_data, chop_loc=chop_loc, return_which_axis='y')
        bkg_part1, bkg_part2 = chop_larger_mosaic(background_info[0], chop_loc=chop_loc, return_which_axis='y')
        bkg_mask_part1, bkg_mask_part2 = chop_larger_mosaic(background_info[1], chop_loc=chop_loc, return_which_axis='y')

        for mosaic, bk_info in zip([mosaic_part1], [(bkg_part1, bkg_mask_part1)]):
            objects, segmap = identify_objects(mosaic, bkg_data= bk_info, params=sep_kwarg_dict)

            objects_table = Table(objects, names=list(objects.dtype.names))
            objects_table['object_id'] = np.arange(len(objects_table))+1
            chip_gap_mask = make_chip_gap_mask(mosaic)

            dilated_chip_gap_mask = dilate_mask(chip_gap_mask, number_of_dilations=25)
            cleaned_mask = clean_artefacts_from_segmap(segmap,dilated_chip_gap_mask)
            processed_table = process_object_table(objects_table, segmap, dilated_chip_gap_mask)
    else:
        mosaic = image_data
        objects, segmap = identify_objects(mosaic, bkg_data= background_info, params=sep_kwarg_dict)

        objects_table = Table(objects, names=list(objects.dtype.names))
        print(f'Source extraction found {len(objects_table)} Objects!')
        
        objects_table['object_id'] = np.arange(len(objects_table))+1
        chip_gap_mask = make_chip_gap_mask(mosaic)

        dilated_chip_gap_mask = dilate_mask(chip_gap_mask, number_of_dilations=25)
        cleaned_mask = clean_artefacts_from_segmap(segmap,dilated_chip_gap_mask)
        processed_table = process_object_table(objects_table, segmap, dilated_chip_gap_mask)
        ra, dec = image_wcs.wcs_pix2world(np.array(processed_table['x']), np.array(processed_table['y']), 0)
        processed_table['ra'] = ra
        processed_table['dec'] = dec
        print(f'Processed the detected objects... now has {len(processed_table)} objects!')
        
        if write_info:
            save_path = os.path.dirname(filename)
            save_file_basename = os.path.basename(filename).strip('.fits.gz')
            ascii.write(processed_table, f'{save_path}/{save_file_basename}_catalog.csv', format='csv', overwrite=overwrite_existing_catalogs)
            fits.writeto(f'{save_path}/{save_file_basename}_segmap.fits', cleaned_mask, image_header, overwrite=overwrite_existing_catalogs)
            print(f'Written the source catalog and segmentation map at {save_path}')
                    
    return processed_table, cleaned_mask


def chop_larger_mosaic(mosaic_data, chop_loc, return_which_axis='y'):
    if return_which_axis == 'x':
        return mosaic_data[:chop_loc,:], mosaic_data[chop_loc:,:]
    elif return_which_axis == 'y':
        return mosaic_data[:,:chop_loc], mosaic_data[:,chop_loc:]
    else:
        assert return_which_axis in ['x', 'y']



def make_chip_gap_mask(raw_image_data):
    chip_gap_mask = np.zeros_like(raw_image_data)
    chip_gap_mask[np.where(raw_image_data == 0)] = 1
    return chip_gap_mask

def dilate_mask(initial_mask, number_of_dilations=25):
    dilation_structure = ndimage.generate_binary_structure(2, 2)
    dilated_mask = ndimage.binary_dilation(initial_mask, structure=dilation_structure, iterations=number_of_dilations, border_value=1).astype(int)
    return dilated_mask


def clean_artefacts_from_segmap(original_segmap, chip_mask):
    processed_segmask = deepcopy(original_segmap)

    for each_artefact in np.unique(original_segmap*chip_mask):
        processed_segmask[np.where(processed_segmask == each_artefact)] = 0
    return processed_segmask


def process_object_table(raw_object_table, original_segmap, chip_mask):
    copied_table = deepcopy(raw_object_table)
    copied_table.remove_rows(np.unique(original_segmap*chip_mask)[1:].astype(int)-1)
    return copied_table

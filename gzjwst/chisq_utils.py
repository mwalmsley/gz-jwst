# use chi-squared stacked detection image and photutils to create JWST catalog of sources
import logging
import tqdm
import numpy as np
from scipy.special import factorial
from astropy.io import fits
from scipy.stats import chi2


def get_detection_image(fits_locs: list) -> np.ndarray:
    
    normalised_science_images = []
    dof_counter = []
    
    for loc in tqdm.tqdm(fits_locs):
      logging.info(f'loading {loc}')
      science_image = fits.getdata(loc, memmap=False, ext=1, header=False)  # background-subtracted science image

      if isinstance(science_image, tuple):
        science_image = science_image[0]
      background_source_mask = fits.getdata(loc, memmap=False, ext=10, header=False)
      assert science_image.shape == background_source_mask.shape, f"science image and mask shape mismatch"
    
      random_image = get_alleged_random_background_image(science_image, background_source_mask)
      loc = np.mean(random_image) 
      assert np.abs(loc) < 1e-4, f"mean of random image is {loc}"
      # don't use loc further

      scale = std_estimator(random_image)
      logging.info(f'scale is {scale}')
      # assumes the images are already lined up in WCS and pixel scale
      normalised_science_images.append((science_image / scale).astype(np.float32))
      dof_counter.append((~np.isnan(science_image) | (science_image != 0)).astype(int))

      del science_image
      del background_source_mask

    chisq_image = np.sum(np.stack(normalised_science_images, axis=-1), axis=-1)
    dof = np.sum(np.stack(dof_counter, axis=-1), axis=-1)
    
    pixel_prob_of_source = get_pixel_prob_of_source(chisq_image, dof)
    pixel_prob_of_source = np.where(dof == 0, 0, pixel_prob_of_source)
    return pixel_prob_of_source, dof

    # using the R image convention from szalay et al 1999, 3.2?
    # return np.sqrt(np.sum(np.stack(normalised_science_images, axis=-1), axis=-1))

def get_alleged_random_background_image(science_image_bkg_subtracted, background_source_mask):
    alleged_random_image = science_image_bkg_subtracted[~background_source_mask.astype(bool)]
    # handful of values set to crazy values
    alleged_random_image = alleged_random_image[alleged_random_image < np.percentile(alleged_random_image, 99.9)]
    alleged_random_image = alleged_random_image[alleged_random_image > np.percentile(alleged_random_image, 0.01)]
    return alleged_random_image

def std_estimator(x):
    # return np.sqrt(np.sum((x - np.mean(x))**2) / (len(x) - 1))
    # no, general biased esimator
    # https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation#Rule_of_thumb_for_the_normal_distribution
    return np.sqrt(np.sum((x - np.mean(x))**2) / (len(x)-1.5))  # biased estimator


def get_pixel_prob_of_source(y, dof_counter):
   
  #  szalay eqn 3
  # return 1/(2**(dof_counter/2) * factorial(dof_counter/2)-1) * np.exp(-y/2) * y**(dof_counter/2-1)
  return chi2.cdf(y, dof_counter)  # 1 - cdf

    # return 0.5 * r_image**3 * np.exp(-r_image**2 / 2)
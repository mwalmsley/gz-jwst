#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 10:14:43 2024

@author: husmak
"""


import os
import sep
import sys
import fitsio
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from fitsio import FITS,FITSHDR
from matplotlib import rcParams
from astropy.nddata import Cutout2D
from matplotlib.patches import Ellipse
from astropy.coordinates import SkyCoord
from astropy.visualization import make_lupton_rgb

'''
all are 5-sigma detections

'''



'''
Reads in the fits files of different image bands from JWST NIRCam.

They correspond to the red, green, and blue channels for rgb image respectively. f115w and f150w are
combined to make blue. f277w is green, and f444w is red. Here, make_lupton_rgb is used to combine and 
make into colour images.

i_name is f150w
i_name_2 is f115w
g_name is f227w
r_name is f444w

'''

i_name = '/Users/husmak/Downloads/nircam1/hlsp_ceers_jwst_nircam_nircam1_f150w_dr0.5_i2d.fits'
i_name_2 = '/Users/husmak/Downloads/nircam1/hlsp_ceers_jwst_nircam_nircam1_f115w_dr0.5_i2d.fits'
g_name = '/Users/husmak/Downloads/nircam1/hlsp_ceers_jwst_nircam_nircam1_f277w_dr0.5_i2d.fits'
r_name = '/Users/husmak/Downloads/nircam1/hlsp_ceers_jwst_nircam_nircam1_f444w_dr0.5_i2d.fits'




'''
Data is being taken from each of the channels respectively.
*blue* is the f115w and f150w bands combined.
'''

i1 = fits.open(i_name)[1].data
i2 = fits.open(i_name_2)[1].data
blue = i1 + i2
g = fits.open(g_name)[1].data
r = fits.open(r_name)[1].data




'''
Figure and axis are made, then make_lupton_rgb is used to create a colour image using the
colour channels, r, g, and blue. The colour image is assigned to the variable 'rgb'.

Full rbg image is plotted.

variable *minimum* adjusts the noise floor.

*filename* saves the image in the directory of this python file.
'''
plt.figure()
fig, ax = plt.subplots(figsize=(8, 8))
rgb = make_lupton_rgb(r*1.1, g*0.75, blue*0.4,  minimum=-0.025, Q=10, stretch=0.5, filename="CEERS 1 Colour.png") #add noise
im = ax.imshow(rgb, interpolation = 'nearest', origin = 'lower')






'''
gets r, g, and b channels individually from the rgb image.
'''

r_fits, g_fits, b_fits = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
#total_fits = r_fits + g_fits + b_fits
#total_objects = sep.extract(total_fits, 100)



'''
creates a .fits file from the 'total_fits' image which is all the channels laid on each other.
Places the images on the desktop.
'''

rgb_fits = r_fits + g_fits + b_fits

# Create a simple 2x2 array as data
total_data = rgb_fits

# Create a Primary HDU (Header Data Unit)
hdu = fits.PrimaryHDU(total_data)

# Create an HDU list
hdul = fits.HDUList([hdu])

# Define the path to your desktop
desktop_path = os.path.join(os.path.expanduser('~'), 'Desktop')

# Define the filename for your FITS file
fits_filename = 'RGB_CEERS1_FITS_IMAGE.fits'

# Combine the desktop path and filename to get the full path
fits_filepath = os.path.join(desktop_path, fits_filename)

# Write the FITS file to the specified path
hdul.writeto(fits_filepath, overwrite=True)

print(f"FITS file saved to: {fits_filepath}")





'''
opens csv files of the sources that have been detected from the full size COSMOS image with a 5-sigma
detection threshold. open_matched_sources is all the sources that are in ALL 4 bands
'''

open_f115w_sources = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/match1_only_CEERS1.csv") # 
open_f150w_sources = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/match2_only_CEERS1.csv") #
open_f277w_sources = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/match3_only_CEERS1.csv") # 
open_f444w_sources = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/match4_only_CEERS1.csv") # 
open_matched_sources = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/ALL_CEERS1_SOURCES.csv") # 




#%%

'''
full sized colour image with the raised noise floor.
'''
plt.figure()
fig, ax = plt.subplots(figsize=(10, 10)) #make bigger
rgb_default = make_lupton_rgb(r*1.1, g*0.75, blue*0.4, minimum=-0.025, Q=10, stretch=0.5, filename="CEERS 1 colour (raised noise floor).png")
im = ax.imshow(rgb_default, interpolation = 'nearest', origin = 'lower')





'''
plots ellipses onto full rgb colour image. Each band corresponds to a specific colour of ellipse.
So we know which source was detected in each band.

... expand

'''
for i in range(len(open_f115w_sources)):
    e = Ellipse(xy=(open_f115w_sources['x_f115w'][i], open_f115w_sources['y_f115w'][i]),
                width = 6*open_f115w_sources['a_f115w'][i],
                height = 6*open_f115w_sources['b_f115w'][i],
                angle = open_f115w_sources['theta_f115w'][i]* 180. / np.pi, lw=1)
    e.set_facecolor('none')
    e.set_edgecolor('teal')
    ax.add_artist(e)
    
for i in range(len(open_f150w_sources)):
    e = Ellipse(xy=(open_f150w_sources['x_f150w'][i], open_f150w_sources['y_f150w'][i]),
                width = 6*open_f150w_sources['a_f150w'][i],
                height = 6*open_f150w_sources['b_f150w'][i],
                angle = open_f150w_sources['theta_f150w'][i]* 180. / np.pi, lw=1)

    e.set_facecolor('none')
    e.set_edgecolor('blue')
    ax.add_artist(e)
    

for i in range(len(open_f277w_sources)):
    e = Ellipse(xy=(open_f277w_sources['x_f277w'][i], open_f277w_sources['y_f277w'][i]),
                width = 6*open_f277w_sources['a_f277w'][i],
                height = 6*open_f277w_sources['b_f277w'][i],
                angle = open_f277w_sources['theta_f277w'][i]* 180. / np.pi, lw=1)

    e.set_facecolor('none')
    e.set_edgecolor('green')
    ax.add_artist(e)



for i in range(len(open_f444w_sources)):
    e = Ellipse(xy=(open_f444w_sources['x_f444w'][i], open_f444w_sources['y_f444w'][i]),
                width = 6*open_f444w_sources['a_f444w'][i],
                height = 6*open_f444w_sources['b_f444w'][i],
                angle = open_f444w_sources['theta_f444w'][i]* 180. / np.pi, lw=1)

    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)
    

ax.set_axis_off()
fig.savefig("CEERS 1 RGB Image with ellipses.png", bbox_inches='tight', pad_inches = 0, dpi=2384.2)



#%%

'''
This plots the sources from the CEERS 1 pointing that were found in the f115w band that
are 32 px in size at least.

The sources are made into rgb colour images.
'''

'''
This section is for all sources detected in the f115w band that are 32px in size at least

make a max and min
'''


plt.figure()
fig, ax = plt.subplots(figsize=(10, 10)) #make bigger


f115w_list = []

eligible_f115w = np.array(np.where(open_matched_sources['npix_f115w'] >= 32))



for i in range(len(eligible_f115w[0])): # run for the length of sources being looked at

    if open_matched_sources['npix_f115w'][i] >= 32:
        
        width_f115w = (open_matched_sources['xmax_f115w'][i] - open_matched_sources['xmin_f115w'][i]) * 1.5
        width_f150w = (open_matched_sources['xmax_f150w'][i] - open_matched_sources['xmin_f150w'][i]) * 1.5
        width_f277w = (open_matched_sources['xmax_f277w'][i] - open_matched_sources['xmin_f277w'][i]) * 1.5
        width_f444w = (open_matched_sources['xmax_f444w'][i] - open_matched_sources['xmin_f444w'][i]) * 1.5

        height_f115w = (open_matched_sources['ymax_f115w'][i] - open_matched_sources['ymin_f115w'][i]) * 1.5
        height_f150w = (open_matched_sources['ymax_f150w'][i] - open_matched_sources['ymin_f150w'][i]) * 1.5
        height_f277w = (open_matched_sources['ymax_f277w'][i] - open_matched_sources['ymin_f277w'][i]) * 1.5
        height_f444w = (open_matched_sources['ymax_f444w'][i] - open_matched_sources['ymin_f444w'][i]) * 1.5
        
        width = max(width_f115w, width_f150w, width_f277w, width_f444w)
        height = max(height_f115w, height_f150w, height_f277w, height_f444w)
        

        cutout_width = np.sqrt(width**2 + height**2)
   
        f115w_list.append("image " + str(i) + ". no. of pixels: " + str(open_matched_sources['npix_f115w'][i]))
        
        small_cutout_r = Cutout2D(r, (open_matched_sources['x_f115w'][i], open_matched_sources['y_f115w'][i]), size=(cutout_width,cutout_width))
        small_cutout_g = Cutout2D(g, (open_matched_sources['x_f115w'][i], open_matched_sources['y_f115w'][i]), size=(cutout_width,cutout_width))
        small_cutout_b = Cutout2D(blue, (open_matched_sources['x_f115w'][i], open_matched_sources['y_f115w'][i]), size=(cutout_width,cutout_width))
        
        small_r_name = small_cutout_r.data
        small_g_name = small_cutout_g.data
        small_i_name = small_cutout_b.data
        
        small_r_name = small_r_name.astype('float64')
        small_g_name = small_g_name.astype('float64')
        small_i_name = small_i_name.astype('float64')
        
        
        plt.figure()
        fig, ax = plt.subplots(figsize=(8, 8)) #make bigger
        small_rgb_default = make_lupton_rgb(small_r_name*1.1, small_g_name*0.75, small_i_name*0.4, minimum=-0.025, Q=10, stretch=0.5, filename=("CEERS1 RGB Cutout f115w" + str(i) + ".png"  ))
        im = ax.imshow(small_rgb_default,  interpolation = 'nearest', origin = 'lower')
        



#%%
'''
This plots the sources from the CEERS 1 pointing that were found in the f150w band that
are 32 px in size at least, and aren't duplicates of those seen in the f115w band, 
according to their x and y coordinates.


should maybe avg height or pick the largest of the bands??
'''


plt.figure()
fig, ax = plt.subplots(figsize=(10, 10)) #make bigger


f150w_list = []

eligible_f150w = np.array(np.where(open_matched_sources['npix_f150w'] >= 32))


for i in range(len(eligible_f150w[0])): # run for the length of sources being looked at

    if open_matched_sources['npix_f150w'][i] >= 32:
        
        width_f115w = (open_matched_sources['xmax_f115w'][i] - open_matched_sources['xmin_f115w'][i]) * 1.5
        width_f150w = (open_matched_sources['xmax_f150w'][i] - open_matched_sources['xmin_f150w'][i]) * 1.5
        width_f277w = (open_matched_sources['xmax_f277w'][i] - open_matched_sources['xmin_f277w'][i]) * 1.5
        width_f444w = (open_matched_sources['xmax_f444w'][i] - open_matched_sources['xmin_f444w'][i]) * 1.5

        height_f115w = (open_matched_sources['ymax_f115w'][i] - open_matched_sources['ymin_f115w'][i]) * 1.5
        height_f150w = (open_matched_sources['ymax_f150w'][i] - open_matched_sources['ymin_f150w'][i]) * 1.5
        height_f277w = (open_matched_sources['ymax_f277w'][i] - open_matched_sources['ymin_f277w'][i]) * 1.5
        height_f444w = (open_matched_sources['ymax_f444w'][i] - open_matched_sources['ymin_f444w'][i]) * 1.5
        
        width = max(width_f115w, width_f150w, width_f277w, width_f444w)
        height = max(height_f115w, height_f150w, height_f277w, height_f444w)
        

        cutout_width = np.sqrt(width**2 + height**2)

                
        f150w_list.append("image " + str(i) + ". no. of pixels: " + str(open_matched_sources['npix_f150w'][i]))
        
        small_cutout_r = Cutout2D(r, (open_matched_sources['x_f150w'][i], open_matched_sources['y_f150w'][i]), size=(cutout_width,cutout_width))
        small_cutout_g = Cutout2D(g, (open_matched_sources['x_f150w'][i], open_matched_sources['y_f150w'][i]), size=(cutout_width,cutout_width))
        small_cutout_b = Cutout2D(blue, (open_matched_sources['x_f150w'][i], open_matched_sources['y_f150w'][i]), size=(cutout_width,cutout_width))
        
        small_r_name = small_cutout_r.data
        small_g_name = small_cutout_g.data
        small_i_name = small_cutout_b.data
        
        small_r_name = small_r_name.astype('float64')
        small_g_name = small_g_name.astype('float64')
        small_i_name = small_i_name.astype('float64')
        
        
        plt.figure()
        fig, ax = plt.subplots(figsize=(8, 8)) #make bigger
        small_rgb_default = make_lupton_rgb(small_r_name*1.1, small_g_name*0.75, small_i_name*0.4, minimum=-0.025, Q=10, stretch=0.5, filename=("CEERS1 RGB Cutout f150w" + str(i) + ".png"  ))
        im = ax.imshow(small_rgb_default,  interpolation = 'nearest', origin = 'lower')
        




 
#%%
'''
This plots the sources from the CEERS 1 pointing that were found in the f277w band that
are 32 px in size at least, and aren't duplicates of those seen in the f115w and f150w
bands, according to their x and y coordinates.
'''


plt.figure()
fig, ax = plt.subplots(figsize=(10, 10)) #make bigger


f277w_list = []

eligible_f277w = np.array(np.where(open_matched_sources['npix_f277w'] >= 32))


for i in range(len(eligible_f277w[0])): # run for the length of sources being looked at

    if open_matched_sources['npix_f277w'][i] >= 32:
        
        width_f115w = (open_matched_sources['xmax_f115w'][i] - open_matched_sources['xmin_f115w'][i]) * 1.5
        width_f150w = (open_matched_sources['xmax_f150w'][i] - open_matched_sources['xmin_f150w'][i]) * 1.5
        width_f277w = (open_matched_sources['xmax_f277w'][i] - open_matched_sources['xmin_f277w'][i]) * 1.5
        width_f444w = (open_matched_sources['xmax_f444w'][i] - open_matched_sources['xmin_f444w'][i]) * 1.5

        height_f115w = (open_matched_sources['ymax_f115w'][i] - open_matched_sources['ymin_f115w'][i]) * 1.5
        height_f150w = (open_matched_sources['ymax_f150w'][i] - open_matched_sources['ymin_f150w'][i]) * 1.5
        height_f277w = (open_matched_sources['ymax_f277w'][i] - open_matched_sources['ymin_f277w'][i]) * 1.5
        height_f444w = (open_matched_sources['ymax_f444w'][i] - open_matched_sources['ymin_f444w'][i]) * 1.5
        
        width = max(width_f115w, width_f150w, width_f277w, width_f444w)
        height = max(height_f115w, height_f150w, height_f277w, height_f444w)
        

        cutout_width = np.sqrt(width**2 + height**2)


                
        f277w_list.append("image " + str(i) + ". no. of pixels: " + str(open_matched_sources['npix_f277w'][i]))
        
        small_cutout_r = Cutout2D(r, (open_matched_sources['x_f277w'][i], open_matched_sources['y_f277w'][i]), size=(cutout_width,cutout_width))
        small_cutout_g = Cutout2D(g, (open_matched_sources['x_f277w'][i], open_matched_sources['y_f277w'][i]), size=(cutout_width,cutout_width))
        small_cutout_b = Cutout2D(blue, (open_matched_sources['x_f277w'][i], open_matched_sources['y_f277w'][i]), size=(cutout_width,cutout_width))
        
        small_r_name = small_cutout_r.data
        small_g_name = small_cutout_g.data
        small_i_name = small_cutout_b.data
        
        small_r_name = small_r_name.astype('float64')
        small_g_name = small_g_name.astype('float64')
        small_i_name = small_i_name.astype('float64')
        
        
        plt.figure()
        fig, ax = plt.subplots(figsize=(8, 8)) #make bigger
        small_rgb_default = make_lupton_rgb(small_r_name*1.1, small_g_name*0.75, small_i_name*0.4, minimum=-0.025, Q=10, stretch=0.5, filename=("CEERS1 RGB Cutout f277w" + str(i) + ".png"  ))
        im = ax.imshow(small_rgb_default,  interpolation = 'nearest', origin = 'lower')
        





#%%
'''
This plots the sources from the CEERS 1 pointing that were found in the f444w band that
are 32 px in size at least, and aren't duplicates of those seen in the f115w, f277w,
and f150w bands, according to their x and y coordinates.
'''


plt.figure()
fig, ax = plt.subplots(figsize=(10, 10)) #make bigger


f444w_list = []


eligible_f444w = np.array(np.where(open_matched_sources['npix_f444w'] >= 32))


for i in range(len(eligible_f444w[0])): # run for the length of sources being looked at

    if open_matched_sources['npix_f444w'][i] >= 32:
        
        width_f115w = (open_matched_sources['xmax_f115w'][i] - open_matched_sources['xmin_f115w'][i]) * 1.5
        width_f150w = (open_matched_sources['xmax_f150w'][i] - open_matched_sources['xmin_f150w'][i]) * 1.5
        width_f277w = (open_matched_sources['xmax_f277w'][i] - open_matched_sources['xmin_f277w'][i]) * 1.5
        width_f444w = (open_matched_sources['xmax_f444w'][i] - open_matched_sources['xmin_f444w'][i]) * 1.5

        height_f115w = (open_matched_sources['ymax_f115w'][i] - open_matched_sources['ymin_f115w'][i]) * 1.5
        height_f150w = (open_matched_sources['ymax_f150w'][i] - open_matched_sources['ymin_f150w'][i]) * 1.5
        height_f277w = (open_matched_sources['ymax_f277w'][i] - open_matched_sources['ymin_f277w'][i]) * 1.5
        height_f444w = (open_matched_sources['ymax_f444w'][i] - open_matched_sources['ymin_f444w'][i]) * 1.5
        
        width = max(width_f115w, width_f150w, width_f277w, width_f444w)
        height = max(height_f115w, height_f150w, height_f277w, height_f444w)
        

        cutout_width = np.sqrt(width**2 + height**2)
            
                
        f444w_list.append("image " + str(i) + ". no. of pixels: " + str(open_matched_sources['npix_f444w'][i]))
        
        small_cutout_r = Cutout2D(r, (open_matched_sources['x_f444w'][i], open_matched_sources['y_f444w'][i]), size=(cutout_width,cutout_width))
        small_cutout_g = Cutout2D(g, (open_matched_sources['x_f444w'][i], open_matched_sources['y_f444w'][i]), size=(cutout_width,cutout_width))
        small_cutout_b = Cutout2D(blue, (open_matched_sources['x_f444w'][i], open_matched_sources['y_f444w'][i]), size=(cutout_width,cutout_width))
        
        small_r_name = small_cutout_r.data
        small_g_name = small_cutout_g.data
        small_i_name = small_cutout_b.data
        
        small_r_name = small_r_name.astype('float64')
        small_g_name = small_g_name.astype('float64')
        small_i_name = small_i_name.astype('float64')
        
        
        plt.figure()
        fig, ax = plt.subplots(figsize=(8, 8)) #make bigger
        small_rgb_default = make_lupton_rgb(small_r_name*1.1, small_g_name*0.75, small_i_name*0.4, minimum=-0.025, Q=10, stretch=0.5, filename=("CEERS1 RGB Cutout f444w" + str(i) + ".png"  ))
        im = ax.imshow(small_rgb_default,  interpolation = 'nearest', origin = 'lower')
 

#%%








   

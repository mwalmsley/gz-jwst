#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 14:37:00 2023

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
Reads in the fits files of different image bands from JWST NIRCam.

They correspond to the red, green, and blue channels for rgb image respectively. f115w and f150w are
combined to make blue. f277w is green, and f444w is red. Here, make_lupton_rgb is used to combine and 
make into colour images.

i_name is f150w
i_name_2 is f115w
g_name is f227w
r_name is f444w

'''

i_name = '/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/COSMOS-Web_NIRCam_DR0.2/mosaic_nircam_f150w_COSMOS-Web_60mas_v0_2_i2d.fits'
i_name_2 = '/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/COSMOS-Web_NIRCam_DR0.2/mosaic_nircam_f115w_COSMOS-Web_60mas_v0_2_i2d.fits'
g_name = '/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/COSMOS-Web_NIRCam_DR0.2/mosaic_nircam_f277w_COSMOS-Web_60mas_v0_2_i2d.fits'
r_name = '/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/COSMOS-Web_NIRCam_DR0.2/mosaic_nircam_f444w_COSMOS-Web_60mas_v0_2_i2d.fits'




'''
Data is being taken from each of the channels respectively.
*blue* is the f1154 and f150w bands combined.
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

*filename *saves the image in the directory of this code file.
'''
fig, ax = plt.subplots(figsize=(8, 8))
rgb = make_lupton_rgb(r*1.1, g*0.75, blue*0.4,  minimum=-0.025, Q=10, stretch=0.5, filename="COSMOS Colour.png") #add noise
im = ax.imshow(rgb, interpolation = 'nearest', origin = 'lower')


#%%



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
fits_filename = 'RGB_COSMOS_FITS_IMAGE.fits'

# Combine the desktop path and filename to get the full path
fits_filepath = os.path.join(desktop_path, fits_filename)

# Write the FITS file to the specified path
hdul.writeto(fits_filepath, overwrite=True)

print(f"FITS file saved to: {fits_filepath}")





'''
opens csv files of the sources that have been detected from the fulss size COSMOS image with a 5-sigma
detection threshold.
'''

open_f115w_sources = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/match1_only.csv") # 
open_f150w_sources = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/match2_only.csv") #
open_f277w_sources = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/match3_only.csv") # 
open_f444w_sources = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/match4_only.csv") # 
open_matched_sources = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/matched_all.csv") # 






'''
full sized colour image with the raised noise floor.
'''
plt.figure()
fig, ax = plt.subplots(figsize=(10, 10)) #make bigger
rgb_default = make_lupton_rgb(r*1.1, g*0.75, blue*0.4, minimum=-0.025, Q=10, stretch=0.5, filename="COSMOS colour (raised noise floor).png")
im = ax.imshow(rgb_default, interpolation = 'nearest', origin = 'lower')





#start running code in chunks





'''
plots ellipses onto full rgb colour image. Each band corresponds to a specific colour of ellipse.
So we know which source was detected in one band, or in all the bands for example.

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
    
    
for i in range(len(open_matched_sources)):
    e = Ellipse(xy=(open_matched_sources['x_all'][i], open_matched_sources['y_all'][i]),
                width = 6*open_matched_sources['a_all'][i],
                height = 6*open_matched_sources['b_all'][i],
                angle = open_matched_sources['theta_all'][i]* 180. / np.pi, lw=1)

    e.set_facecolor('none')
    e.set_edgecolor('white')
    ax.add_artist(e)
    

ax.set_axis_off()
fig.savefig("RGB Image with ellipses.png", bbox_inches='tight', pad_inches = 0, dpi=2384.2)






'''
This is plotiing the first 100 images from the table *open matched sources*, whihch are sources that appear 
in all four JWST bands. The images are rgb, and are made from each individual colour band for each image.
'''



plt.figure()
fig, ax = plt.subplots(figsize=(10, 10)) #make bigger


for i in range(2): # run for the length of sources being looked at
    small_cutout_r = Cutout2D(r, (open_matched_sources['x_all'][i], open_matched_sources['y_all'][i]), size=(200,200))
    small_cutout_g = Cutout2D(g, (open_matched_sources['x_all'][i], open_matched_sources['y_all'][i]), size=(200,200))
    small_cutout_b = Cutout2D(blue, (open_matched_sources['x_all'][i], open_matched_sources['y_all'][i]), size=(200,200))
    
    small_r_name = small_cutout_r.data
    small_g_name = small_cutout_g.data
    small_i_name = small_cutout_b.data
    
    small_r_name = small_r_name.astype('float64')
    small_g_name = small_g_name.astype('float64')
    small_i_name = small_i_name.astype('float64')
    
    
    plt.figure()
    fig, ax = plt.subplots(figsize=(8, 8)) #make bigger
    small_rgb_default = make_lupton_rgb(small_r_name*1.1, small_g_name*0.75, small_i_name*0.4, minimum=-0.025, Q=10, stretch=0.5, filename=("RGB Cutout " + str(i) + ".png"  ))
    im = ax.imshow(small_rgb_default,  interpolation = 'nearest', origin = 'lower')
    
    
    
#%%


'''
Open table that contains sources from objects that are exclusively in their respective bands.
'''

all_f115w_sources = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/F115w Objects.csv") # 
all_f150w_sources = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/F150w Objects.csv") #
all_f277w_sources = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/F277w Objects.csv") # 
all_f444w_sources = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/F444w Objects.csv") # 


all_sources_total = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/All_Detected_Sources.csv")
#%%

'''
gets a list of the ra and dec of the sources in each band separately.
the x and y coordinates are then converted to ra and dec and are added to arrays blue_band_x_coords, and 
blue_band_y_coords, (and their respective variables).
'''

#get x and y coords of of all sources in each band separately
#convert these x and y into ra and dec and add them to the table... 
#match all the tables based on the ra and dec of each said table
#then match that table with the ra and dec of the x-ray sources tables.

#append the x and y to array/list and add to table column

blue_band_x_coords = []
blue_band_y_coords = []

for i in range(len(all_f115w_sources)):
    blue_band_image = fits.open("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/COSMOS-Web_NIRCam_DR0.2/mosaic_nircam_f115w_COSMOS-Web_60mas_v0_2_i2d.fits")
    blue_band_image_coordinates = WCS(blue_band_image[1].header)
    blue_band_ra_and_dec = blue_band_image_coordinates.pixel_to_world(all_f115w_sources['x_f115w'][i], all_f115w_sources['y_f115w'][i])
    blue_band_x_coords.append(blue_band_ra_and_dec.ra.degree)
    blue_band_y_coords.append(blue_band_ra_and_dec.dec.degree)

    #print("PIXEL TO RA and DEC")
    #print(blue_band_ra_and_dec)

blue_band_x_coords = np.array(blue_band_x_coords)
blue_band_y_coords = np.array(blue_band_y_coords)




bluer_band_x_coords = []
bluer_band_y_coords = []

for i in range(len(all_f150w_sources)):
    bluer_band_image = fits.open("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/COSMOS-Web_NIRCam_DR0.2/mosaic_nircam_f150w_COSMOS-Web_60mas_v0_2_i2d.fits")
    bluer_band_image_coordinates = WCS(bluer_band_image[1].header)
    bluer_band_ra_and_dec = bluer_band_image_coordinates.pixel_to_world(all_f150w_sources['x_f150w'][i], all_f150w_sources['y_f150w'][i])
    bluer_band_x_coords.append(bluer_band_ra_and_dec.ra.degree)
    bluer_band_y_coords.append(bluer_band_ra_and_dec.dec.degree)

    #print("PIXEL TO RA and DEC")
    #print(str(bluer_band_ra_and_dec))  

bluer_band_x_coords = np.array(bluer_band_x_coords)
bluer_band_y_coords = np.array(bluer_band_y_coords)






green_band_x_coords = []
green_band_y_coords = []

for i in range(len(all_f277w_sources)):
    green_band_image = fits.open("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/COSMOS-Web_NIRCam_DR0.2/mosaic_nircam_f277w_COSMOS-Web_60mas_v0_2_i2d.fits")
    green_band_image_coordinates = WCS(green_band_image[1].header)
    green_band_ra_and_dec = green_band_image_coordinates.pixel_to_world(all_f277w_sources['x_f277w'][i], all_f277w_sources['y_f277w'][i])
    green_band_x_coords.append(green_band_ra_and_dec.ra.degree)
    green_band_y_coords.append(green_band_ra_and_dec.dec.degree)

    #print("PIXEL TO RA and DEC")
    #print(str(green_band_ra_and_dec))  

green_band_x_coords = np.array(green_band_x_coords)
green_band_y_coords = np.array(green_band_y_coords)






#%%
red_band_x_coords = []
red_band_y_coords = []


for i in range(len(all_f444w_sources)):
    red_band_image = fits.open("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/COSMOS-Web_NIRCam_DR0.2/mosaic_nircam_f444w_COSMOS-Web_60mas_v0_2_i2d.fits")
    red_band_image_coordinates = WCS(red_band_image[1].header)
    red_band_ra_and_dec = red_band_image_coordinates.pixel_to_world(all_f444w_sources['x_f444w'][i], all_f444w_sources['y_f444w'][i])    
    red_band_x_coords.append(red_band_ra_and_dec.ra.degree)
    red_band_y_coords.append(red_band_ra_and_dec.dec.degree)

    #print("PIXEL TO RA and DEC")
    #print(str(red_band_ra_and_dec))
    
red_band_x_coords = np.array(red_band_x_coords)
red_band_y_coords = np.array(red_band_y_coords)



#%%


all_band_x_coords = []
all_band_y_coords = []


for i in range(len(all_sources_total)):
    
    all_band_image = fits.open("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/COSMOS-Web_NIRCam_DR0.2/mosaic_nircam_f277w_COSMOS-Web_60mas_v0_2_i2d.fits")
    all_band_image_coordinates = WCS(all_band_image[1].header)
    all_band_ra_and_dec = all_band_image_coordinates.pixel_to_world(all_sources_total['x_all'][i], all_sources_total['y_all'][i])  
    all_band_x_coords.append(all_band_ra_and_dec.ra.degree)
    all_band_y_coords.append(all_band_ra_and_dec.dec.degree)

    #print("PIXEL TO RA and DEC")
    #print(str(red_band_ra_and_dec))
    
all_band_x_coords = np.array(all_band_x_coords)
all_band_y_coords = np.array(all_band_y_coords)



#%%

'''
now each table with sources in a give band have two new columns appended to said tables. the ra and dec
of each source.
'''

#add new columns to to fits tables.
#This is for blue-band. 

blue_sources = pd.read_csv('/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/F115w Objects.csv') # make false if no headers
blue_df = pd.DataFrame(blue_sources)

# Declare a list that is to be converted into a column
RA_blue = blue_band_x_coords
DEC_blue = blue_band_y_coords

# Using 'Address' as the column name
# and equating it to the list
blue_df['RA_f115w'] = RA_blue
blue_df['DEC_f115w'] = DEC_blue

print(blue_df)

#save csv
blue_df.to_csv('f115w_only_ra_and_dec.csv')






#%%



#add new columns to to fits tables.
#This is for bluer-band. 

bluer_sources = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/F150w Objects.csv") # make false if no headers
bluer_df = pd.DataFrame(bluer_sources)

# Declare a list that is to be converted into a column
RA_bluer = bluer_band_x_coords
DEC_bluer = bluer_band_y_coords

# Using 'Address' as the column name
# and equating it to the list
bluer_df['RA_f150w'] = RA_bluer
bluer_df['DEC_f150w'] = DEC_bluer

print(bluer_df)

#save csv
bluer_df.to_csv('f150w_only_ra_and_dec.csv')





#%%



#add new columns to to fits tables.
#This is for green-band. 

green_sources = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/F277w Objects.csv") # make false if no headers
green_df = pd.DataFrame(green_sources)

# Declare a list that is to be converted into a column
RA_green = green_band_x_coords
DEC_green = green_band_y_coords

# Using 'Address' as the column name
# and equating it to the list
green_df['RA_f277w'] = RA_green
green_df['DEC_f277w'] = DEC_green

print(green_df)

#save csv
green_df.to_csv('f277w_only_ra_and_dec.csv')





#%%



#add new columns to to fits tables.
#This is for red-band. 

red_sources = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/F444w Objects.csv") # make false if no headers
red_df = pd.DataFrame(red_sources)

# Declare a list that is to be converted into a column
RA_red = red_band_x_coords
DEC_red = red_band_y_coords

# Using 'Address' as the column name
# and equating it to the list
red_df['RA_f444w'] = RA_red
red_df['DEC_f444w'] = DEC_red

print(red_df)

#save csv
red_df.to_csv('f444w_only_ra_and_dec.csv')





#%%


#add new columns to to fits tables.
#This is for all-detections. 
all_sources = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/All_Detected_Sources.csv")

all_df = pd.DataFrame(all_sources)

# Declare a list that is to be converted into a column
RA_all = all_band_x_coords
DEC_all = all_band_y_coords

# Using 'Address' as the column name
# and equating it to the list
all_df['RA_all'] = RA_all
all_df['DEC_all'] = DEC_all

print(all_df)

#save csv
all_df.to_csv('all_bands_with_coords.csv')







#%%

'''
Tables are being opened which contain potenatial AGN sources, or more accurately, X-ray sources.

A colour cutout is made for each source. And All x-ray sources are plotted.
'''

pot_xray_sources_green = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Potential_AGN_green.csv")
pot_xray_sources_f150 = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Potential_AGN_f150.csv")
pot_xray_sources_f115 = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Potential_AGN_f115.csv")


all_pot_xray_sources = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/ALL_XRAY_SOURCES.csv")
all_sources = pd.read_csv("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/All_Detected_Sources.csv")


plt.figure()
fig, ax = plt.subplots(figsize=(10, 10)) #make bigger


for i in range(len(all_pot_xray_sources)): # run for the length of sources being looked at
    xray_small_cutout_r = Cutout2D(r, (all_pot_xray_sources['x_all'][i], all_pot_xray_sources['y_all'][i]), size=(200,200))
    xray_small_cutout_g = Cutout2D(g, (all_pot_xray_sources['x_all'][i], all_pot_xray_sources['y_all'][i]), size=(200,200))
    xray_small_cutout_b = Cutout2D(blue, (all_pot_xray_sources['x_all'][i], all_pot_xray_sources['y_all'][i]), size=(200,200))
    
    xray_small_r_name = xray_small_cutout_r.data
    xray_small_g_name = xray_small_cutout_g.data
    xray_small_i_name = xray_small_cutout_b.data
    
    xray_small_r_name = xray_small_r_name.astype('float64')
    xray_small_g_name = xray_small_g_name.astype('float64')
    xray_small_i_name = xray_small_i_name.astype('float64')
    
    
    plt.figure()
    fig, ax = plt.subplots(figsize=(8, 8)) #make bigger
    xray_small_rgb_default = make_lupton_rgb(xray_small_r_name*1.1, xray_small_g_name*0.75, xray_small_i_name*0.4, minimum=-0.025, Q=10, stretch=0.5, filename=("RGB XRAY Cutouts " + str(i) + ".png"  ))
    im = ax.imshow(xray_small_rgb_default,  interpolation = 'nearest', origin = 'lower')
    

#%%


'''
conda cheat sheet.
save parameters of large colour image I have.
For galaxy zoo, brighten images a little (perhaps maybe by 1.1 or 1.2).
Make colour image so that you can ever so slightly see the noise.
scale for whole image rather than cutout.
segmentation maps sep. desaturate sky pixels (which are noise).
Make all objects have pixel value of zero, so just get map of sky.
See source shapes in seg map

Outer-join, include things that are matched and thinngs that aren't matched.'

plot all matches
remove duplicates
plot what band each source shows up in.
'''




'''
whatever is in 1, 2, 3, or 4?
An individual colour for each right?
Make a new column in topcat for this?
'''








#star class source extractor
#segmentation map sep
#make multiple images of cutouts to check which settings are best

'''
select stars and galaxies properly with sep.
No noise selection.
no shredding galaxies/stars into small sources.
3/4 images in each band.
image with bright star, image with bright galaxy, edge image, random location.
subplots.

astropy lupton colour images, "preparing RGB images from CCD data"
make colour images
see detail in faint galaxies, don't saturate bright sources.

ceers data release images, which bands went into which images.
pick released ceers image, try and make my image look like that.

'''




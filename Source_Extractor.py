#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 11:02:08 2023

@author: husmak
"""

'''
sep - Source Extraction and Photometry.

Example of using SEP to detect objects in an image, and perfom
basic aperture photometry.
'''

import sep
import sys
import fitsio
import PIL.Image
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.image as mpimg
from matplotlib.patches import Ellipse
from astropy.nddata import Cutout2D



'''
Reads in the fits files of different image bands from JWST NIRCam.
Band names are prefixes of the following four variables.
'''
f115w_data = fitsio.read("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/COSMOS-Web_NIRCam_DR0.2/mosaic_nircam_f115w_COSMOS-Web_60mas_v0_2_i2d.fits")
f150w_data = fitsio.read("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/COSMOS-Web_NIRCam_DR0.2/mosaic_nircam_f150w_COSMOS-Web_60mas_v0_2_i2d.fits")
f277w_data = fitsio.read("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/COSMOS-Web_NIRCam_DR0.2/mosaic_nircam_f277w_COSMOS-Web_60mas_v0_2_i2d.fits")
f444w_data = fitsio.read("/Users/husmak/iCloud Drive (Archive)/Desktop/Source_Extractor/COSMOS-Web_NIRCam_DR0.2/mosaic_nircam_f444w_COSMOS-Web_60mas_v0_2_i2d.fits")





'''
Plots each of the images that I have read in. 
mean_ and std_ refer to the mean and std of the respective image. These are used for
the vmin and vmax of the image.
'''
mean_f115w, std_f115w = np.mean(f115w_data), np.std(f115w_data)
plt.figure()
plt.imshow(f115w_data, interpolation = 'nearest', cmap = 'gray', vmin = mean_f115w-std_f115w, vmax = mean_f115w+std_f115w, origin = 'lower')
plt.colorbar()

mean_f150w, std_f150w = np.mean(f150w_data), np.std(f150w_data)
plt.figure()
plt.imshow(f150w_data, interpolation = 'nearest', cmap = 'gray', vmin = mean_f150w-std_f150w, vmax = mean_f115w+std_f150w, origin = 'lower')
plt.colorbar()

mean_f277w, std_f277w = np.mean(f277w_data), np.std(f277w_data)
plt.figure()
plt.imshow(f277w_data, interpolation = 'nearest', cmap = 'gray', vmin = mean_f277w-std_f277w, vmax = mean_f277w+std_f277w, origin = 'lower')
plt.colorbar()

mean_f444w, std_f444w = np.mean(f444w_data), np.std(f444w_data)
plt.figure()
plt.imshow(f444w_data, interpolation = 'nearest', cmap = 'gray', vmin = mean_f444w-std_f444w, vmax = mean_f444w+std_f444w, origin = 'lower')
plt.colorbar()


'''
Here we generate the backgrounds of our images, which we will soon remove.
'''
f115w_bkg = sep.Background(f115w_data)
f150w_bkg = sep.Background(f150w_data)
f277w_bkg = sep.Background(f277w_data)
f444w_bkg = sep.Background(f444w_data)


#look at documentation if we want to mask the data and have more options for background subtraction


'''
We print the get global mean and noise of the image background... (units?).
Not sure what the noise is, or why the mean is zero.
'''
print(f115w_bkg.globalback)
print(f115w_bkg.globalrms)

print(f150w_bkg.globalback)
print(f150w_bkg.globalrms)

print(f277w_bkg.globalback)
print(f277w_bkg.globalrms)

print(f444w_bkg.globalback)
print(f444w_bkg.globalrms)




'''
Here we evaluate the backgrounds as 2D arrays, they are same size as original image.
'''
f115w_bkg_array = np.array(f115w_bkg)
f150w_bkg_array = np.array(f150w_bkg)
f277w_bkg_array = np.array(f277w_bkg)
f444w_bkg_array = np.array(f444w_bkg)


'''
The backgrounds are being plotted, just for illustration purposes.
'''

plt.figure()
plt.imshow(f115w_bkg_array, interpolation = 'nearest', cmap = 'gray', origin = 'lower')
plt.colorbar()

plt.figure()
plt.imshow(f150w_bkg_array, interpolation = 'nearest', cmap = 'gray', origin = 'lower')
plt.colorbar()

plt.figure()
plt.imshow(f277w_bkg_array, interpolation = 'nearest', cmap = 'gray', origin = 'lower')
plt.colorbar()

plt.figure()
plt.imshow(f444w_bkg_array, interpolation = 'nearest', cmap = 'gray', origin = 'lower')
plt.colorbar()



'''
Now the background *noises* are being evaluated as 2D arrays, also same size as original images.
'''
f115w_bkg_rms = f115w_bkg.rms()
f150w_bkg_rms = f150w_bkg.rms()
f277w_bkg_rms = f277w_bkg.rms()
f444w_bkg_rms = f444w_bkg.rms()


'''
Plotting background noises, again, just for illustration purposes.
'''
plt.figure()
plt.imshow(f115w_bkg_rms, interpolation = 'nearest', cmap='gray', origin='lower')
plt.colorbar()

plt.figure()
plt.imshow(f150w_bkg_rms, interpolation = 'nearest', cmap='gray', origin='lower')
plt.colorbar()

plt.figure()
plt.imshow(f277w_bkg_rms, interpolation = 'nearest', cmap='gray', origin='lower')
plt.colorbar()

plt.figure()
plt.imshow(f444w_bkg_rms, interpolation = 'nearest', cmap='gray', origin='lower')
plt.colorbar()




'''
background is being subtracted from the image.
'''
f115w_data_sub = f115w_data - f115w_bkg
f150w_data_sub = f150w_data - f150w_bkg
f277w_data_sub = f277w_data - f277w_bkg
f444w_data_sub = f444w_data - f444w_bkg





'''
Background has been subtracted, so now sources are being extracted. Here it is all sources that are 
5-sigma detections and above.
'''

f115w_objects = sep.extract(f115w_data_sub, 5, err = f115w_bkg.globalrms)
f150w_objects = sep.extract(f150w_data_sub, 5, err = f150w_bkg.globalrms)
f277w_objects = sep.extract(f277w_data_sub, 5, err = f277w_bkg.globalrms)
f444w_objects = sep.extract(f444w_data_sub, 5, err = f444w_bkg.globalrms)


'''
These are segmentation maps that show all the sources detected in their respective bands.
'''
f115w_segmap = sep.extract(f115w_data_sub, 5, err = f115w_bkg.globalrms, segmentation_map=True)
f150w_segmap = sep.extract(f150w_data_sub, 5, err = f150w_bkg.globalrms, segmentation_map=True)
f277w_segmap = sep.extract(f277w_data_sub, 5, err = f277w_bkg.globalrms, segmentation_map=True)
f444w_segmap = sep.extract(f444w_data_sub, 5, err = f444w_bkg.globalrms, segmentation_map=True)

'''
segmentation maps of each individual band are plotted.

Then a *total* segmentation map is plotted, a segmap with all the sources from all bands on it.
'''

fig, ax = plt.subplots(figsize=(10, 10))
plt.imshow(f115w_segmap[1], origin='lower')
fig.savefig("b_segmap.png", dpi=800)

fig, ax = plt.subplots(figsize=(8, 8))
plt.imshow(f150w_segmap[1], origin='lower')
fig.savefig("b2_segmap.png", dpi=800)

fig, ax = plt.subplots(figsize=(8, 8))
plt.imshow(f277w_segmap[1], origin='lower')
fig.savefig("g_segmap.png", dpi=800)

fig, ax = plt.subplots(figsize=(8, 8))
plt.imshow(f444w_segmap[1], origin='lower')
fig.savefig("r_segmap.png", dpi=800)


fig, ax = plt.subplots(figsize=(8, 8))
plt.imshow(f115w_segmap[1], origin='lower')
plt.imshow(f150w_segmap[1], origin='lower')
plt.imshow(f277w_segmap[1], origin='lower')
plt.imshow(f444w_segmap[1], origin='lower')
fig.savefig("total_segmap.png", dpi=800)


'''
prints how many objects were detected in each band.
'''
print(len(f115w_objects[0]))
print(len(f150w_objects[0]))
print(len(f277w_objects[0]))
print(len(f444w_objects[0]))






'''
background subracted images in each band are being plotted with ellipses over each source that has
been detected in the respective band.

e.g. the f115w band image is plotted, with ellipses over each source that was detected in the f115w band.

This applies to each of the 4 bands.
'''


#plot background subtracted images
fig, ax = plt.subplots()
mean_f115w, std_f115w = np.mean(f115w_data_sub), np.std(f115w_data_sub)
im = ax.imshow(f115w_data_sub, interpolation = 'nearest', cmap = 'gray',
               vmin = mean_f115w-std_f115w, vmax = mean_f115w+std_f115w, origin = 'lower')


#plot an ellipse for each object (what do the letters mean?)
for i in range(len(f115w_objects)):
    e = Ellipse(xy=(f115w_objects['x'][i], f115w_objects['y'][i]),
                width = 6*f115w_objects['a'][i],
                height = 6*f115w_objects['b'][i],
                angle = f115w_objects['theta'][i]* 180. / np.pi)
    
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)


fig.savefig("f115w background subtraction.png")





fig, ax = plt.subplots()
mean_f150w, std_f150w = np.mean(f150w_data_sub), np.std(f150w_data_sub)
im = ax.imshow(f150w_data_sub, interpolation = 'nearest', cmap = 'gray',
               vmin = mean_f150w-std_f150w, vmax = mean_f150w+std_f150w, origin = 'lower')


for i in range(len(f150w_objects)):
    e = Ellipse(xy=(f150w_objects['x'][i], f150w_objects['y'][i]),
                width = 6*f150w_objects['a'][i],
                height = 6*f150w_objects['b'][i],
                angle = f150w_objects['theta'][i]* 180. / np.pi)
    
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)

fig.savefig("f150w background subtraction.png")





fig, ax = plt.subplots()
mean_277w, std_277w = np.mean(f277w_data_sub), np.std(f277w_data_sub)
im = ax.imshow(f277w_data_sub, interpolation = 'nearest', cmap = 'gray',
               vmin = mean_277w-std_f277w, vmax = mean_f277w+std_f277w, origin = 'lower')

for i in range(len(f277w_objects)):
    e = Ellipse(xy=(f277w_objects['x'][i], f277w_objects['y'][i]),
                width = 6*f277w_objects['a'][i],
                height = 6*f277w_objects['b'][i],
                angle = f277w_objects['theta'][i]* 180. / np.pi)
    
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)


fig.savefig("f277w background subtraction.png", dpi=800)





fig, ax = plt.subplots()
mean_f444w, std_f444w = np.mean(f444w_data_sub), np.std(f444w_data_sub)
im = ax.imshow(f444w_data_sub, interpolation = 'nearest', cmap = 'gray',
               vmin = mean_f444w-std_f444w, vmax = mean_f444w+std_f444w, origin = 'lower')


for i in range(len(f444w_objects)):
    e = Ellipse(xy=(f444w_objects['x'][i], f444w_objects['y'][i]),
                width = 6*f444w_objects['a'][i],
                height = 6*f444w_objects['b'][i],
                angle = f444w_objects['theta'][i]* 180. / np.pi)
    
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)


fig.savefig("f444w background subtraction.png")




'''
Dataframes of all the objects detected in each band are created. These are then saved
as csv files for use in topcat, or any other table work I need to be done.
'''

f115w_df = pd.DataFrame(f115w_objects)

f115w_colnames = f115w_df.columns
f115w_newcolnames = ["%s_f115w" % q for q in f115w_colnames]
f115w_df.columns = f115w_newcolnames

f115w_df.to_csv("F115w Objects.csv")

read_f115w_df = pd.read_csv("F115w Objects.csv") # make false if no headers
f115w_row_number = 500
f115w_df_row = read_f115w_df.iloc[f115w_row_number,:]
print(f115w_df_row['x_f115w'], f115w_df_row['y_f115w'])




f150w_df = pd.DataFrame(f150w_objects)

f150w_colnames = f150w_df.columns
f150w_newcolnames = ["%s_f150w" % q for q in f150w_colnames]
f150w_df.columns = f150w_newcolnames

f150w_df.to_csv("F150w Objects.csv")

read_f150w_df = pd.read_csv("F150w Objects.csv") # make false if no headers
f150w_row_number = 500
f150w_df_row = read_f150w_df.iloc[f150w_row_number,:]
print(f150w_df_row['x_f150w'], f150w_df_row['y_f150w'])




f277w_df = pd.DataFrame(f277w_objects)

f277w_colnames = f277w_df.columns
f277w_newcolnames = ["%s_f277w" % q for q in f150w_colnames]
f277w_df.columns = f277w_newcolnames

f277w_df.to_csv("F277w Objects.csv")

read_f277w_df = pd.read_csv("F277w Objects.csv") # make false if no headers
f277w_row_number = 500
f277w_df_row = read_f277w_df.iloc[f277w_row_number,:]
print(f277w_df_row['x_f277w'], f277w_df_row['y_f277w'])




f444w_df = pd.DataFrame(f444w_objects)

f444w_colnames = f444w_df.columns
f444w_newcolnames = ["%s_f444w" % q for q in f150w_colnames]
f444w_df.columns = f444w_newcolnames

f444w_df.to_csv("F444w Objects.csv")

read_f444w_df = pd.read_csv("F444w Objects.csv") # make false if no headers
f444w_row_number = 500
f444w_df_row = read_f444w_df.iloc[f444w_row_number,:]
print(f444w_df_row['x_f444w'], f444w_df_row['y_f444w'])



#%%

'''
source extractor on all 4 bands individually.
want a table that is the union of all 4 bands
then tables that are in each band individually, and in each band but not others, 
plot ellipses that show in each band
'''


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







maybe open all files at once? DONE
export to csv with column headers (probably need pandas) DONE
focus on NirCam data rather than Miri. Work on one band at a time first, then combine later. DONE
select source across each band. DONE
zoom into image section, one section on edge, one random, one with big source. TO Do WITH TOPCAT coords. DONE
change thresholds. DONE
sources identifiable by ra and dec in source extractor. SEMI-DONE (only have pixels, then converted to ra and dec)

data reduction...

change size of ellipses.
check email from brooke about this.
galaxy zoo talk. JWST CEERS metadat information
see what to reproduce from brooke's post.
find sources that are detectable in each band, if detected in each band, add fluxes. sep DOCUMENTATION. 


topcat


what we said we're going to do
what I've done
what we're going to do next
'''














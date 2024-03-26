#! /usr/bin/python3

## References:
# https://prancer.physics.louisville.edu/astrowiki/index.php/Image_processing_with_Python_and_SciPy
# https://learn.astropy.org/tutorials/FITS-images.html

## Import modules
import numpy as np
import sys
import os
from astropy.io import fits
from matplotlib import pyplot as plt

## Import files
image_file = sys.argv[1]
dark_lib = '/home/gergoe/Dokumentumok/supplementary_images/dark'
bias_lib = '/home/gergoe/Dokumentumok/supplementary_images/bias'
flat_lib = '/home/gergoe/Dokumentumok/supplementary_images/flat'
image_limit = 20000

## Read image
raw_image = fits.open(image_file)[0].data
image_data = np.array(raw_image, dtype=np.float64)
print(np.info(image_data))

## Build marster dark
dark_dir = os.listdir(dark_lib)
dark_files = [f for f in dark_dir if os.path.isfile(dark_lib+'/'+f) and (f.endswith('.fit') or f.endswith('.fits'))]
dark_list = []
for file in dark_files:
    image = fits.open(dark_lib+'/'+file)
    #print(np.info(image[0].data))
    dark_list.append(image[0].data)
dark_stack = np.array(dark_list)
master_dark = np.median(dark_stack, axis=0)

## Build marster bias
bias_dir = os.listdir(bias_lib)
bias_files = [f for f in bias_dir if os.path.isfile(bias_lib+'/'+f) and (f.endswith('.fit') or f.endswith('.fits'))]
bias_list = []
for file in bias_files:
    image = fits.open(bias_lib+'/'+file)
    #print(np.info(image[0].data))
    bias_list.append(image[0].data)
bias_stack = np.array(bias_list)
master_bias = np.median(bias_stack, axis=0)

## Build marster flat
flat_dir = os.listdir(flat_lib)
flat_files = [f for f in flat_dir if os.path.isfile(flat_lib+'/'+f) and (f.endswith('.fit') or f.endswith('.fits'))]
flat_list = []
for file in flat_files:
    image = fits.open(flat_lib+'/'+file)
    #print(np.info(image[0].data))
    flat_list.append(image[0].data)
flat_stack = np.array(flat_list)
master_flat = np.median(flat_stack, axis=0)

## Remove dark
# dark_corrected_image = raw_image - master_dark
dark_corrected_image = image_data - master_dark

## Remove bias
# bias_corrected_image = raw_image - master_bias
bias_corrected_image = dark_corrected_image - master_bias

## Correct flat
# final_image = bias_corrected_image / master_flat
final_image = bias_corrected_image / master_flat

## Troubleshooting
plt.imshow(final_image, aspect='equal') # , vmax=image_limit, vmin=0, origin='lower',cmap='Greys', 
plt.savefig('final_image.jpg',dpi=150.0, bbox_inches='tight', pad_inches=0.2)
plt.show()
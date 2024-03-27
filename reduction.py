#! /usr/bin/python3

## References:
# https://prancer.physics.louisville.edu/astrowiki/index.php/Image_processing_with_Python_and_SciPy
# https://learn.astropy.org/tutorials/FITS-images.html

## Import modules
import numpy as np
import sys
import os
import ccdproc
from astropy.io import fits
from astropy import units as u
from astropy.nddata import CCDData
from matplotlib import pyplot as plt

## Import files
image_file = sys.argv[1]
dark_lib = '/home/gergoe/Dokumentumok/supplementary_images/dark'
bias_lib = '/home/gergoe/Dokumentumok/supplementary_images/bias'
flat_lib = '/home/gergoe/Dokumentumok/supplementary_images/flat'

## Read image
raw_image = fits.open(image_file)[0].data
image_data = np.array(raw_image, dtype=np.float64)
print(np.info(image_data))

## Define a function for making a linear gray scale
def lingray(image_data, image_min=None, image_max=None):
    """
    Auxiliary function that specifies the linear gray scale.
    a and b are the cutoffs : if not specified, min and max are used
    """
    if image_min == None:
        image_min = np.min(image_data)
    if image_max == None:
        image_max = np.max(image_data)
    return 255.0 * (image_data - float(image_min)) / (image_max - image_min)

## Define a function for making a logarithmic gray scale
def loggray(x, a=None, b=None):
    """
    Auxiliary function that specifies the logarithmic gray scale.
    a and b are the cutoffs : if not specified, min and max are used
    """
    if a == None:
        a = np.min(x)
    if b == None:
        b = np.max(x)          
    linval = 10.0 + 990.0 * (x-float(a))/(b-a)
    return (np.log10(linval)-1.0) * 0.5 * 255.0

#grayscaled_image_data = lingray(image_data)

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
#master_flat = CCDData(master_flat, unit=u.electron)

## Remove dark
dark_corrected_image = image_data - master_dark

## Remove bias
bias_corrected_image = dark_corrected_image - master_bias

## Correct flat
flat_corrected_image = bias_corrected_image / master_flat

## Create final image by correctiong the grayscale
final_image = flat_corrected_image
#final_image = lingray(flat_corrected_image)
testitng_image = final_image
testing_name = 'final_image'
new_image_min = 0.
new_image_max = np.max(testitng_image)

## Troubleshooting
plt.imshow(testitng_image, aspect='equal', vmax=new_image_max, vmin=new_image_min) # ,  , , , , , origin='lower',  cmap='Greys'
plt.savefig(testing_name + '.jpg',dpi=150.0, bbox_inches='tight', pad_inches=0.2)
plt.show()
print(image_data)
#print(grayscaled_image_data)
print(master_dark)
print(master_bias)
print(master_flat)
print(final_image)

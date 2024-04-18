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
import math

## Import files
image_file = sys.argv[1]
dark_lib = '/home/gergoe/Képek/Atik/DARK'
bias_lib = '/home/gergoe/Képek/Atik/BIAS'
flat_lib = '/home/gergoe/Képek/Atik/FLAT'

## Read image
raw_image = fits.open(image_file)[0].data
raw_image_header = fits.open(image_file)[0].header
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
#master_flat = CCDData(master_flat, unit=u.electron)

## Remove dark
dark_corrected_image = image_data - master_dark

## Remove bias
bias_corrected_image = dark_corrected_image - master_bias

## Remove negative values from the image
#bias_corrected_image[bias_corrected_image<0] = 0

## Correct flat
flat_average = np.mean(master_flat)
flat_corrected_image = dark_corrected_image / master_flat

## Create final image by correctiong the grayscale
#final_image = flat_corrected_image * flat_average
#final_image = flat_corrected_image + math.fabs(np.min(flat_corrected_image))
#new_image_min = 0.
#new_image_max = np.max(final_image)
#gray_scaling = 65535 / np.max(final_image)
#corrected_image = gray_scaling * final_image

## Troubleshooting
plt.imshow(flat_corrected_image, aspect='equal', vmax=np.max(flat_corrected_image), vmin=np.min(flat_corrected_image), cmap='Greys', origin='lower') # 
plt.savefig('final_image.jpg',dpi=150.0, bbox_inches='tight', pad_inches=0.2)
plt.show()
print(flat_corrected_image)
print('final image min: ', np.min(flat_corrected_image))
print('final image max: ', np.max(flat_corrected_image))

outhdu = fits.PrimaryHDU(flat_corrected_image)
outhdr = outhdu.header
outhdr.extend(raw_image_header)
outhdu.writeto('new_file.fits', overwrite=True)

#fits.writeto(image_file, outhdulist, outhdr, overwrite=True)

'''# Chat GPT suggestion

from astropy.io import fits
import numpy as np
import os
import sys
from matplotlib import pyplot as plt

def process_astronomical_images(light_folder, dark_folder, bias_folder, flat_folder):
    """
    Process astronomical images including light, dark, bias, and flat frames.

    Args:
        light_folder (str): Path to the folder containing light images.
        dark_folder (str): Path to the folder containing dark images.
        bias_folder (str): Path to the folder containing bias images.
        flat_folder (str): Path to the folder containing flat images.

    Returns:
        numpy.ndarray: Processed light image.
    """
    # Function to combine images
    def combine_images(folder):
        images = []
        for filename in os.listdir(folder):
            if filename.endswith('.fits' or '.fit' or 'new'):
                with fits.open(os.path.join(folder, filename)) as hdul:
                    image_data = hdul[0].data
                    images.append(image_data)
        combined_image = np.median(images, axis=0)
        return combined_image

    # Combine dark, bias, and flat images
    dark_image = combine_images(dark_folder)
    bias_image = combine_images(bias_folder)
    flat_image = combine_images(flat_folder)

    # Process light images
    for filename in os.listdir(light_folder):
        if filename.endswith('.fits'):
            with fits.open(os.path.join(light_folder, filename)) as hdul:
                light_image = hdul[0].data
                # Subtract bias
                light_image -= bias_image
                # Subtract dark
                light_image -= dark_image
                # Divide by flat
                light_image /= flat_image
                # Additional processing steps can be added here
                # Save or return processed light image
                

                return light_image  # or whatever you want to return

# Example usage:
light_folder = sys.argv[1]
dark_folder = '/home/gergoe/Képek/Atik/DARK'
bias_folder = '/home/gergoe/Képek/Atik/BIAS'
flat_folder = '/home/gergoe/Képek/Atik/FLAT'

processed_light_image = process_astronomical_images(light_folder, dark_folder, bias_folder, flat_folder)

print(process_astronomical_images)'''
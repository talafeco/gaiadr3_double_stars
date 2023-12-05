#!/usr/local/bin/python3

"""

  Find stars and do aperture photometry in a FITS image

    Input a fits image and a base for output file names
    Find stars in the image
    Perform relative photometry on the stars
    Optionally output a CSV file with the photometry
    Optionally display the image with the stars marked    
    Ouput an AstroImageJ aperture file with FITS image coordinates
    Options are set by flags in the following code
    
  Defaults
  
    Will create a csv file of photometry 
    Will show the image with stars marked
    Seeks stars with FWHM ~ 8 pixels
    Applies a threshold of 8 sigma of the background to identify stars

  Dependencies
  
    Requires astropy io.fits, io.asci, stats
    Uses photutils IRAFStarFinder to find stars and perform photometry    

"""

import os
import sys
import numpy as np
import astropy.io.fits as pyfits
import astropy.io.ascii
from astropy.stats import sigma_clipped_stats
from photutils import IRAFStarFinder

# For visualization

import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils import CircularAperture


if len(sys.argv) == 1:
  print("Find stars in a fits image and create an AIJ apertures file.")
  sys.exit("Usage: fits_find_stars.py infile.fits output_basename")
  exit()
elif len(sys.argv) == 3:
  infits = sys.argv[1]
  outbase = sys.argv[2]
  ap_file = outbase+".apertures"
  csv_file = outbase+".csv"
else:
  sys.exit("Usage: fits_find_stars.py infile.fits output_basename")
  exit() 
  
# Set true for feedback
verbose_flag = True

# Set true for diagnostics
diagnostics_flag = False

# Set true for visualization
visuals_flag = True  

# Set true for csv speadsheet output
csv_flag = True

# Open the fits file readonly by default and create an input hdulist
inlist = pyfits.open(infits) 

# Assign the input header 
inhdr = inlist[0].header

# Assign image data to numpy array 
indata =  inlist[0].data.astype('float32')

# Get the dimensions for possible use later
xsize, ysize = indata.shape

# Image statistics
mean, median, sigma = sigma_clipped_stats(indata, sigma=3.0)

if verbose_flag:
  print("\n")
  print("Starting aperture photometry --\n")
  print("Image "+infits+" : ", xsize, " x ", ysize)
  print("Mean: ", mean)
  print("Median: ", median)
  print("Sigma: ", sigma)
  print("\n")
  
# Find stars based on photutils IRAFstarfinder  

# Search image for local density maxima 
# Peak amplitude greater than threshold above the local background
# Seek PSF full-width at half-maximum similar to the input fwhm 
# Calculate centroid, ellipticity, and sharpness from image moments
# Estimates background using starfind method

starfinder = IRAFStarFinder(fwhm=8.0, threshold=8.*sigma)  

# Identify sources meeting the criteria
# Subtract the median before applying the sigma threshold

sources = starfinder(indata - median) 
n_sources = len(sources['id'])
 
# Optionally print the table

if diagnostics_flag:
  print("\n",sources.colnames,"\n")
  for col in sources.colnames:  
    sources[col].info.format = '%.8g' 
  print(sources) 
  print(sources[1])
  print("\n")  
  print("Sources found: ",n_sources,"\n")

# Output the stars in a spreadsheet
# These data are in numpy array format
# Directions x and y are transposed from FITS coordinates
# Photutils coordinates may have offsets from origin

if csv_flag:
  print("\n")
  print("Writing photometry for ",n_sources," stars to ", csv_file,"\n")
  astropy.io.ascii.write(sources, csv_file, format='csv', overwrite=True)


# Extract pixel coordinates for sources
# Sources have ['mag'] relative to background that can discriminate 
# For example use "if (sources[i]['mag'] < -1.) :" 
# Transpose x and y for numpy to fits mapping

pixelxy = np.transpose((sources['xcentroid'], sources['ycentroid']))

# Optionally show image and apertures

if visuals_flag:
  pixelrad = 7.
  apertures = CircularAperture(pixelxy, r=pixelrad)
  norm = ImageNormalize(stretch=SqrtStretch())
  plt.imshow(indata, cmap='Greys', origin='lower', norm=norm)
  apertures.plot(color='blue', lw=1.5, alpha=0.5)
  plt.show()

# Close image file
 
inlist.close()


# Format the aij aperture file

# #AstroImageJ Saved Apertures
# #Sat Jun 22 22:19:34 EDT 2013
# .multiaperture.naperturesmax=500
# .multiaperture.isrefstar=false,false,true,true,true,true
# .multiaperture.xapertures=1473.4099,1580.1968,1539.5604,1593.1986,1511.0461,1547.8082
# .multiaperture.yapertures=1231.65,1272.3362,1156.1813,1132.1183,1280.981,1323.3054

# Open an file for AstroImageJ apertures
ap_fp = open(ap_file, 'w')
  
# Write some useful global parameters that might differ from the defaults
ap_fp.write(".multiaperture.naperturesmax=1000\n")

# Write the target or calibration flags for each star in the list
# Here we make them all targets

# Include offset between image FITS coordinates and internal photutils coordinates
x_offset = 1.
y_offset = 1.

ap_fp.write(".multiaperture.isrefstar=") 
for i in range(n_sources - 1):
  ap_fp.write("false,")
ap_fp.write("false\n")

# Write the x apertures for FITS coordinates in AstroImageJ
ap_fp.write(".multiaperture.xapertures=FITS")
for i in range(n_sources - 1):   
  x = sources[i]['xcentroid'] + x_offset
  x_line = "%7.2f," % (x,)
  ap_fp.write(x_line)
x = sources[n_sources - 1]['xcentroid']
x_line = "%7.2f\n" % (x,)
ap_fp.write(x_line)

# Write the y apertures for FITS coordinates in AstroImageJ
ap_fp.write(".multiaperture.yapertures=FITS")
for i in range(n_sources - 1):   
  y = sources[i]['ycentroid'] + y_offset
  y_line = "%7.2f," % (y,)
  ap_fp.write(y_line)
y = sources[n_sources - 1]['ycentroid'] 
y_line = "%7.2f\n" % (y,)
ap_fp.write(y_line)

# Close the apertures  file

ap_fp.close()
exit()

'''
This script opens the original FITS file, retrieves its header, creates a new FITS file with a primary HDU, updates the header of the new FITS file with the header from the original file, saves the new FITS file, and finally closes both files.

You can customize this script according to your needs, such as specifying different HDUs or modifying header keywords before saving the new FITS file. Additionally, you may want to handle errors and exceptions appropriately depending on your specific requirements and input files.
'''

from astropy.io import fits

# Open the FITS file from which you want to copy the header
original_hdulist = fits.open('original.fits')
original_header = original_hdulist[0].header

# Create a new FITS file with a primary HDU
new_hdulist = fits.PrimaryHDU()
new_header = new_hdulist.header

# Update the header of the new FITS file with the header from the original file
new_header.extend(original_header.cards)

# Save the new FITS file
new_hdulist.writeto('new.fits', overwrite=True)

# Close the FITS files
original_hdulist.close()
import rawpy
import numpy as np
import cv2
from astropy.io import fits
from astropy.wcs import WCS
from exif import Image as ExifImage
import os
import subprocess

def process_raw_image(file_path):
    # Step 1: Read Canon RAW file
    with rawpy.imread(file_path) as raw:
        rgb_image = raw.postprocess(output_bps=16)
    
    # Step 2: Convert to Grayscale (16-bit)
    grayscale_image = cv2.cvtColor(rgb_image, cv2.COLOR_RGB2GRAY)
    
    # Step 3: Star Detection
    # Use OpenCV SimpleBlobDetector or other algorithms for star detection
    params = cv2.SimpleBlobDetector_Params()
    params.filterByArea = True
    params.minArea = 5
    detector = cv2.SimpleBlobDetector_create(params)
    keypoints = detector.detect(grayscale_image)
    stars = np.array([[kp.pt[0], kp.pt[1]] for kp in keypoints])

    # Step 4: Call astrometry.net for plate-solving
    star_coords_path = 'temp_stars.txt'
    np.savetxt(star_coords_path, stars)
    output_wcs_file = 'temp.wcs'
    subprocess.run([
        'solve-field', file_path, '--no-plots', '--overwrite', 
        '--scale-units', 'degwidth', '--scale-low', '0.1', '--scale-high', '180',
        '--ra', '0', '--dec', '0', '--radius', '180', '--wcs', output_wcs_file,
        '--match', star_coords_path
    ])

    # Step 5: Read EXIF Data
    with open(file_path, 'rb') as img_file:
        exif_image = ExifImage(img_file)
        exif_data = {tag: exif_image.get(tag) for tag in exif_image.list_all()}

    # Step 6: Create a new FITS file
    fits_filename = os.path.splitext(os.path.basename(file_path))[0] + '.fit'
    hdu = fits.PrimaryHDU(grayscale_image)
    hdul = fits.HDUList([hdu])

    # Step 7: Write WCS and EXIF data to FITS
    try:
        wcs = WCS(output_wcs_file)
        hdu.header.update(wcs.to_header())
    except Exception as e:
        print(f"Error applying WCS: {e}")
    
    # Add EXIF data to FITS header
    for tag, value in exif_data.items():
        if isinstance(value, str) and len(value) <= 68:  # FITS header value length limit
            hdu.header[tag[:8].upper()] = value
    
    hdul.writeto(fits_filename, overwrite=True)
    print(f"Saved processed FITS file as {fits_filename}")

# Define the directory containing Canon RAW files
raw_image_directory = '/path/to/raw_images'

# Process each RAW image in the specified directory
for filename in os.listdir(raw_image_directory):
    if filename.lower().endswith('.cr2'):
        file_path = os.path.join(raw_image_directory, filename)
        process_raw_image(file_path)

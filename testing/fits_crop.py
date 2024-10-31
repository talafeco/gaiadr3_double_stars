from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
from PIL import Image, ImageDraw
import os

def crop_double_star_to_jpg_with_markers(fits_path, output_directory, star1_ra, star1_dec, star2_ra, star2_dec):
    # Load the FITS file to extract WCS and image data
    with fits.open(fits_path) as hdul:
        wcs = WCS(hdul[0].header)
        image_data = hdul[0].data
    
    # Define the coordinates for each star
    star1 = SkyCoord(ra=star1_ra, dec=star1_dec, unit="deg")
    star2 = SkyCoord(ra=star2_ra, dec=star2_dec, unit="deg")
    
    # Calculate the midpoint and separation between the stars
    midpoint = star1.midpoint(star2)
    separation = star1.separation(star2).arcsec
    
    # Define the frame size as three times the star separation (in arcseconds)
    frame_size = 3 * separation
    
    # Convert the midpoint to pixel coordinates
    midpoint_pix = wcs.world_to_pixel(midpoint)
    
    # Convert the frame size from arcseconds to pixels
    pixel_scale = np.mean(np.abs(wcs.pixel_scale_matrix)) * 3600  # arcsec/pixel
    frame_size_pix = int(frame_size / pixel_scale)
    
    # Define the cropping box in pixel coordinates
    half_size = frame_size_pix // 2
    left = int(midpoint_pix[0] - half_size)
    right = int(midpoint_pix[0] + half_size)
    bottom = int(midpoint_pix[1] - half_size)
    top = int(midpoint_pix[1] + half_size)
    
    # Crop the image data
    cropped_data = image_data[bottom:top, left:right]
    
    # Normalize the cropped data for better visualization in JPG
    cropped_data = cropped_data - np.min(cropped_data)
    cropped_data = (cropped_data / np.max(cropped_data) * 255).astype(np.uint8)
    
    # Convert the cropped image data to a PIL image
    cropped_image = Image.fromarray(cropped_data)
    
    # Rotate the image to make north-up and east-left if needed
    cropped_image = cropped_image.rotate(90, expand=True)  # Adjust rotation as needed
    
    # Recalculate WCS for the cropped image
    cropped_wcs = wcs.deepcopy()
    cropped_wcs.wcs.crpix[0] -= left
    cropped_wcs.wcs.crpix[1] -= bottom

    # Convert star positions to pixel coordinates in the cropped frame
    star1_pix = cropped_wcs.world_to_pixel(star1)
    star2_pix = cropped_wcs.world_to_pixel(star2)

    # Draw cross markers at each star position
    draw = ImageDraw.Draw(cropped_image)
    marker_size = 10  # Length of the cross lines in pixels
    
    # Draw markers for Star 1
    x1, y1 = star1_pix
    draw.line((x1 - marker_size, y1, x1 + marker_size, y1), fill="red", width=1)
    draw.line((x1, y1 - marker_size, x1, y1 + marker_size), fill="red", width=1)
    
    # Draw markers for Star 2
    x2, y2 = star2_pix
    draw.line((x2 - marker_size, y2, x2 + marker_size, y2), fill="blue", width=1)
    draw.line((x2, y2 - marker_size, x2, y2 + marker_size), fill="blue", width=1)
    
    # Save the cropped image in JPG format to the specified directory
    output_path = os.path.join(output_directory, "double_star_cropped_with_markers.jpg")
    cropped_image.save(output_path, format="JPEG")
    
    return output_path

# Example usage:
# output_path = crop_double_star_to_jpg_with_markers(
#     fits_path="path/to/your_image.fits",
#     output_directory="path/to/output_directory",
#     star1_ra=10.684, star1_dec=41.269,
#     star2_ra=10.686, star2_dec=41.270
# )
# print(f"Cropped image with markers saved to: {output_path}")

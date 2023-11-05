from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from matplotlib import pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import utils
from astropy.wcs import WCS

image_file = 'C:\Astro\\2023-09-05\Gaia_DR3_2869377048124227200\90s\Gaia_DR3_2869377048124227200_2023-09-06_03-05-47_90.00s_0005.new'
image_data = fits.open(image_file)

image_center_ra = 352.8750000
image_center_dec = 28.9247222

wcs_helix = WCS(image_data[0].header)

image = image_data[0].data

star_a = SkyCoord(352.861358469928 * u.deg, 28.975806180801026 * u.deg, frame='icrs')
star_b = SkyCoord(352.8493865739493 * u.deg, 28.981565231115656 * u.deg, frame='icrs')

star_a_pix = utils.skycoord_to_pixel(star_a, wcs_helix)
star_b_pix = utils.skycoord_to_pixel(star_b, wcs_helix)

fig = plt.figure(figsize=(10, 10), frameon=False) # 
ax = plt.subplot(projection=wcs_helix)
ax.arrow(image_center_ra, image_center_dec, 0, 0.016, 
         head_width=0, head_length=0, 
         fc='white', ec='white', width=0.0003, 
         transform=ax.get_transform('icrs'))
plt.text(image_center_ra, image_center_dec, '1 arcmin', 
         color='white',
         transform=ax.get_transform('icrs'))
plt.scatter(star_a_pix[0] + 30, star_a_pix[1], marker="_", s=50, color="grey", label='Main star')
plt.scatter(star_a_pix[0], star_a_pix[1] + 30, marker="|", s=50, color="grey", label='Main star')
plt.scatter(star_b_pix[0] + 30, star_b_pix[1], marker="_", s=50, color="grey", label='Companion star')
plt.scatter(star_b_pix[0], star_b_pix[1] + 30, marker="|", s=50, color="grey", label='Companion star')
plt.xlabel(image_center_ra)
plt.ylabel(image_center_dec)

overlay = ax.get_coords_overlay('icrs')
overlay.grid(color='grey', ls='dotted')
plt.imshow(image, origin='lower',cmap='grey', aspect='equal', vmax=2500, vmin=0) # , cmap='cividis'

plt.show()

plt.savefig('image.png', bbox_inches='tight')
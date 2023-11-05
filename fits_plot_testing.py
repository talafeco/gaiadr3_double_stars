import matplotlib.pyplot as plt

from astropy.visualization import astropy_mpl_style

plt.style.use(astropy_mpl_style)

from astropy.io import fits


image_data = fits.getdata('C:\Astro\\2023-09-05\Gaia_DR3_2869377048124227200\90s\Gaia_DR3_2869377048124227200_2023-09-06_03-01-11_90.00s_0002.new', ext=0)

print(image_data.shape)

plt.figure()


plt.imshow(image_data, cmap='gray', vmin=0, vmax=255)

plt.colorbar()

plt.show()
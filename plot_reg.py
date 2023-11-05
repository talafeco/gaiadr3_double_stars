# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from matplotlib import pyplot as plt

from regions import Regions


image_file = 'C:\Astro\\2023-09-05\Gaia_DR3_2869377048124227200\90s\Gaia_DR3_2869377048124227200_2023-09-06_03-01-11_90.00s_0002.new'
image_data = fits.getdata(image_file, ext=0, memmap=False)

fig, ax = plt.subplots()
ax.imshow(image_data, cmap='gray', vmin=0, vmax=2500)
ax.set_ylim([-0.5, 892.5])

region_file = get_pkg_data_filename('data/plot_image.reg',
                                    package='regions.io.ds9.tests')
regions = Regions.read(region_file, format='ds9')
for i, region in enumerate(regions):
    region.plot(ax=ax)

plt.show()
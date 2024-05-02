#! /usr/bin/python3

from matplotlib import pyplot as plt
from astropy.table import Table
import numpy as np

hipparcos_file = Table.read(f"/usr/share/dr3map/hipparcos/I_239_selection.csv", format='ascii')

# Create HRD plot of the double stars based on Hipparcos
hipparcos_abs_mag = hipparcos_file['Abs_mag']
hipparcos_bv_index = hipparcos_file['B-V']

colors = (hipparcos_bv_index)
plt.scatter(hipparcos_bv_index, hipparcos_abs_mag, c=colors, s=0.5, alpha=0.1, cmap='RdYlBu_r', vmax=1.9, vmin=-0.4) #, 
#plt.scatter(bv_a, mag_abs_a, s=14, color="blue", label='Main star') # s= 1 / mag_abs_a
#plt.scatter(bv_b, mag_abs_b, s=7, color="red", label='Companion star') # s= 1 / mag_abs_a
plt.legend(loc="upper left")
plt.axis((-0.4,1.9,15,-10))
plt.title('Double Star H-R Diagram')
plt.xlabel('B-V index')
plt.ylabel('Absolute magnitude')
plt.gca().set_aspect(0.1)
#savename = str(workingDirectory + '/' + pairname + '_hrd.jpg').replace(' ', '')
#plt.savefig(savename, bbox_inches='tight', dpi=150.0)
plt.show()
plt.close()

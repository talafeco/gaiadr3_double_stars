import sys
from matplotlib import pyplot as plt
from astropy.table import Table

hipparcos_file = Table.read(sys.argv[1], format='ascii')

abs_mag = [5, 7]
bv_index = [1, 1.2]

hipparcos_abs_mag = hipparcos_file['Abs_mag']
hipparcos_bv_index = hipparcos_file['B-V']

#plt.scatter(bv_index, abs_mag)
plt.scatter(hipparcos_bv_index, hipparcos_abs_mag, s=0.5, alpha=0.015, color="grey") #, 
plt.scatter(bv_index[0], abs_mag[0], s=6, color="blue", label='Main star')
plt.scatter(bv_index[1], abs_mag[1], s=3, color="red", label='Companion star')
plt.legend(loc="upper left")
plt.axis((-0.6,2.1,21,-16))
plt.title('Double Star Hertzsprung-Russell Diagram')
plt.xlabel('B-V index')
plt.ylabel('Absolute magnitude')

plt.savefig('ds_hrd.png', bbox_inches='tight')

plt.show()
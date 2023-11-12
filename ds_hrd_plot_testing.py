import sys
from matplotlib import pyplot as plt
from astropy.table import Table

pairname = 'HJ 3060'
gaia_file = Table.read(sys.argv[1], format='ascii')

abs_mag = [4.489298081260197, 2.746070876244115]
bv_index = [0.30146027, 0.27097988]

hipparcos_abs_mag = gaia_file['Abs_mag']
hipparcos_bv_index = gaia_file['B-V']

#plt.scatter(bv_index, abs_mag)
plt.scatter(hipparcos_bv_index, hipparcos_abs_mag, s=0.5, alpha=0.2, color="grey") #, 
plt.scatter(bv_index[0], abs_mag[0], s= 1 / abs_mag[0] * 40, color="blue", label='Main star')
plt.scatter(bv_index[1], abs_mag[1], s= 1 / abs_mag[1] * 40, color="red", label='Companion star')
plt.legend(loc="upper left")
plt.axis((-0.4,1.9,21,-16))
plt.title('Double Star ' + pairname + ' H-R Diagram')
plt.xlabel('B-V index')
plt.ylabel('Absolute magnitude')

plt.gca().set_aspect(0.07)

plt.savefig('ds_hrd.png', bbox_inches='tight')

plt.show()




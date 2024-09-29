#! /usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from astroquery.gaia import Gaia
from astropy.table import Table
import pandas as pd

gaia_file = Table.read(f'/usr/share/dr3map/gaia/star_catalog.csv', format='ascii')

df = gaia_file.to_pandas()

colors = (df['bp_rp'])
plt.figure(figsize=(10, 10), facecolor='black')  # Set figure background to black
ax = plt.gca()
ax.set_facecolor('black')  # Set axes background to black
ax.spines['bottom'].set_color('white')
ax.spines['top'].set_color('white') 
ax.spines['right'].set_color('white')
ax.spines['left'].set_color('white')
plt.scatter(df['bp_rp'], df['abs_mag'], c=colors, s=0.5, alpha=0.2, cmap='RdYlBu_r') #, 
#plt.scatter(bv_a, mag_abs_a, s=14, color="blue", label='Main star') # s= 1 / mag_abs_a
#plt.scatter(bv_b, mag_abs_b, s=7, color="red", label='Companion star') # s= 1 / mag_abs_a
plt.axis((-0.5,4,12,-2))
#plt.title('Double Star ' + pairname + ' H-R Diagram')
plt.xlabel('G_BP - G_RP index', color='white')
plt.ylabel('Absolute magnitude', color='white')
plt.title('Hertzsprung-Russell Diagram (Color Index vs Absolute Magnitude)', color='white')
#plt.grid(color='gray', linestyle='--', linewidth=0.5)
ax.tick_params(axis='x', colors='white')
ax.tick_params(axis='y', colors='white')
savename = str('_hrd.jpg').replace(' ', '')
plt.savefig(savename, bbox_inches='tight', dpi=300.0)
plt.show()
plt.close()
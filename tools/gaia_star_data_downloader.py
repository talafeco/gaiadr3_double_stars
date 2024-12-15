#! /usr/bin/python3

import numpy as np
from astroquery.gaia import Gaia
from astropy.table import Table
import pandas as pd

# Step 1: Query Gaia DR3 using astroquery to get effective temperature from the database

job = Gaia.launch_job_async("""
SELECT TOP 1000000
    source_id, 
    parallax, 
    parallax_error, 
    phot_g_mean_mag, 
    phot_bp_mean_mag, 
    phot_rp_mean_mag, 
    teff_gspphot
FROM gaiadr3.gaia_source
WHERE parallax > 0 
      AND parallax/parallax_error > 10
      AND phot_g_mean_mag IS NOT NULL
      AND phot_bp_mean_mag IS NOT NULL
      AND phot_rp_mean_mag IS NOT NULL
      AND teff_gspphot IS NOT NULL
ORDER BY parallax DESC
""")

gaia_data = job.get_results()


# Step 2: Convert to Pandas DataFrame for easier manipulation
df = gaia_data.to_pandas()

# Step 3: Calculate the absolute magnitude (M_G)
df['distance_pc'] = 1000.0 / df['parallax']  # Distance in parsecs
df['abs_mag'] = df['phot_g_mean_mag'] - 5 * np.log10(df['distance_pc']) + 5

# Step 4: Calculate the color index (G_BP - G_RP)
df['bp_rp'] = df['phot_bp_mean_mag'] - df['phot_rp_mean_mag']
df['bp_g'] = df['phot_bp_mean_mag'] - df['phot_g_mean_mag']

df.to_csv('star_catalog_1m.csv', index=False)

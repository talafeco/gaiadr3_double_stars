import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from astropy.coordinates import SkyCoord
import numpy as np
from matplotlib.colors import TwoSlopeNorm

def plot_3d_stellar_positions_from_csv(filename):
    """
    Reads Right Ascension (RA) in HMS, Declination (Dec) in DMS, and Distance data from a CSV file
    and generates a 3D plot with Distance replacing Declination on the Y-axis.
    
    Parameters:
    - filename: Name of the CSV file containing 'ra_hms', 'dec_dms', and 'distance' columns.
    
    Returns:
    - A 3D plot of the data
    """
    # Read the CSV file
    try:
        data = pd.read_csv(filename)
    except Exception as e:
        raise FileNotFoundError(f"Error reading file '{filename}': {e}")

    # Extract the data
    ra = data['ra']
    dec = data['dec']
    distance = 1000 / data['parallax']
    magnitude = data['phot_g_mean_mag']
    bp_rp = data['bp_rp']

    sky_coords = SkyCoord(ra=ra, dec=dec, unit=('deg', 'deg'))

    # Get RA and Dec in decimal degrees
    ra_decimal = sky_coords.ra.deg
    dec_decimal = sky_coords.dec.deg

    # Plot with BP-RP as the color
    sc = plt.scatter(ra_decimal, dec_decimal, s=(15 - magnitude)*15, c=bp_rp, cmap='RdYlBu_r') # 
    plt.colorbar(sc, shrink=0.5, aspect=10, label='Blue - Red color index')

    # Labeling the axes
    plt.xlabel('RA (deg)')
    plt.ylabel('DEC (deg)')
    plt.title('Kernya 77')
    plt.show()

# Example Usage
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python script_name.py <csv_filename>")
    else:
        csv_filename = sys.argv[1]
        plot_3d_stellar_positions_from_csv(csv_filename)
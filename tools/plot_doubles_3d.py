import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from astropy.coordinates import SkyCoord
import numpy as np

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
    
    # Ensure required columns exist
    required_columns = {'ra_hms', 'dec_dms', 'distance_pc'}
    if not required_columns.issubset(data.columns):
        raise ValueError(f"The file must contain the following columns: {required_columns}")

    # Extract the data
    ra_hms = data['ra_hms']
    dec_dms = data['dec_dms']
    distance = data['distance_pc']

    # Convert RA and Dec from HMS/DMS to decimal degrees using Astropy
    try:
        sky_coords = SkyCoord(ra=ra_hms, dec=dec_dms, unit=('hourangle', 'deg'))
    except Exception as e:
        raise ValueError(f"Error converting coordinates: {e}")
    
    # Get RA and Dec in decimal degrees
    ra_decimal = sky_coords.ra.deg
    dec_decimal = sky_coords.dec.deg

    # Convert RA and Dec to Cartesian coordinates for the plot
    x = ra_decimal  # RA in decimal degrees
    y = distance  # Distance in parsecs
    z = dec_decimal  # Dec in decimal degrees

    # Create 3D figure
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Scatter plot
    scatter = ax.scatter(x, y, z, c=distance, cmap='viridis', s=50)
    colorbar = fig.colorbar(scatter, ax=ax, shrink=0.5, aspect=10)
    colorbar.set_label('Distance (pc)')

    # Labeling the axes
    ax.set_xlabel('Right Ascension (RA) [°]')
    ax.set_ylabel('Distance [pc]')
    ax.set_zlabel('Declination (Dec) [°]')
    ax.set_title('3D Stellar Positions')
    plt.show()

# Example Usage
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python script_name.py <csv_filename>")
    else:
        csv_filename = sys.argv[1]
        plot_3d_stellar_positions_from_csv(csv_filename)
import sys
import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
from astropy.coordinates import Longitude

# Function to read CSV file and extract relevant data
def read_star_data(filename):
    data = pd.read_csv(filename)
    ra = data['ra'].values * u.deg
    dec = data['dec'].values * u.deg
    parallax = data['parallax'].values * u.mas
    pmra = data['pmra'].values * u.mas / u.yr
    pmdec = data['pmdec'].values * u.mas / u.yr
    designation = data['designation'].values
    return designation, ra, dec, parallax, pmra, pmdec

# Function to calculate star positions for a given time step
def calculate_positions(ra, dec, parallax, pmra, pmdec, time_step):
    # Convert parallax to distance in parsecs
    distance_pc = 1000 / parallax

    # Calculate position offsets from proper motion
    delta_ra = (pmra * time_step).to(u.deg) / np.cos(dec)
    delta_dec = (pmdec * time_step).to(u.deg)

    # Update positions
    updated_ra = ra - delta_ra
    updated_dec = dec - delta_dec

    # Clamp declination to the valid range [-90, 90] degrees
    updated_dec = np.clip(updated_dec, -90 * u.deg, 90 * u.deg)

    # Handle RA adjustment if Dec is clamped
    updated_ra = Longitude(updated_ra)  # Convert to Longitude for wrap_at
    updated_ra = updated_ra.wrap_at(360 * u.deg)

    # Create SkyCoord object with updated positions
    coords = SkyCoord(ra=updated_ra, dec=updated_dec, distance=distance_pc * u.pc, frame='icrs')
    return coords

# Function to create and save a 3D plot
def create_3d_plot(ra, dec, distance_pc, iteration, output_dir):
    # Extract values for plotting
    if hasattr(ra, "value"):
        x = ra.value  # RA in decimal degrees
        y = distance_pc.value  # Distance in parsecs
        z = dec.value  # Dec in decimal degrees
    else:
        x = ra  # RA in decimal degrees
        y = distance_pc  # Distance in parsecs
        z = dec  # Dec in decimal degrees

    # Create 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, s=1, alpha=0.7)

    # Label axes
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('Distance (pc)')
    ax.set_zlabel('Dec (deg)')
    ax.set_title(f"Star Positions at {iteration * -100000} years")

    # Save the plot
    os.makedirs(output_dir, exist_ok=True)
    plot_file = os.path.join(output_dir, f"star_positions_{iteration * -100000}_years.png")
    plt.savefig(plot_file)
    plt.close()
    print(f"Saved 3D plot for {iteration * -100000} years to {plot_file}")

# Main program
if len(sys.argv) < 3:
    print("Usage: python program.py <input_csv_file> <number_of_iterations>")
    sys.exit(1)

filename = sys.argv[1]
iterations = int(sys.argv[2])

# Read star data
designation, ra, dec, parallax, pmra, pmdec = read_star_data(filename)

# Output directory for plots
output_dir = "star_position_plots"

# Time step for each iteration
time_step = 1e6 * u.yr  # 100,000 years

# Iterate and create plots
for iteration in range(1, iterations + 1):
    print(f"Iteration {iteration}: Calculating positions and creating plot...")
    
    # Calculate positions
    coords = calculate_positions(ra, dec, parallax, pmra, pmdec, time_step * iteration)
    
    # Create and save 3D plot
    create_3d_plot(coords.ra, coords.dec, coords.distance, iteration, output_dir)

print(f"All plots saved in {output_dir}.")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import Normalize

def plot_star_data(csv_file):
    # Read the data from the CSV file
    data = pd.read_csv(csv_file)

    # Combine the mass, absolute magnitude, and luminosity columns
    mass = pd.concat([data['mass_a'], data['mass_b']]).values
    abs_mag = pd.concat([data['abs_mag_a'], data['abs_mag_b']]).values
    luminosity = pd.concat([data['lum_a'], data['lum_b']]).values

    # Normalize the mass for color mapping
    norm = Normalize(vmin=np.min(mass), vmax=np.max(mass))
    colors = cm.jet(norm(mass))

    # Create the 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    scatter = ax.scatter(mass, abs_mag, luminosity, c=colors, marker='o')

    # Add color bar to indicate mass values
    sm = plt.cm.ScalarMappable(cmap=cm.jet, norm=norm)
    sm.set_array([])  # Dummy array for the color bar
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Mass (Solar Masses)')

    # Set axis labels
    ax.set_xlabel('Mass (Solar Masses)')
    ax.set_ylabel('Absolute Magnitude (Gaia G band)')
    ax.set_zlabel('Luminosity (Solar Luminosities)')

    plt.title('3D Star Data Plot')
    plt.show()

# Example usage:
# plot_star_data('path_to_your_file.csv')


# Example Usage
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python script_name.py <csv_filename>")
    else:
        csv_filename = sys.argv[1]
        plot_star_data(csv_filename)
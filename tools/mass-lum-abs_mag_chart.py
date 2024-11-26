import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

def plot_star_mass(csv_file):
    """Plots the mass of the stars, colored by mass, with an average trend line."""
    data = pd.read_csv(csv_file)
    mass = np.concatenate([data['mass_a'].values, data['mass_b'].values])
    avg_mass = mass.mean()

    norm = Normalize(vmin=mass.min(), vmax=mass.max())
    cmap = plt.cm.coolwarm

    scatter = plt.scatter(range(len(mass)), mass, c=mass, cmap=cmap, norm=norm)
    plt.axhline(y=avg_mass, color='black', linestyle='--', linewidth=2, label=f'Average: {avg_mass:.2f}')
    plt.colorbar(scatter, label='Mass (Solar Masses)')
    plt.xlabel('Star Index')
    plt.ylabel('Mass (Solar Masses)')
    plt.title('Mass of Stars')
    plt.legend()
    plt.show()

def plot_mass_luminosity(csv_file):
    """Plots the mass-luminosity relationship."""
    data = pd.read_csv(csv_file)
    mass = data[['mass_a', 'mass_b']].mean(axis=1)
    luminosity = data[['lum_a', 'lum_b']].mean(axis=1)

    plt.scatter(mass, luminosity, c=mass, cmap=plt.cm.coolwarm, norm=Normalize(vmin=mass.min(), vmax=mass.max()))
    plt.colorbar(label='Mass (Solar Masses)')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Mass (Solar Masses)')
    plt.ylabel('Luminosity (Solar Luminosities)')
    plt.title('Mass-Luminosity Relationship')
    plt.show()

def plot_hr_diagram(csv_file):
    """Plots the Hertzsprung-Russell diagram."""
    data = pd.read_csv(csv_file)
    abs_mag = np.concatenate([data['abs_mag_a'].values, data['abs_mag_b'].values])
    luminosity = np.concatenate([data['lum_a'].values, data['lum_b'].values])

    plt.scatter(abs_mag, luminosity, c=luminosity, cmap=plt.cm.coolwarm.reversed(), norm=Normalize(vmin=luminosity.min(), vmax=luminosity.max()))
    plt.colorbar(label='Luminosity (Solar Luminosities)')
    #plt.gca().invert_xaxis()  # HR diagrams often have inverted magnitude axes
    plt.gca()  # HR diagrams often have inverted magnitude axes
    plt.yscale('log')
    plt.xlabel('Absolute Magnitude')
    plt.ylabel('Luminosity (Solar Luminosities based on Gaia G band)')
    plt.title('Hertzsprung-Russell Diagram')
    plt.show()

def plot_luminosity_distribution(csv_file):
    """Plots the distribution of stellar luminosities."""
    data = pd.read_csv(csv_file)
    luminosity = np.concatenate([data['lum_a'].values, data['lum_b'].values])

    plt.hist(luminosity, bins=30, log=True, color='blue', edgecolor='black', alpha=0.7)
    plt.xlabel('Luminosity (Solar Luminosities based on Gaia G band)')
    plt.ylabel('Number of Stars')
    plt.title('Luminosity Distribution')
    plt.show()

def plot_mass_distribution(csv_file):
    """Plots the distribution of stellar masses."""
    data = pd.read_csv(csv_file)
    mass = np.concatenate([data['mass_a'].values, data['mass_b'].values])

    plt.hist(mass, bins=30, color='green', edgecolor='black', alpha=0.7)
    plt.xlabel('Mass (Solar Masses)')
    plt.ylabel('Number of Stars')
    plt.title('Mass Distribution')
    plt.show()

def plot_absolute_magnitude_luminosity(csv_file):
    """Plots the absolute magnitude vs luminosity."""
    data = pd.read_csv(csv_file)
    abs_mag = np.concatenate([data['abs_mag_a'].values, data['abs_mag_b'].values])
    luminosity = np.concatenate([data['lum_a'].values, data['lum_b'].values])

    plt.scatter(abs_mag, luminosity, c=luminosity, cmap=plt.cm.coolwarm.reversed(), norm=Normalize(vmin=luminosity.min(), vmax=luminosity.max()))
    plt.colorbar(label='Luminosity (Solar Luminosities)')
    plt.gca().invert_xaxis()  # Inverted for absolute magnitude
    #plt.gca()  # Inverted for absolute magnitude
    plt.yscale('log')
    plt.xlabel('Absolute Magnitude')
    plt.ylabel('Luminosity (Solar Luminosities)')
    plt.title('Absolute Magnitude vs Luminosity')
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Generate star data plots from a CSV file.")
    parser.add_argument('csv_file', type=str, help='Path to the CSV file containing star data.')
    parser.add_argument('--plot', type=str, required=True, choices=[
        'mass', 'mass_luminosity', 'hr_diagram', 'luminosity_distribution',
        'mass_distribution', 'absolute_magnitude_luminosity'
    ], help="Choose the type of plot to generate.")

    args = parser.parse_args()

    if args.plot == 'mass':
        plot_star_mass(args.csv_file)
    elif args.plot == 'mass_luminosity':
        plot_mass_luminosity(args.csv_file)
    elif args.plot == 'hr_diagram':
        plot_hr_diagram(args.csv_file)
    elif args.plot == 'luminosity_distribution':
        plot_luminosity_distribution(args.csv_file)
    elif args.plot == 'mass_distribution':
        plot_mass_distribution(args.csv_file)
    elif args.plot == 'absolute_magnitude_luminosity':
        plot_absolute_magnitude_luminosity(args.csv_file)
    else:
        print("Invalid plot type selected.")

if __name__ == '__main__':
    main()

'''For spectral type distribution: python star_plots.py stars.csv spectral_type_distribution
For HR diagram: python star_plots.py stars.csv hr_diagram'''

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def calculate_spectral_type(teff):
    """
    Calculate spectral type based on effective temperature (Teff).
    Spectral classification is simplified and mapped to temperature ranges.
    """
    if teff >= 30000:
        return 'O'
    elif teff >= 10000:
        return 'B'
    elif teff >= 7500:
        return 'A'
    elif teff >= 6000:
        return 'F'
    elif teff >= 5200:
        return 'G'
    elif teff >= 3700:
        return 'K'
    elif teff >= 1700:
        return 'M'

def plot_spectral_type_distribution(df):
    """
    Plot the distribution of spectral types based on Teff values.
    """
    df['Spectral Type'] = df['teff_gspphot'].apply(calculate_spectral_type)
    spectral_type_counts = df['Spectral Type'].value_counts().sort_index()

    plt.figure(figsize=(10, 6))
    spectral_type_counts.plot(kind='bar', color='skyblue')
    plt.title('Spectral Type Distribution')
    plt.xlabel('Spectral Type')
    plt.ylabel('Number of Stars')
    plt.xticks(rotation=0)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()

def calculate_absolute_magnitude(parallax, phot_g_mean_mag):
    """
    Calculate absolute magnitude from parallax and apparent magnitude.
    """
    # Avoid division by zero in parallax
    parallax = np.where(parallax > 0, parallax, np.nan)
    distance = 1000 / parallax  # Distance in parsecs
    abs_magnitude = phot_g_mean_mag - 5 * np.log10(distance) + 5
    return abs_magnitude

def plot_hr_diagram(df):
    spectral_boundaries = [30000, 10000, 7500, 6000, 5200, 3700, 2400]
    spectral_types = ['O', 'B', 'A', 'F', 'G', 'K', 'M']
    
    df['Spectral Type'] = df['teff_gspphot'].apply(calculate_spectral_type)
    df['Absolute Magnitude'] = calculate_absolute_magnitude(df['parallax'], df['phot_g_mean_mag'])

    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Set the background color to black
    ax.set_facecolor('black')
    fig.patch.set_facecolor('black')
    
    colors = ['blue', 'cyan', 'white', 'yellow', 'orange', 'red', 'darkred']
    
    for spectral_type, color in zip(spectral_types, colors):
        subset = df[df['Spectral Type'] == spectral_type]
        ax.scatter(subset['teff_gspphot'], subset['Absolute Magnitude'], 
                   label=spectral_type, color=color, s=10, alpha=0.7)

    ax.invert_yaxis()  # HR diagrams have brighter stars at the top
    ax.invert_xaxis()  # Higher temperature stars on the left
    ax.set_title('Hertzsprung-Russell Diagram', color='white')
    ax.set_xlabel('Effective Temperature (K)', color='white')
    ax.set_ylabel('Absolute Magnitude (Gaia G band)', color='white')
    ax.legend(title="Spectral Type", facecolor='black', edgecolor='white', labelcolor='white')
    ax.grid(alpha=0.3, color='gray')  # Adjust grid color for visibility
    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')

    # Add a secondary x-axis for spectral types
    def teff_to_spectral_type(teff):
        """
        Convert effective temperature to a spectral type index for labeling.
        """
        for i, boundary in enumerate(spectral_boundaries):
            if teff >= boundary:
                return spectral_types[i]
        return 'M'  # Default to 'M' for low temperatures

    def spectral_type_formatter(teff, pos):
        """Formatter for the spectral type secondary axis."""
        return teff_to_spectral_type(teff)

    secax = ax.secondary_xaxis('top')
    secax.set_xlabel('Spectral Type', color='white')
    secax.set_xticks(spectral_boundaries, labels=spectral_types)
    secax.tick_params(axis='x', colors='white')

    plt.tight_layout()
    plt.show()

def main():
    if len(sys.argv) < 3:
        print("Usage: python program.py <csv_file> <plot_type>")
        print("Plot types: spectral_type_distribution, hr_diagram")
        sys.exit(1)

    csv_file = sys.argv[1]
    plot_type = sys.argv[2]

    try:
        df = pd.read_csv(csv_file)
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        sys.exit(1)

    if plot_type == "spectral_type_distribution":
        plot_spectral_type_distribution(df)
    elif plot_type == "hr_diagram":
        plot_hr_diagram(df)
    else:
        print(f"Unknown plot type: {plot_type}")
        print("Available plot types: spectral_type_distribution, hr_diagram")

if __name__ == "__main__":
    main()

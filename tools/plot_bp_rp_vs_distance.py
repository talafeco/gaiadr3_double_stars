import pandas as pd
import matplotlib.pyplot as plt

def plot_bp_rp_distance_from_csv(filename):
    """
    Plots Gaia BP-RP magnitude index against distance, coloring points based on BP-RP values.
    
    Parameters:
    filename (str): Path to the CSV file containing 'bp_rp' and 'parallax' columns.
    """
    try:
        # Load the data from CSV
        data = pd.read_csv(filename)
    except Exception as e:
        raise FileNotFoundError(f"Error reading file '{filename}': {e}")
    
    # Ensure required columns exist
    if 'bp_rp' not in data or 'parallax' not in data:
        raise ValueError("The input CSV file must contain 'bp_rp' and 'parallax' columns.")

    # Extract the data
    bp_rp = data['bp_rp']
    # Avoid division by zero
    data['parallax'] = data['parallax'].replace(0, float('nan'))
    distance = 1000 / data['parallax']

    # Plot with BP-RP as the color
    sc = plt.scatter(distance, bp_rp, c=bp_rp, cmap='RdYlBu_r', s=50, edgecolor='k', alpha=0.7)

    # Labeling the axes and the plot
    plt.xlabel('Distance [pc]')
    plt.ylabel('Gaia BP-RP magnitude index')
    plt.title('Gaia BP-RP Index vs Distance')
    plt.colorbar(sc, label='BP-RP Index')  # Add colorbar for BP-RP values
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.show()

# Example Usage
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python script_name.py <csv_filename>")
    else:
        csv_filename = sys.argv[1]
        plot_bp_rp_distance_from_csv(csv_filename)

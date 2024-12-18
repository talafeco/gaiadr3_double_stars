import numpy as np
import pandas as pd
import sys
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.spatial import distance
from collections import defaultdict
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt

def load_data(input_file):
    """
    Load and validate the CSV file containing stellar data.
    """
    try:
        data = pd.read_csv(input_file)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    # Check required columns
    required_columns = ['designation', 'ra', 'dec', 'parallax', 'pmra', 'pmdec', 'radial_velocity']
    for col in required_columns:
        if col not in data.columns:
            print(f"Error: Missing column '{col}'")
            sys.exit(1)

    # Drop rows with missing required data
    data = data.dropna(subset=required_columns)
    return data

def convert_to_skycoord(data):
    """
    Convert input data to Astropy SkyCoord objects.
    """
    ra = data['ra'].values * u.deg
    dec = data['dec'].values * u.deg
    parallax = data['parallax'].values * u.mas
    pmra = data['pmra'].values * u.mas / u.yr
    pmdec = data['pmdec'].values * u.mas / u.yr
    radial_velocity = data['radial_velocity'].values * u.km / u.s

    # Convert parallax to distance (parallax > 0 to avoid invalid distances)
    distance = 1000.0 / parallax.value  # Convert to parsecs
    distance = distance * u.pc

    # Create SkyCoord object
    coords = SkyCoord(ra=ra, dec=dec, distance=distance,
                      pm_ra_cosdec=pmra, pm_dec=pmdec,
                      radial_velocity=radial_velocity, frame='icrs')

    return coords, data['designation'].values

def trace_back_positions(coords, time):
    """
    Trace back the positions of stars based on velocities.
    """
    velocities = coords.velocity.d_xyz.to(u.lyr / u.yr)
    positions = coords.cartesian.xyz.to(u.lyr)
    origin_positions = positions - velocities * time
    return positions, velocities, origin_positions

def calculate_distances(positions, indices):
    """
    Calculate minimal pairwise distance within a cluster.
    """
    cluster_positions = positions[:, indices].T.value  # Strip units and transpose
    dists = distance.pdist(cluster_positions, metric='euclidean')
    return np.min(dists) if len(dists) > 0 else 0.0

def analyze_clusters(positions, velocities, designations, labels):
    """
    Analyze clusters, compute convergence times, and report results.
    """
    clusters = defaultdict(list)
    for i, label in enumerate(labels):
        if label != -1:
            clusters[label].append(i)

    print(f"Found {len(clusters)} clusters of stars.")
    for cluster_id, indices in clusters.items():
        indices = np.array(indices)
        if np.any(indices >= positions.shape[1]):
            print(f"Error: Invalid indices in cluster {cluster_id}")
            continue

        # Calculate average positions and velocities
        cluster_positions = positions[:, indices]
        cluster_velocities = velocities[:, indices]
        avg_position = np.mean(cluster_positions, axis=1)
        avg_velocity = np.mean(cluster_velocities, axis=1)

        # Compute convergence time
        t_convergence = -np.dot(avg_position, avg_velocity) / np.dot(avg_velocity, avg_velocity)

        # Get star designations
        cluster_designations = designations[indices]

        # Calculate distances
        min_dist_at_convergence = calculate_distances(cluster_positions, range(len(indices)))
        current_positions = positions[:, indices].T.value
        current_min_dist = np.min(distance.pdist(current_positions, metric='euclidean'))

        # Report results
        print(f"\nCluster {cluster_id}:")
        print(f" - Number of stars: {len(indices)}")
        print(f" - Star designations: {', '.join(cluster_designations)}")
        if t_convergence > 0:
            print(f" - Convergence time: {t_convergence:.2f} years ago")
        else:
            print(" - Convergence time: Not valid (not in the past).")
        print(f" - Minimal distance at convergence: {min_dist_at_convergence:.4f} light-years")
        print(f" - Current minimal separation: {current_min_dist:.4f} light-years")

def plot_origin_positions(origin_positions):
    """
    Plot the traced-back positions of stars.
    """
    plt.figure(figsize=(8, 6))
    plt.scatter(origin_positions[0].value, origin_positions[1].value, s=5, color='blue', alpha=0.5)
    plt.xlabel("X Position (light-years)")
    plt.ylabel("Y Position (light-years)")
    plt.title("Traced-back Star Positions (Origin)")
    plt.grid(True)
    plt.show()

def plot_current_positions(positions):
    """
    Plot the current positions of stars.
    """
    plt.figure(figsize=(8, 6))
    plt.scatter(positions[0].value, positions[1].value, s=5, color='red', alpha=0.5)
    plt.xlabel("X Position (light-years)")
    plt.ylabel("Y Position (light-years)")
    plt.title("Current Star Positions")
    plt.grid(True)
    plt.show()

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_file.csv>")
        sys.exit(1)

    input_file = sys.argv[1]
    print(f"Loading data from {input_file}...")
    data = load_data(input_file)

    print("Converting to SkyCoord objects...")
    coords, designations = convert_to_skycoord(data)

    # Parameters
    time = 1e5 * u.yr  # 100,000 years back
    eps = 10.0         # Clustering radius in light-years

    print(f"Tracing back positions over {time}...")
    positions, velocities, origin_positions = trace_back_positions(coords, time)

    print("Clustering stars...")
    origin_array = np.column_stack([origin_positions[0].value,
                                    origin_positions[1].value,
                                    origin_positions[2].value])
    dbscan = DBSCAN(eps=eps, min_samples=3, metric='euclidean')
    labels = dbscan.fit_predict(origin_array)

    analyze_clusters(positions, velocities, designations, labels)

    # Plotting
    print("Plotting traced-back positions...")
    plot_origin_positions(origin_positions)

    print("Plotting current positions...")
    plot_current_positions(positions)

if __name__ == "__main__":
    main()

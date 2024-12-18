import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.spatial import distance

# Example Data: Replace these with your actual inputs
num_stars = 10000
ra = np.random.uniform(0, 360, num_stars) * u.deg          # Right Ascension
dec = np.random.uniform(-90, 90, num_stars) * u.deg       # Declination
distance_pc = np.random.uniform(10, 1000, num_stars) * u.pc  # Distance in parsecs
pm_ra_cosdec = np.random.uniform(-10, 10, num_stars) * u.mas/u.yr  # Proper Motion RA
pm_dec = np.random.uniform(-10, 10, num_stars) * u.mas/u.yr  # Proper Motion Dec
radial_velocity = np.random.uniform(-100, 100, num_stars) * u.km/u.s  # Radial Velocity

# Step 1: Convert input data to SkyCoord objects
coords = SkyCoord(ra=ra, dec=dec, distance=distance_pc, 
                  pm_ra_cosdec=pm_ra_cosdec, pm_dec=pm_dec, 
                  radial_velocity=radial_velocity, frame='icrs')

# Step 2: Calculate velocities and trace back positions
velocities = coords.velocity.d_xyz.to(u.lyr/u.yr)  # Convert to light-years/year
positions = coords.cartesian.xyz.to(u.lyr)  # Positions in light-years

# Trace back trajectories
time = 1e6 * u.yr  # Large backward time
origin_positions = positions - velocities * time  # Projected origins

# Step 3: Calculate pairwise distances between origins
origin_array = np.column_stack([origin_positions[0].value,
                                origin_positions[1].value,
                                origin_positions[2].value])

# Compute pairwise distances
pairwise_dist = distance.pdist(origin_array, metric='euclidean')
distance_matrix = distance.squareform(pairwise_dist)

# Step 4: Find stars within 2 light-years of each other
threshold_distance = 2.0  # 2 light-years
close_pairs = np.argwhere((distance_matrix < threshold_distance) & (distance_matrix > 0))

# Step 5: Group stars into clusters manually
def group_clusters(close_pairs, num_stars):
    visited = set()
    clusters = []

    def dfs(node, cluster):
        if node not in visited:
            visited.add(node)
            cluster.append(node)
            neighbors = np.unique(close_pairs[np.where(close_pairs[:, 0] == node)][:, 1])
            for neighbor in neighbors:
                dfs(neighbor, cluster)

    for i in range(num_stars):
        if i not in visited:
            cluster = []
            dfs(i, cluster)
            if len(cluster) > 1:  # Only consider clusters with > 1 star
                clusters.append(cluster)
    return clusters

clusters = group_clusters(close_pairs, num_stars)

# Step 6: Print the results
print(f"Number of clusters found: {len(clusters)}")
for idx, cluster in enumerate(clusters[:10]):  # Print first 10 clusters
    print(f"Cluster {idx}: {len(cluster)} stars")
    print(f"Indices: {cluster}")

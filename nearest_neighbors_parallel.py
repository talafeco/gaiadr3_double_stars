import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from joblib import Parallel, delayed
import math

def find_nearest_chunk(tree, chunk, offset):
    distances, indices = tree.query(chunk)
    matched_coords = tree.data[indices]
    return pd.DataFrame({
        "original_x": chunk[:, 0],
        "original_y": chunk[:, 1],
        "matched_x": matched_coords[:, 0],
        "matched_y": matched_coords[:, 1],
        "index": indices,
        "distance": distances
    })

def find_nearest_neighbors_parallel(large_coords_flat, small_coords_flat, n_jobs=4):
    large_coords = np.array(large_coords_flat, dtype=np.float64).reshape(-1, 2)
    small_coords = np.array(small_coords_flat, dtype=np.float64).reshape(-1, 2)
    tree = cKDTree(large_coords)
    chunk_size = math.ceil(len(small_coords) / n_jobs)
    chunks = [small_coords[i:i + chunk_size] for i in range(0, len(small_coords), chunk_size)]
    results = Parallel(n_jobs=n_jobs)(delayed(find_nearest_chunk)(tree, chunk, i * chunk_size) 
                                      for i, chunk in enumerate(chunks))
    return pd.concat(results, ignore_index=True)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Find nearest neighbors using KDTree and multiprocessing.")
    parser.add_argument("--large_input", required=True, help="Path to binary .npy file with large coordinate set")
    parser.add_argument("--small_input", required=True, help="Path to binary .npy file with small coordinate set")
    parser.add_argument("--output", required=True, help="Path to save output CSV file")
    parser.add_argument("--n_jobs", type=int, default=4, help="Number of parallel jobs")
    args = parser.parse_args()
    large_coords = np.load(args.large_input)
    small_coords = np.load(args.small_input)
    results_df = find_nearest_neighbors_parallel(large_coords, small_coords, n_jobs=args.n_jobs)
    results_df.to_csv(args.output, index=False)
    print(f"Results saved to {args.output}")

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from joblib import Parallel, delayed
import math

def find_nearest_chunk(tree, chunk, original_names):
    distances, indices = tree.query(chunk)
    matched_coords = tree.data[indices]
    return pd.DataFrame({
        "original_name": original_names,
        "original_ra": chunk[:, 0],
        "original_dec": chunk[:, 1],
        "matched_ra": matched_coords[:, 0],
        "matched_dec": matched_coords[:, 1],
        "index": indices,
        "distance": distances
    })

def find_nearest_neighbors_parallel(large_df, small_df, n_jobs=4):
    large_coords = large_df[["ra", "dec"]].to_numpy(dtype=np.float64)
    small_coords = small_df[["ra", "dec"]].to_numpy(dtype=np.float64)
    small_names = small_df["name"].to_numpy()

    tree = cKDTree(large_coords)
    chunk_size = math.ceil(len(small_coords) / n_jobs)
    chunks = [(small_coords[i:i + chunk_size], small_names[i:i + chunk_size]) 
              for i in range(0, len(small_coords), chunk_size)]

    results = Parallel(n_jobs=n_jobs)(
        delayed(find_nearest_chunk)(tree, chunk_coords, chunk_names)
        for chunk_coords, chunk_names in chunks
    )

    return pd.concat(results, ignore_index=True)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Find nearest neighbors from CSV files with RA/Dec coordinates.")
    parser.add_argument("--large_csv", required=True, help="Path to CSV file with large coordinate set (columns: name, ra, dec)")
    parser.add_argument("--small_csv", required=True, help="Path to CSV file with small coordinate set (columns: name, ra, dec)")
    parser.add_argument("--output", required=True, help="Path to save output CSV file")
    parser.add_argument("--n_jobs", type=int, default=4, help="Number of parallel jobs")
    args = parser.parse_args()

    large_df = pd.read_csv(args.large_csv)
    small_df = pd.read_csv(args.small_csv)

    results_df = find_nearest_neighbors_parallel(large_df, small_df, n_jobs=args.n_jobs)
    results_df.to_csv(args.output, index=False)
    print(f"Results saved to {args.output}")

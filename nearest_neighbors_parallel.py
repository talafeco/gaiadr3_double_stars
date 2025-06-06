
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from joblib import Parallel, delayed
import math

def normalize_columns(df):
    df.columns = df.columns.str.strip().str.lower()
    print("Normalized columns:", df.columns.tolist())  # Add this line for debug
    return df

# Angular distance in arcseconds
def angular_distance_arcsec(ra1, dec1, ra2, dec2):
    ra1_rad, dec1_rad = math.radians(ra1), math.radians(dec1)
    ra2_rad, dec2_rad = math.radians(ra2), math.radians(dec2)
    cos_angle = (math.sin(dec1_rad) * math.sin(dec2_rad) +
                 math.cos(dec1_rad) * math.cos(dec2_rad) * math.cos(ra1_rad - ra2_rad))
    cos_angle = min(1.0, max(-1.0, cos_angle))
    angle_rad = math.acos(cos_angle)
    return math.degrees(angle_rad) * 3600

def find_nearest_chunk(tree, chunk, original_names, large_names, large_meta_df, indices_offset):
    distances, indices = tree.query(chunk)
    matched_coords = tree.data[indices]
    matched_meta = large_meta_df.iloc[indices].reset_index(drop=True)
    angular_distances = [angular_distance_arcsec(r1, d1, r2, d2)
                         for (r1, d1), (r2, d2) in zip(chunk, matched_coords)]

    result_df = pd.DataFrame({
        "original_name": original_names,
        "original_ra": chunk[:, 0],
        "original_dec": chunk[:, 1],
        "wdss_name": matched_meta["name"].values,
        "wdss_ra": matched_coords[:, 0],
        "wdss_dec": matched_coords[:, 1],
        "index": indices,
        "distance_arcsec": angular_distances
    })

    for col in matched_meta.columns:
        if col not in ["name", "ra", "dec"]:
            result_df["wdss_" + col] = matched_meta[col].values

    return result_df

def find_nearest_neighbors_parallel(large_df, small_df, n_jobs=4):
    large_df = normalize_columns(large_df)
    small_df = normalize_columns(small_df)

    required_cols = {"name", "ra", "dec"}
    for label, df in [("large", large_df), ("small", small_df)]:
        if not required_cols.issubset(df.columns):
            raise ValueError(f"{label}_df must have columns: {required_cols}")

    large_coords = large_df[["ra", "dec"]].to_numpy(dtype=np.float64)
    small_coords = small_df[["ra", "dec"]].to_numpy(dtype=np.float64)
    small_names = small_df["name"].to_numpy()
    large_names = large_df["name"].to_numpy()

    tree = cKDTree(large_coords)
    chunk_size = math.ceil(len(small_coords) / n_jobs)
    chunks = [(small_coords[i:i + chunk_size], small_names[i:i + chunk_size])
              for i in range(0, len(small_coords), chunk_size)]

    results = Parallel(n_jobs=n_jobs)(
        delayed(find_nearest_chunk)(
            tree, chunk_coords, chunk_names, large_names, large_df, i * chunk_size
        )
        for i, (chunk_coords, chunk_names) in enumerate(chunks)
    )

    return pd.concat(results, ignore_index=True)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Find nearest neighbors from CSV files with RA/Dec coordinates.")
    parser.add_argument("--large_csv", required=True, help="Path to CSV file with large coordinate set (columns: name, ra, dec, ...)")
    parser.add_argument("--small_csv", required=True, help="Path to CSV file with small coordinate set (columns: name, ra, dec)")
    parser.add_argument("--output", required=True, help="Path to save output CSV file")
    parser.add_argument("--n_jobs", type=int, default=4, help="Number of parallel jobs")
    args = parser.parse_args()

    large_df = pd.read_csv(args.large_csv)
    small_df = pd.read_csv(args.small_csv)

    results_df = find_nearest_neighbors_parallel(large_df, small_df, n_jobs=args.n_jobs)
    results_df.to_csv(args.output, index=False)
    print(f"Results saved to {args.output}")

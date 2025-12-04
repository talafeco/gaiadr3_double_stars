# isotool.py
import os
import numpy as np
import pandas as pd
import math
import sys
from pathlib import Path
from calculator import compute_lum_from_MG0

# --- Import Custom Modules ---
from queries import (
    query_gaia_dr3_data,
    get_ebv_inf_from_dr3,
)
from calculator import (
    compute_MG_obs, compute_MG0,
    compute_BP_RP0, estimate_spectral_type_bprp,
    compute_logTeff, compute_distance_and_error, 
    compute_EBV, compute_galactic_xyz, compute_pm_total,
    compute_radius_from_logL_Teff, classify_hrd_position,
    get_teff
)
from writer import write_structured_output
# We import load_isochrones directly to optimize I/O
from isochrone_tool import find_best_isochrone, load_isochrones

# --- CONFIGURATION ---
ISO_PATH = Path(os.environ.get("ISO_DIR", "./izokrona_rendezett.txt"))

# Default fitting tolerances (formerly asked via input)
FIT_CONFIG = {
    "teff_tol": 0.05,   # +/- in log scale
    "logl_tol": 0.2,    # +/- in log scale
    "min_mass": 0.1     # Solar masses
}

# --- HELPER FUNCTIONS ---
def is_missing(val):
    return val is None or val == "--" or (hasattr(val, "mask") and val is np.ma.masked)

def safe_format(val, digits=3):
    try:
        if val is None or (hasattr(val, "mask") and val is np.ma.masked):
            return "--"
        val = float(val)
        return f"{val:.{digits}f}"
    except (ValueError, TypeError):
        return "--"

def print_vertical(table, title="Adatok"):
    """ Prints table data nicely without user interaction. """
    print("\n🔹 " + title + ":")
    if table is None or len(table) == 0:
        print("  --")
        return

    row = table[0]
    for col in table.colnames:
        # Skip error columns for clean display
        if any(x in col for x in ["_error", "err", "Error"]):
            continue
            
        value = row[col]
        if is_missing(value): value = '--'
        print(f"{col:>20}: {value}")

def main():
    print("🔭 IsoTool - Automated Stellar Parameter Engine")
    
    # --- 1. SINGLE USER INPUT ---
    gaia_id = input("Please enter the Gaia DR3 Source ID: ").strip()

    if not gaia_id.isdigit():
        print("❌ Error: Gaia ID must be numeric.")
        sys.exit(1)

    # --- 2. DATA ACQUISITION ---
    print(f"🔍 Querying Gaia DR3 for ID: {gaia_id}...")
    dr3_result = query_gaia_dr3_data(gaia_id)
    
    if dr3_result is None or len(dr3_result) == 0:
        print("❌ ID not found in Gaia DR3.")
        sys.exit(1)
        
    print_vertical(dr3_result, title=f"Gaia Data Found")

    # RUWE Check (Non-blocking warning)
    if "ruwe" in dr3_result.colnames:
        ruwe = dr3_result["ruwe"][0]
        if not is_missing(ruwe) and float(ruwe) > 1.4:
            print(f"⚠️  WARNING: RUWE = {ruwe:.2f} (> 1.4). Result may be unreliable (Binary?).")

    # --- 3. ASTROPHYSICAL CALCULATIONS ---
    print("\n⚙️  Calculating physical parameters...")
    
    # Extract basic data
    G_mag = dr3_result["phot_g_mean_mag"][0]
    bp_rp_obs = dr3_result["bp_rp"][0]
    plx = dr3_result["parallax"][0]
    plx_err = dr3_result["parallax_error"][0]
    
    # Kinematics
    distance, distance_err = compute_distance_and_error(plx, plx_err)
    
    ra = float(dr3_result["ra"][0])
    dec = float(dr3_result["dec"][0])
    l, b, x, y, z = compute_galactic_xyz(ra, dec, distance) if distance else (None,)*5
    
    pmra = dr3_result["pmra"][0] if "pmra" in dr3_result.colnames else None
    pmdec = dr3_result["pmdec"][0] if "pmdec" in dr3_result.colnames else None
    pm_total = compute_pm_total(pmra, pmdec)

    # Extinction & Intrinsic Props
    ebv_inf = get_ebv_inf_from_dr3(dr3_result)
    try:
        ebv_inf_val = float(ebv_inf)
    except (TypeError, ValueError):
        ebv_inf_val = None

    EBV = compute_EBV(ebv_inf_val, distance, b)
    Ag = 2.0 * EBV if EBV is not None else None
    
    MG_obs = compute_MG_obs(G_mag, distance)
    MG0 = compute_MG0(MG_obs, Ag) if Ag is not None else None
    bp_rp0 = compute_BP_RP0(bp_rp_obs, 1.3 * EBV) if EBV is not None else None

    # Temperature & Luminosity
    # Note: We pass None for 'spt_for_teff' since we removed manual SIMBAD input
    teff, teff_source = get_teff(None, dr3_result, bp_rp0)
    logTeff = compute_logTeff(teff)
    
    Lum = compute_lum_from_MG0(MG0, bprp=bp_rp0) if MG0 is not None else None
    logL = math.log10(Lum) if (Lum is not None and Lum > 0) else None

    # Spectral Estimation
    spec, hrd_class = estimate_spectral_type_bprp(
        bp_rp0 if bp_rp0 is not None else bp_rp_obs,
        MG0 if MG0 is not None else MG_obs
    )

    # --- 4. ISOCHRONE FITTING (Automated) ---
    print("\n📉 Running Isochrone Fit...")
    
    # OPTIMIZATION: Check if we have required params
    interpolated = None
    R_sun, Age, Teff_iso = None, None, None
    iso_vals = {} # To hold raw iso results

    if logTeff is not None and logL is not None:
        try:
            # Load DF only if we have data to fit
            # Ideally, pass this 'df' into find_best_isochrone if you refactor that too. 
            # For now, we assume standard usage.
            
            # Note: We pass the filename string as requested by the original tool, 
            # but ideally you'd modify isochrone_tool to accept a DF.
            interpolated = find_best_isochrone(
                logTeff, logL, 
                FIT_CONFIG["teff_tol"], 
                FIT_CONFIG["logl_tol"], 
                FIT_CONFIG["min_mass"], 
                str(ISO_PATH)
            )
            
            if interpolated:
                iso_vals = interpolated
                Teff_iso = 10**iso_vals.get("logTe") if iso_vals.get("logTe") else None
                Age = (10**iso_vals.get("logAge"))/1e9 if iso_vals.get("logAge") else None
                R_sun = compute_radius_from_logL_Teff(iso_vals.get("logL"), Teff_iso)
                print("✅ Fit successful.")
            else:
                print("❌ No fit found within default tolerances.")

        except Exception as e:
            print(f"❌ Error during isochrone fitting: {e}")
    else:
        print("⚠️  Insufficient data (Teff/Lum) to run isochrone fit.")

    # --- 5. DATA PACKAGING & SAVING ---
    filename_id = gaia_id  # Use Gaia ID as filename

    data_groups = {}
    
    # Placeholder for SIMBAD (Empty since we removed input)
    data_groups["SIMBAD"] = {k: ("--", "") for k in ["MAIN_ID", "OTYPE", "RA", "DEC", "SP_TYPE", "RV"]}

    data_groups["DR3"] = {
        "Source_ID": (f'"{gaia_id}"', ""),
        "RA": (safe_format(ra, 5), "deg"),
        "DEC": (safe_format(dec, 5), "deg"),
        "Parallax": (safe_format(plx, 3), f"mas ± {safe_format(plx_err, 3)}"),
        "G mag": (safe_format(G_mag, 3), ""),
        "BP-RP": (safe_format(bp_rp_obs, 3), ""),
        "RUWE": (safe_format(dr3_result["ruwe"][0] if "ruwe" in dr3_result.colnames else None, 3), ""),
        "Teff_DR3": (safe_format(teff, 0), "K")
    }

    data_groups["Caltech"] = {
        "E(B–V)_inf": (ebv_inf, "")
    }

    data_groups["Számított"] = {
        "PM_total": (safe_format(pm_total, 2), "mas/yr"),
        "Távolság": (safe_format(distance, 1), f"pc ± {safe_format(distance_err, 1)}"),
        "M_G_obs": (safe_format(MG_obs, 2), ""),
        "M_G0": (safe_format(MG0, 2), ""),
        "BP-RP0": (safe_format(bp_rp0, 3), ""),
        "E(B-V)": (safe_format(EBV, 3), ""),
        "A_G": (safe_format(Ag, 3), ""),
        "logL": (safe_format(logL, 3), ""),
        "logTeff": (safe_format(logTeff, 4), ""),
        "Est. SpT": (spec, ""),
        "HRD Class": (hrd_class, ""),
        "X": (safe_format(x, 1), "pc"),
        "Y": (safe_format(y, 1), "pc"),
        "Z": (safe_format(z, 1), "pc"),
    }

    data_groups["Izokrón"] = {
        "Illesztési Δ": (safe_format(iso_vals.get("fit_delta"), 4), ""),
        "Tömeg": (safe_format(iso_vals.get("mass"), 3), "M\u2299"),
        "Age": (safe_format(Age, 3), "Gyr"),
        "R_sun": (safe_format(R_sun, 2), "R\u2299"),
        "logg": (safe_format(iso_vals.get("logg"), 3), ""),
        "[Fe/H]": (safe_format(iso_vals.get("feh"), 3), "")
    }

    print(f"\n💾 Saving results to {filename_id}_output.txt ...")
    write_structured_output(filename_id, data_groups)
    print("✅ Done.")

if __name__ == "__main__":
    main()
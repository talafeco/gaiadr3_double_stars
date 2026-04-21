import argparse
import warnings
import astropy.units as u
from astropy.coordinates import SkyCoord, Distance
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
import numpy as np

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore', category=UserWarning)

def estimate_star_color(bp_mag, rp_mag):
    """
    Calculates the Gaia BP-RP color index and provides a rough
    estimate of the visual color and spectral class.
    """
    # Check for missing or invalid magnitude data
    if bp_mag is None or rp_mag is None or np.isnan(bp_mag) or np.isnan(rp_mag):
        return None

    color_index = bp_mag - rp_mag

    # Rough approximation of visual color based on Gaia BP-RP
    if color_index < 0.0:
        color_text = "Blue (O/B-type)"
    elif 0.0 <= color_index < 0.5:
        color_text = "Blue-White to White (A-type)"
    elif 0.5 <= color_index < 0.9:
        color_text = "Yellow-White to Yellow (F/G-type)" # The Sun is ~0.82
    elif 0.9 <= color_index < 1.5:
        color_text = "Orange (K-type)"
    else:
        color_text = "Red (M-type or later)"

    return {
        'bp_rp_index': color_index,
        'visual_estimate': color_text
    }

def calculate_stellar_parameters(parallax_mas, apparent_mag):
    """
    Calculates distance (pc and ly), absolute magnitude, and luminosity
    based on parallax and apparent magnitude.
    """
    # Check for invalid, masked, or negative/zero parallax
    if parallax_mas is None or np.isnan(parallax_mas) or parallax_mas <= 0:
        return None

    # 1. Distance Calculation (using Astropy's Distance object)
    distance = Distance(parallax=parallax_mas * u.mas)
    dist_pc = distance.pc
    dist_ly = distance.lyr

    # 2. Absolute Magnitude Calculation
    # M = m - distance_modulus
    dist_mod = distance.distmod.value
    absolute_mag = apparent_mag - dist_mod

    # 3. Luminosity Calculation (L / L_sun)
    # Using the standard Gaia G-band absolute magnitude of the Sun
    M_sun_G = 4.66
    luminosity_Lsun = 10 ** (-0.4 * (absolute_mag - M_sun_G))

    return {
        'distance_pc': dist_pc,
        'distance_ly': dist_ly,
        'absolute_magnitude': absolute_mag,
        'luminosity_Lsun': luminosity_Lsun
    }

def fetch_star_data(gaia_dr3_source_id):
    """
    Downloads star data from Gaia DR3, SIMBAD, 2MASS, and AllWISE
    based on a Gaia DR3 source_id.
    """
    print(f"--- Initiating search for Gaia DR3 {gaia_dr3_source_id} ---")
    data_collection = {}

    # 1. Gaia DR3 Query
    print("Querying Gaia DR3 archive...")
    adql_query = f"""
        SELECT * FROM gaiadr3.gaia_source
        WHERE source_id = {gaia_dr3_source_id}
    """
    job = Gaia.launch_job(adql_query)
    gaia_table = job.get_results()

    if len(gaia_table) == 0:
        print("Error: No source found in Gaia DR3 with that ID.")
        return None

    data_collection['Gaia'] = gaia_table
    print(f"✓ Found Gaia data (G mag: {gaia_table['phot_g_mean_mag'][0]:.2f})")

    ra = gaia_table['ra'][0]
    dec = gaia_table['dec'][0]
    coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')

    # 2. SIMBAD Query
    print("Querying SIMBAD database...")
    custom_simbad = Simbad()
    custom_simbad.add_votable_fields('flux(V)', 'flux(B)', 'sp', 'ids')

    simbad_id = f"Gaia DR3 {gaia_dr3_source_id}"
    simbad_table = custom_simbad.query_object(simbad_id)

    if simbad_table is None:
        print("  Exact ID not found in SIMBAD, attempting coordinate search...")
        simbad_table = custom_simbad.query_region(coord, radius=2*u.arcsec)

    if simbad_table is not None:
        data_collection['SIMBAD'] = simbad_table
        print(f"✓ Found SIMBAD object: {simbad_table['MAIN_ID'][0]}")
    else:
        print("x No SIMBAD data found.")
        data_collection['SIMBAD'] = None

    # 3. VizieR Query (Expanded)
    print("Querying VizieR for additional catalog data...")
    vizier = Vizier(columns=['*'])

    # We can expand the row limit if we expect crowded fields, but for a single star 50 is fine
    vizier.ROW_LIMIT = 50

    # Just add the catalog IDs strings to this list
    catalogs_to_query = [
        "II/246/out",      # 2MASS
        "II/328/allwise",  # AllWISE
        "II/349/ps1",      # Pan-STARRS1
        "II/312/ais",      # GALEX
        "IV/39/tic82"      # TESS Input Catalog
    ]

    vizier_tables = vizier.query_region(
        coord,
        radius=2*u.arcsec,
        catalog=catalogs_to_query
    )

    if len(vizier_tables) > 0:
        if 'II/246/out' in vizier_tables.keys():
            data_collection['2MASS'] = vizier_tables['II/246/out']
            print("✓ Found 2MASS data")
        if 'II/328/allwise' in vizier_tables.keys():
            data_collection['AllWISE'] = vizier_tables['II/328/allwise']
            print("✓ Found AllWISE data")
    else:
         print("x No VizieR catalog matches found within 2 arcseconds.")

    print("--- Search Complete ---\n")
    return data_collection

# ==========================================
# Command Line Interface Setup & Data Output
# ==========================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch cross-catalog data for a Gaia DR3 star.")
    parser.add_argument("target_id", type=int, help="The 19-digit Gaia DR3 source_id of the target star.")
    args = parser.parse_args()

    results = fetch_star_data(args.target_id)

    if results:
        print("=" * 60)
        print("FULL CATALOG DATA EXPORT")
        print("=" * 60)

        # Iterate over every catalog we successfully fetched
        for catalog_name, table_data in results.items():
            if table_data is not None:
                print(f"\n>>> [{catalog_name} Catalog Data] <<<")
                print("-" * 60)

                # Iterate through rows (usually just 1 for exact matches)
                for row_idx, row in enumerate(table_data):
                    print(f"--- Record {row_idx + 1} ---")
                    # Iterate through all columns in the row
                    for col_name in table_data.colnames:
                        # Extract value and handle masked/missing data gracefully
                        val = row[col_name]
                        print(f"{col_name}: {val}")
            else:
                print(f"\n>>> [{catalog_name} Catalog Data] <<<")
                print("No data available.")

        print("\n" + "=" * 60)

# Calculate and print derived stellar physics
        if results.get('Gaia') is not None:
            # Extract raw values
            parallax = results['Gaia']['parallax'][0]
            g_mag = results['Gaia']['phot_g_mean_mag'][0]

            # Pass to our new function
            physics = calculate_stellar_parameters(parallax, g_mag)

            print("\n>>> [Derived Stellar Physics] <<<")
            print("-" * 60)
            if physics:
                print(f"Distance:           {physics['distance_pc']:.2f} pc")
                print(f"Distance:           {physics['distance_ly']:.2f} ly")
                print(f"Absolute Mag (G):   {physics['absolute_magnitude']:.2f}")
                print(f"Luminosity:         {physics['luminosity_Lsun']:.4f} L_sun")
            else:
                print("Cannot calculate physics: Invalid, negative, or missing parallax data.")

        print("\n" + "=" * 60)

# Calculate and print derived stellar physics
        if results.get('Gaia') is not None:
            # Extract raw values for physics
            parallax = results['Gaia']['parallax'][0]
            g_mag = results['Gaia']['phot_g_mean_mag'][0]

            # Extract raw values for color
            bp_mag = results['Gaia']['phot_bp_mean_mag'][0]
            rp_mag = results['Gaia']['phot_rp_mean_mag'][0]

            # Pass to our calculation functions
            physics = calculate_stellar_parameters(parallax, g_mag)
            color_data = estimate_star_color(bp_mag, rp_mag)

            print("\n>>> [Derived Stellar Physics] <<<")
            print("-" * 60)
            if physics:
                print(f"Distance:           {physics['distance_pc']:.2f} pc")
                print(f"Distance:           {physics['distance_ly']:.2f} ly")
                print(f"Absolute Mag (G):   {physics['absolute_magnitude']:.2f}")
                print(f"Luminosity:         {physics['luminosity_Lsun']:.4f} L_sun")
            else:
                print("Cannot calculate physics: Invalid, negative, or missing parallax data.")

            print("-" * 60)
            if color_data:
                print(f"BP-RP Color Index:  {color_data['bp_rp_index']:.3f}")
                print(f"Estimated Color:    {color_data['visual_estimate']}")
            else:
                print("Cannot calculate color: Missing BP or RP magnitude data.")

        print("\n" + "=" * 60)

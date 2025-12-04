from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
from astroquery.irsa_dust import IrsaDust
from astropy.coordinates import SkyCoord
import astropy.units as u

# SIMBAD testreszabása: B, V magni hibák, spektrum, típus, RV + hiba, pozíció szövegként
custom_simbad = Simbad()
custom_simbad.TIMEOUT = 30
custom_simbad.add_votable_fields(
    "main_id", "otype", "ra(d)", "dec(d)",
    "ra", "dec", "sp", "rv_value", "rvz_error",
    "flux(B)", "flux_error(B)", "flux(V)", "flux_error(V)"
)

# --- GAIA ID ellenőrzése ---
def check_gaia_id(gaia_id):
    query = f"SELECT TOP 1 source_id FROM gaiadr3.gaia_source WHERE source_id = {gaia_id}"
    try:
        job = Gaia.launch_job(query)
        result = job.get_results()
        return len(result) > 0
    except Exception as e:
        print(f"Hiba a Gaia lekérdezés során: {e}")
        return False

# --- SIMBAD ID ellenőrzése ---
def check_simbad_id(simbad_id):
    try:
        result = Simbad.query_object(simbad_id)
        return result is not None
    except Exception as e:
        print(f"Hiba a SIMBAD lekérdezés során: {e}")
        return False

def query_simbad_data(simbad_id):
    try:
        result = custom_simbad.query_object(simbad_id)
        if result is not None:
            return result
        else:
            return None
    except Exception as e:
        print(f"Hiba a SIMBAD lekérdezésben: {e}")
        return None

# --- Gaia DR3 adatlekérdezés ---
def query_gaia_dr3_data(gaia_id):
    query = f"""
        SELECT 
            gs.ra, gs.dec, gs.parallax, gs.parallax_error,
            gs.pmra, gs.pmra_error, gs.pmdec, gs.pmdec_error,
            gs.radial_velocity, gs.radial_velocity_error,
            gs.phot_g_mean_mag, gs.bp_rp, gs.ruwe,
            ap.teff_gspphot, ap.logg_gspphot,
            ap.ag_gspphot, ap.mh_gspphot
        FROM gaiadr3.gaia_source AS gs
        LEFT JOIN gaiadr3.astrophysical_parameters AS ap
        ON gs.source_id = ap.source_id
        WHERE gs.source_id = {gaia_id}
    """
    try:
        # job = Gaia.launch_job(query)
        job = Gaia.launch_job_async(query)
        result = job.get_results()
        return result
    except Exception as e:
        print(f"Hiba a DR3 lekérdezésben: {e}")
        return None
        
def get_ebv_inf_from_dr3(dr3_result):
    try:
        if dr3_result is None or 'ra' not in dr3_result.colnames or 'dec' not in dr3_result.colnames:
            return None
        ra_value = dr3_result['ra'][0]; dec_value = dr3_result['dec'][0]
        if ra_value is None or dec_value is None:
            return None
        coord = SkyCoord(ra_value, dec_value, unit="deg")
        table = IrsaDust.get_query_table(coord, section='ebv')
        ebv = float(table['ext SandF mean'][0])
        return ebv
    except Exception as e:
        print(f"Hiba az E(B–V)_inf lekérdezésnél: {e}")
        return None


# def get_ebv_inf_from_dr3(dr3_result):
    # """
    # Lekérdezi az E(B–V)_inf értéket a Caltech vörösödési adatbázisból a DR3 RA/DEC alapján.
    # """
    # try:
        # if dr3_result is None:
            # return '--'
        # if 'ra' not in dr3_result.colnames or 'dec' not in dr3_result.colnames:
            # return '--'

        # ra_value = dr3_result['ra'][0]
        # dec_value = dr3_result['dec'][0]

        # if ra_value is None or dec_value is None:
            # return '--'

        # coord = SkyCoord(ra_value, dec_value, unit="deg")
        # table = IrsaDust.get_query_table(coord, section='ebv')
        # ebv = table['ext SandF mean'][0]
        # return ebv
    # except Exception as e:
        # print(f"Hiba az E(B–V)_inf lekérdezésnél: {e}")
        # return '--'
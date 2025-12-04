import math
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import pandas as pd
import re

mamajek_teff_dict = {
    "O3V": 47860,
    "O4V": 44500,
    "O5V": 41950,
    "O6V": 40000,
    "O7V": 37450,
    "O8V": 35800,
    "O9V": 33700,
    "B0V": 30000,
    "B1V": 25400,
    "B2V": 21400,
    "B3V": 18700,
    "B5V": 15700,
    "B8V": 11900,
    "A0V": 9600,
    "A1V": 9200,
    "A2V": 8840,
    "A3V": 8600,
    "A5V": 8080,
    "F0V": 7420,
    "F2V": 7050,
    "F5V": 6530,
    "F6V": 6400,
    "F8V": 6200,
    "G0V": 5940,
    "G2V": 5770,
    "G5V": 5660,
    "G8V": 5470,
    "K0V": 5250,
    "K2V": 4960,
    "K5V": 4400,
    "K7V": 4060,
    "M0V": 3850,
    "M1V": 3700,
    "M2V": 3550,
    "M3V": 3410,
    "M4V": 3240,
    "M5V": 3090,
    "M6V": 2960,
    "M7V": 2850,
    "M8V": 2710
}

import math

def compute_lum_from_MG0(MG0, bprp=None):
    """
    Calculates Luminosity (in Solar units) from Absolute G Magnitude.
    
    Parameters:
    -----------
    MG0 : float
        Absolute Magnitude in G-band (extinction corrected).
    bprp : float, optional
        The BP-RP color index. Used to estimate Bolometric Correction (BC).
        If None, BC is assumed to be 0 (Sun-like approximation).
        
    Returns:
    --------
    float : Luminosity (L_sun)
    """
    if MG0 is None:
        return None

    try:
        MG0 = float(MG0)
        
        # 1. Solar Bolometric Magnitude (IAU 2015 resolution)
        M_bol_sun = 4.74
        
        # 2. Estimate Bolometric Correction (BC_G)
        # Empirical polynomial approximation for Gaia G-band based on BP-RP color.
        # This is a simplified fit suitable for Main Sequence stars.
        BC_G = 0.0
        
        if bprp is not None:
            try:
                c = float(bprp)
                # Simple polynomial fit region for Main Sequence
                # Very Cool stars (Red Dwarfs, BP-RP > 1.5): Significant negative correction
                if c > 1.5:
                    BC_G = -0.05 * (c - 1.5)**2 - 0.1
                # Hot stars (BP-RP < 0): Correction for UV output
                elif c < 0:
                     BC_G = -1.5 * c**2
                # Sun-like / Intermediate: BC is small, close to 0 or slightly negative
                else:
                    BC_G = -0.05 * c  # Gentle slope
            except ValueError:
                BC_G = 0.0
        
        # 3. Calculate Bolometric Magnitude of the star
        M_bol_star = MG0 + BC_G

        # 4. Calculate Luminosity relative to Sun
        # Formula: M_bol_star - M_bol_sun = -2.5 * log10(L_star / L_sun)
        # Rearranged:
        luminosity = 10 ** (0.4 * (M_bol_sun - M_bol_star))
        
        return luminosity

    except (ValueError, TypeError):
        return None

def get_clean_spt_from_simbad(simbad_data, mamajek_dict=None):
    """
    A SIMBAD SP_TYPE mező megtisztítása és egyszerűsítése Mamajek-kulcshoz.
    Példák:
        'G2V+WD'     → 'G2V'
        'B9.5Ve'     → 'B95V'
        'A1Vn'       → 'A1V'
        'M1IIIe'     → 'M1III'
    """
    if simbad_data is None or "SP_TYPE" not in simbad_data.colnames:
        return None

    raw_spt = str(simbad_data["SP_TYPE"][0]).strip().upper()

    # Csak az első komponens (ha többes rendszer)
    spt = raw_spt.split("+")[0]

    # Tisztítás
    spt = (
        spt.replace("VE", "V")   # emission típusok normalizálása
           .replace("VZ", "V")   # nagyon fiatal szubfőág
           .replace(".", "")    # pont eltávolítása (pl. B9.5 → B95)
           .replace(" ", "")
    )
    
    if mamajek_dict:
        for n in range(4, 1, -1):  # Próbálkozás 4→3→2 karakterrel
            prefix = spt[:n]
            matches = [k for k in mamajek_dict if k.startswith(prefix)]
            if matches:
                return matches[0]  # legelső találat
        return None  # nincs megfelelő
    else:
        return spt
    
# def load_mamajek_table(filepath="D:/Astro/catalogs/mamajek.txt"):
    # """
    # Mamajek-tábla beolvasása, visszaadja a Teff és BCv értékeket szótárként.
    # """
    # try:
        # df = pd.read_csv(filepath, sep=r"\s+", comment="#", skip_blank_lines=True)
        
        # Oszlopok kiválasztása és normalizálás
        # df = df.rename(columns={"SpT": "SpT", "Teff": "Teff", "BCv": "BCv"})
        # df = df[["SpT", "Teff", "BCv"]].dropna()
        # df["SpT"] = df["SpT"].str.upper()
        
        # Két szótárba rendezés
        # teff_dict = df.set_index("SpT")["Teff"].to_dict()
        # bc_dict = df.set_index("SpT")["BCv"].to_dict()

        # return teff_dict, bc_dict
    # except Exception as e:
        # print(f"⚠️ Mamajek-tábla beolvasási hiba: {e}")
        # return {}, {}

def get_teff(SpT, dr3_result, bp_rp0):
    """
    Effektív hőmérséklet meghatározása prioritási sorrendben:
    1️⃣ Mamajek SpT alapján (ha van)
    2️⃣ Gaia DR3 teff_gspphot (ha van)
    3️⃣ BP–RP0 színindex alapján számított érték
    """
    if isinstance(SpT, str) and SpT.strip():
        SpT_clean = SpT.strip().upper().replace("VE", "V").replace("VZ", "V").replace(".", "").replace(" ", "")
        for n in range(4, 1, -1):
            prefix = SpT_clean[:n]
            matches = [k for k in mamajek_teff_dict.keys() if k.startswith(prefix)]
            if matches:
                match = matches[0]
                return mamajek_teff_dict[match], f"mamajek ({match})"

    # 2️⃣ Gaia DR3 Teff
    if dr3_result is not None and "teff_gspphot" in dr3_result.colnames:
        teff_raw = dr3_result["teff_gspphot"][0]
        if teff_raw is not np.ma.masked and not np.isnan(teff_raw):
            return float(teff_raw), "dr3_teff"
            
    # 3️⃣ BP–RP0 színindex alapján becslés
    if bp_rp0 not in [None, "--"]:
        teff = estimate_teff_from_bp_rp(bp_rp0)
        if teff:
            return teff, "bp_rp0"

    return None, "none"
    
def estimate_teff_from_bp_rp(bp_rp0):
    """
    BP–RP színindex alapján becsült effektív hőmérséklet (K).
    Jordi et al. 2010 – Gaiához illeszkedő polinomiális becslés
    Csak 0.5 < BP–RP < 2.0 tartományban megbízható.
    """
    if bp_rp0 is None or bp_rp0 == "--":
        return None
    try:
        x = float(bp_rp0)
        if x < -0.5 or x > 4.0:
            return None  # abszurd érték
        if x < 0.5 or x > 2.0:
            print(f"⚠️ Figyelem: a BP–RP₀ = {x:.3f} kívül esik a megbízható 0.5–2.0 tartományon, becslés csak tájékoztató jellegű.")
        
        # Példa: polinomiális illesztés (csak illusztráció!)
        teff = 9110 - 6820 * x + 2630 * x**2 - 400 * x**3
        return round(teff, 1)
    except:
        return None

def compute_logTeff(teff):
    if teff and teff > 0:
        return math.log10(teff)
    return None

def compute_distance_and_error(parallax, parallax_error):
    if parallax and parallax > 0:
        dist = 1000 / parallax
        dist_err = 1000 * parallax_error / (parallax ** 2) if parallax_error else None
        return dist, dist_err
    return None, None

def compute_MG_obs(G_mag, distance):
    if G_mag is not None and distance:
        return G_mag - 5 * math.log10(distance) + 5
    return None

def compute_EBV(ebv_inf, distance, gal_lat_deg):
    if ebv_inf is None or distance is None or gal_lat_deg is None:
        return None
    sinb = math.sin(math.radians(abs(gal_lat_deg)))
    factor = 1 - math.exp(-0.008 * distance * sinb)
    return ebv_inf * factor

    
# def compute_Ag(MG_obs, ebv):
    # """
    # Gaia G sávra becsült extinkció E(B–V) alapján.
    # Ez egy analóg képlet az Av-hez, de MG-re alkalmazva.
    # """
    # if ebv is not None and MG_obs is not None:
        # return (2.72 + 0.15 * MG_obs + 0.02 * ebv) * ebv
    # return None
   
def compute_MG0(MG_obs, Ag):
    """
    Korrigált abszolút G magnitúdó számítása az AG extinkció alapján.
    """
    if MG_obs is not None and Ag is not None:
        return MG_obs - Ag
    return None

def compute_BP_RP0(bp_rp_obs, ebp_rp):
    if bp_rp_obs is not None and ebp_rp is not None:
        return bp_rp_obs - ebp_rp
    return None
    


def estimate_spectral_type_bprp(bp_rp0=None, MG0=None, bp_rp_obs=None, MG_obs=None, spec_hint=None):
    """
    Spektráltípus + LC becslése dereddenelt adatokból (elsőbbségben BP–RP0, MG0).
    Ha ezek nem elérhetők, fallback az obs értékekre.
    Vissza: pl. "G5 IV*", "K0 V*", vagy "ismeretlen".
    """

    # --- Dereddenelt prioritás ---
    bp_rp = bp_rp0 if bp_rp0 is not None else bp_rp_obs
    MG    = MG0    if MG0    is not None else MG_obs

    if bp_rp is None:
        return "ismeretlen", "ismeretlen"

    try:
        bp_rp = float(bp_rp)
        MG    = float(MG)
    except Exception:
        return "ismeretlen", "ismeretlen"

    # Konzervatív tartomány – ha nagyon kék/piros, nincs becslés
    # if bp_rp < 0.5 or bp_rp > 2.0:
        # return "ismeretlen", "ismeretlen"

    # --- SpT becslés a BP–RP alapján ---
    bp_rp_sp = [
        (-0.50, -0.44, "O3"),
        (-0.44, -0.41, "O4"),
        (-0.41, -0.38, "O5"),
        (-0.38, -0.35, "O6"),
        (-0.35, -0.32, "O7"),
        (-0.32, -0.28, "O8"),
        (-0.28, -0.20, "O9"),
        (-0.20, -0.10, "B0"),
        (-0.10, -0.02, "B1"),
        (-0.02,  0.05, "B2"),
        ( 0.05,  0.10, "B3"),
        ( 0.10,  0.15, "B4"),
        ( 0.15,  0.23, "B5"),
        ( 0.23,  0.26, "B7"),
        ( 0.26,  0.29, "B8"),
        ( 0.29,  0.32, "B9"),
        ( 0.32,  0.38, "A0"),
        ( 0.38,  0.44, "A1"),
        ( 0.44,  0.50, "A2"),
        ( 0.50,  0.56, "A5"),
        ( 0.56,  0.62, "F0"),
        ( 0.62,  0.70, "F2"),
        ( 0.70,  0.78, "F5"),
        ( 0.78,  0.86, "F8"),
        ( 0.86,  0.94, "G0"),
        ( 0.94,  1.00, "G2"),
        ( 1.00,  1.06, "G5"),
        ( 1.06,  1.12, "G8"),
        ( 1.12,  1.18, "K0"),
        ( 1.18,  1.24, "K2"),
        ( 1.24,  1.30, "K3"),
        ( 1.30,  1.36, "K5"),
        ( 1.36,  1.42, "K7"),
        ( 1.42,  1.48, "M0"),
        ( 1.48,  1.54, "M1"),
        ( 1.54,  1.60, "M2"),
        ( 1.60,  1.66, "M3"),
        ( 1.66,  1.72, "M4"),
        ( 1.72,  1.78, "M5"),
        ( 1.78,  1.84, "M6"),
        ( 1.84,  1.90, "M7"),
        ( 1.90,  2.00, "M8"),
    ]

    spt = "?"
    for lo, hi, sp in bp_rp_sp:
        if lo <= bp_rp < hi:
            spt = sp
            break
    if spt == "?":
        return "ismeretlen", "ismeretlen"

    # --- LC becslés ---
    lc = None

    # 1) Regex a spec_hintből
    if spec_hint:
        m = re.search(r"\b(I|II|III|IV|V|WD)\b", str(spec_hint))
        if m:
            lc = m.group(1)

    # 2) Ha nincs LC a hintben, MG0 és BP–RP alapján
    if lc is None:
        if bp_rp < 0.8:  # forróbb (A–F típus)
            if MG < 0.0:
                lc = "III"
            elif MG < 3.5:
                lc = "V"
            else:
                lc = "WD"
        elif bp_rp < 1.5:  # közepes (G–K típus)
            if MG < 1.0:
                lc = "III"
            elif MG < 4.5:
                lc = "V"
            else:
                lc = "WD"
        else:  # vörös (K–M típus)
            if MG < -0.5:
                lc = "II"
            elif MG < 2.0:
                lc = "III"
            elif MG < 7.5:
                lc = "V"
            else:
                lc = "WD"

    hrd_class = classify_hrd_position(lc) if lc else "ismeretlen, ismeretlen"
    return f"{spt} {lc}*", hrd_class
    
def estimate_BCg_from_bprp(bp_rp):
    """
    Bolometrikus korrekció (BCg) becslése Gaia BP–RP színindex alapján.
    Durva közelítés!
    """
    if bp_rp is None:
        return None
    try:
        bp_rp = float(bp_rp)
        if bp_rp < -0.2:
            return -1.5
        elif bp_rp < 0.0:
            return -1.0
        elif bp_rp < 0.5:
            return -0.4
        elif bp_rp < 1.0:
            return -0.1
        elif bp_rp < 1.5:
            return 0.2
        elif bp_rp < 2.0:
            return 0.5
        else:
            return 0.8
    except:
        return None

def compute_lum_from_MG0(MG0, BCg=None, bprp=None):
    """
    Gaia G-sáv bolometrikus korrekció alapján számítja a luminozitást.
    Ha a BCg nincs megadva, becsli BP–RP alapján.
    """
    if BCg is None:
        BCg = estimate_BCg_from_bprp(bprp)

    if MG0 is not None and BCg is not None:
        try:
            return 10 ** ((4.74 - (MG0 + BCg)) / 2.5)
        except Exception as e:
            print(f"⚠️ Hiba a Lum számításnál (G): {e}")
    return None

def compute_radius_from_logL_Teff(logL, teff):
    """
    Kiszámítja a csillag sugarát a logL és a Teff értékek alapján.
    R / R\u2299 = sqrt(L / L\u2299) / (Teff / 5772)^2
    vagy logR = 0.5 * logL - 2 * log(Teff / 5772)
    """
    if logL is not None and teff is not None and teff > 0:
        try:
            logR = 0.5 * logL - 2 * math.log10(teff / 5772)
            return round(10**logR, 3)
        except:
            return None
    return None

def classify_hrd_position(spectral_classification):
    """
    HRD-besorolás szöveges spektráltípus alapján (pl. 'G5 IV*' → 'alóóriás').
    """
    if not spectral_classification:
        return "ismeretlen"
    if isinstance(spectral_classification, tuple):
        spectral_classification = spectral_classification[0]

    spectral_classification = spectral_classification.upper()

    # Vegyes (IV-V) eset
    if "IV-V" in spectral_classification or "V-IV" in spectral_classification:
        return "fősorozat"  # vegyes típus, főághoz közeli

    # FONTOS: az IV ellenőrzése jöjjön a V ELÉ!
    if "IV" in spectral_classification:
        return "alóriás"
    elif "III" in spectral_classification:
        return "óriás"
    elif "II" in spectral_classification or "I" in spectral_classification:
        return "szuperóriás"
    elif "V" in spectral_classification:
        return "fősorozat"

    return "ismeretlen"

        
def compute_galactic_xyz(ra_deg, dec_deg, distance_pc):
    """
    RA, DEC és távolság alapján kiszámítja a galaktikus koordinátákat (l, b)
    és az x, y, z térbeli pozíciót parsecben.
    """
    coord = SkyCoord(ra=ra_deg*u.degree, dec=dec_deg*u.degree, distance=distance_pc*u.pc, frame='icrs')
    gal = coord.galactic

    l = gal.l.deg
    b = gal.b.deg

    l_rad = np.radians(l)
    b_rad = np.radians(b)

    x = distance_pc * np.cos(b_rad) * np.cos(l_rad)
    y = distance_pc * np.cos(b_rad) * np.sin(l_rad)
    z = distance_pc * np.sin(b_rad)

    return round(l, 4), round(b, 4), round(x, 3), round(y, 3), round(z, 3)

def compute_pm_total(pmra, pmdec):
    """
    Kiszámítja a teljes sajátmozgást (PM) milliarcmásodperc/év (mas/yr) egységben,
    a pmRA és pmDEC komponensekből.
    """
    if pmra is not None and pmdec is not None:
        try:
            return math.sqrt(pmra**2 + pmdec**2)
        except:
            return None
    return None

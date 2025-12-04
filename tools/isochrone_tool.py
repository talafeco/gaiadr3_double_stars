# isochrone_tool.py
import pandas as pd
import numpy as np

def load_isochrones(filepath):
    """
    Betölti a formázott izokron fájlt szóköz alapon.
    """
    df = pd.read_csv(filepath, sep=r"\s+")
    print(f"\n📊 Izokron fájl statisztikák:")
    print(f"🔸 Tömegtartomány:     {df['star_mass'].min():.2f} – {df['star_mass'].max():.2f} M\u2299")
    print(f"🔸 log_Teff tartomány: {df['log_Teff'].min():.3f} – {df['log_Teff'].max():.3f}")
    print(f"🔸 log_L tartomány:    {df['log_L'].min():.3f} – {df['log_L'].max():.3f}")
    return df
    
def get_filter_parameters(default_teff_tol=0.1, default_logl_tol=0.4, default_min_mass=1.0):
    """
    Szűrési paraméterek interaktív bekérése a dinamikus szűréshez.
    Ha nem módosítasz, az alapértékek maradnak.
    """
    print("\n🎛 Szűrési beállítások:")
    print(f"   1. Teff tolerancia: ±{default_teff_tol}")
    print(f"   2. Luminozitás tolerancia: ±{default_logl_tol}")
    print(f"   3. Minimális csillagtömeg: > {default_min_mass} M\u2299")

    modify = input("🔧 Szeretnéd módosítani ezeket? (i/n): ").strip().lower()
    if modify != "i":
        return default_teff_tol, default_logl_tol, default_min_mass

    # Új értékek bekérése
    try:
        teff_tol = float(input("   Új Teff tolerancia (pl. 0.3): ").strip())
    except ValueError:
        print("⚠️ Hibás érték, marad az alapértelmezett.")
        teff_tol = default_teff_tol

    try:
        logl_tol = float(input("   Új Lum. tolerancia (pl. 0.2): ").strip())
    except ValueError:
        print("⚠️ Hibás érték, marad az alapértelmezett.")
        logl_tol = default_logl_tol

    try:
        min_mass = float(input("   Új minimális tömeg (pl. 1.0): ").strip())
    except ValueError:
        print("⚠️ Hibás érték, marad az alapértelmezett.")
        min_mass = default_min_mass

    return teff_tol, logl_tol, min_mass
    
def dynamic_filter(df, logTeff, logL, teff_tol, logl_tol, min_mass, ms_only=True):

    print("🔍 Dinamikus szűrés: előnyben a nagy tömegű, forró csillagok.")
    print(f"   Teff ±{teff_tol}, Lum ±{logl_tol}, Tömeg > {min_mass} M\u2299")
    
    filtered = df[
        (df["star_mass"] > min_mass) &                 # ↓ lazább tömegküszöb
        (df["log_Teff"] > logTeff - teff_tol) &
        (df["log_Teff"] < logTeff + teff_tol) &
        (df["log_L"] > logL - logl_tol) &
        (df["log_L"] < logL + logl_tol)
    ]

    print(f"📉 Szűrés után {len(filtered)} izokrón pont maradt.")
    return filtered
   
def match_isochrone_point(df, logTeff, logL, return_neighbors=False, n_neighbors=5):
    df = df.copy()
    df["distance"] = ((df["log_Teff"] - logTeff)**2 + (df["log_L"] - logL)**2)**0.5
    df_sorted = df.sort_values(by="distance")

    best_row = df_sorted.iloc[0]
    min_distance = best_row["distance"]

    if return_neighbors:
        neighbor_df = df_sorted.iloc[:n_neighbors]
        return best_row, min_distance, neighbor_df
    else:
        return best_row
        
def find_best_isochrone(logTeff, logL, teff_tol=0.08, logl_tol=0.4, min_mass=0.1, filepath=None):
    if filepath is None:
        raise ValueError("Az izokrón fájl útvonala nincs megadva. Add át paraméterként vagy állítsd be az ISO_DIR-t.")
    print(f"\n💫 Csillag paraméterei: logTeff = {logTeff:.4f}, logL = {logL:.4f}")

    df = load_isochrones(filepath)

    df_filtered = dynamic_filter(df, logTeff, logL, teff_tol, logl_tol, min_mass)
    if df_filtered.empty:
        print("❌ Nincs illeszkedő izokrón pont.")
        return None

    best_point, best_distance, neighbors = match_isochrone_point(df_filtered, logTeff, logL, return_neighbors=True)

    # Interpolációs súlyok: 1/d, biztosítva hogy legyen 'distance' az adatokban
    distances = np.sqrt((neighbors["log_Teff"]-logTeff)**2 + (neighbors["log_L"]-logL)**2)
    neighbors = neighbors.assign(distance=distances.values)

    interpolated = interpolate_isochrone_values(best_point, neighbors, distances)
    interpolated["fit_delta"] = float(best_distance)  # <-- EGY név

    return interpolated

def interpolate_isochrone_values(center_point, neighbor_df, distances):
    eps = 1e-6
    w = 1.0 / (neighbor_df["distance"].to_numpy() + eps)
    W = w.sum()
    return {
        "mass":   float((neighbor_df["star_mass"].to_numpy() * w).sum() / W),
        "logAge": float((neighbor_df["logAge"].to_numpy()   * w).sum() / W),
        "logg":   float((neighbor_df["log_g"].to_numpy()    * w).sum() / W),
        "logTe":  float((neighbor_df["log_Teff"].to_numpy() * w).sum() / W),
        "logL":   float((neighbor_df["log_L"].to_numpy()    * w).sum() / W),
        "feh":    float((neighbor_df["[Fe/H]"].to_numpy()   * w).sum() / W),
    }


# def find_best_isochrone(logTeff, logL, teff_tol=0.08, logl_tol=0.4, min_mass=0.1, filepath="d:\\Astro\\catalogs\\izokrona_rendezett.txt"):
    # """
    # Komplett izokron-illesztő folyamat:
    # - fősorozat-ellenőrzés
    # - izokron betöltése
    # - szűrés a fizikai tartományra
    # - legjobb pont kiválasztása
    # - 3D interpoláció a legközelebbi N pontból
    # A visszatérés egy dict az interpolált értékekkel.
    # """

    # print(f"\n💫 Csillag paraméterei: logTeff = {logTeff:.4f}, logL = {logL:.4f}")

    # # 💡 HRD besorolás és fősorozat meghatározása
    # # hrd_class = classify_hrd_position(lc)
    # # main_sequence = hrd_class == "fősorozat"

    # # 📂 Izokrón fájl betöltése
    # df = load_isochrones(filepath)

    # # 🔍 Dinamikus szűrés a megadott toleranciákkal
    # df_filtered = dynamic_filter(
        # df, logTeff, logL,
        # teff_tol, logl_tol, min_mass,
    # )

    # if df_filtered.empty:
        # print("❌ Nincs illeszkedő izokrón pont.")
        # return None

    # # 🔎 Legjobb pont kiválasztása
    # best_point, best_distance, neighbors = match_isochrone_point(df_filtered, logTeff, logL, return_neighbors=True)

    # # 📏 Távolságok újraszámítása az interpolációhoz
    # distances = np.sqrt(
        # (neighbors["log_Teff"] - logTeff) ** 2 +
        # (neighbors["log_L"] - logL) ** 2
    # )

    # # 🔄 Interpoláció
    # interpolated = interpolate_isochrone_values(best_point, neighbors, distances)
    # interpolated["Δ"] = best_distance  # az illesztési hiba

    # return interpolated

# def interpolate_isochrone_values(center_point, neighbor_df, distances):
    # epsilon = 1e-6
    # weights = 1 / (neighbor_df["distance"] + epsilon)
    # total_weight = weights.sum()

    # interpolated = {
        # "mass": (neighbor_df["star_mass"] * weights).sum() / total_weight,
        # "logAge": (neighbor_df["logAge"] * weights).sum() / total_weight,
        # "logg": (neighbor_df["log_g"] * weights).sum() / total_weight,
        # "logTe": (neighbor_df["log_Teff"] * weights).sum() / total_weight,
        # "logL": (neighbor_df["log_L"] * weights).sum() / total_weight,
        # "feh": (neighbor_df["[Fe/H]"] * weights).sum() / total_weight,
        # "delta": np.mean(distances[:len(weights)]),
    # }

    # return interpolated


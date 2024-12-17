def stellar_classification(bp_magnitude, rp_magnitude):
    """
    Classifies a star based on Gaia DR3 blue (G_BP) and red (G_RP) magnitudes.

    Parameters:
        bp_magnitude (float): Gaia DR3 G_BP magnitude.
        rp_magnitude (float): Gaia DR3 G_RP magnitude.

    Returns:
        str: Stellar classification (O, B, A, F, G, K, or M).
    """
    # Calculate the color index
    color_index = bp_magnitude - rp_magnitude

    spectral_type = "Unknown"
    star_color = "Unknown"

    # Stellar classification based on color index
    if color_index < -0.4:
        spectral_type = "O"
        star_color = "Blue"
    elif -0.4 <= color_index < 0.0:
        spectral_type = "B"
        star_color = "Blue-White"
    elif 0.0 <= color_index < 0.4:
        spectral_type = "A"
        star_color = "White"
    elif 0.4 <= color_index < 0.8:
        spectral_type = "F"
        star_color = "Yellow-White"
    elif 0.8 <= color_index < 1.15:
        spectral_type = "G"
        star_color = "Yellow"
    elif 1.15 <= color_index < 1.6:
        spectral_type = "K"
        star_color = "Orange"
    elif color_index >= 1.6:
        spectral_type = "M"
        star_color = "Red"

    return spectral_type, star_color, color_index

def calculate_temperature(bp_magnitude, rp_magnitude):
    """
    Estimates the surface temperature of a star based on Gaia DR3 blue (G_BP) and red (G_RP) magnitudes.

    Parameters:
        bp_magnitude (float): Gaia DR3 G_BP magnitude.
        rp_magnitude (float): Gaia DR3 G_RP magnitude.

    Returns:
        float: Estimated surface temperature in Kelvin (K).
    """
    # Calculate the color index
    color_index = bp_magnitude - rp_magnitude

    # Approximate constants for the color-temperature relation
    # Derived from empirical stellar data
    if color_index < 0.0:
        # Blue stars (hot)
        a, b = -0.25, 3.9
    elif 0.0 <= color_index < 1.5:
        # Yellow/white stars (medium temperature)
        a, b = -0.15, 3.8
    else:
        # Red stars (cool)
        a, b = -0.10, 3.7

    # Calculate log10(T)
    log_temperature = a * color_index + b

    # Convert to temperature in Kelvin
    temperature = 10 ** log_temperature

    return temperature

def classify_star_by_temperature(temperature):
    """
    Classifies a star's spectral type and color based on its surface temperature.

    Parameters:
        temperature (float): Surface temperature of the star in Kelvin (K).

    Returns:
        dict: Dictionary with spectral class and color.
    """
    if temperature > 30000:
        spectral_class = "O"
        color = "Blue"
    elif 10000 <= temperature <= 30000:
        spectral_class = "B"
        color = "Blue-white"
    elif 7500 <= temperature < 10000:
        spectral_class = "A"
        color = "White"
    elif 6000 <= temperature < 7500:
        spectral_class = "F"
        color = "Yellow-white"
    elif 5200 <= temperature < 6000:
        spectral_class = "G"
        color = "Yellow"
    elif 3700 <= temperature < 5200:
        spectral_class = "K"
        color = "Orange"
    elif temperature < 3700:
        spectral_class = "M"
        color = "Red"
    else:
        spectral_class = "Unknown"
        color = "Unknown"

    return {
        "Spectral Class": spectral_class,
        "Color": color
    }


# Example Usage
# Example Gaia DR3 magnitudes
bp_magnitude = 12.557045
g_magnitude = 12.163539
rp_magnitude = 11.595571
surface_temp = 6508.499

#stellar_class, star_color, color_index = stellar_classification(bp_magnitude, rp_magnitude)
stellar_class, star_color, color_index = stellar_classification(bp_magnitude, g_magnitude)
print(f"Stellar Classification: {stellar_class}")
print(f"Star Color: {star_color}")
print(f"Star Index: {color_index}")


#temperature = calculate_temperature(bp_magnitude, rp_magnitude)
temperature = calculate_temperature(bp_magnitude, g_magnitude)
print(f"The star's approximate surface temperature is: {temperature:.2f} K")

# Example Usage
classification = classify_star_by_temperature(surface_temp)
print(f"Surface Temperature: {surface_temp} K")
print(f"Spectral Class: {classification['Spectral Class']}")
print(f"Color: {classification['Color']}")
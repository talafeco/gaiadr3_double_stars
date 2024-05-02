#! /usr/bin/python3

import numpy as np

def gravitational_bond(distance_au, mass1_solar_mass, mass2_solar_mass):
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
    M_sun = 1.989e30  # Mass of the Sun in kg
    
    distance_m = distance_au * 1.496e11  # Convert AU to meters
    distance_pc = distance_m / 3.086e16  # Convert meters to parsecs
    
    mass1_kg = mass1_solar_mass * M_sun  # Convert solar masses to kg
    mass2_kg = mass2_solar_mass * M_sun
    
    gravitational_force = G * mass1_kg * mass2_kg / distance_m**2  # Calculate gravitational force
    
    return gravitational_force

def velocity_difference(proper_motion_ra1_mas, proper_motion_dec1_mas, proper_motion_ra2_mas, proper_motion_dec2_mas, distance_pc, radial_velocity1_km_s, radial_velocity2_km_s):
    mas_to_km = distance_pc * 3.086e13  # Convert milliarcseconds to kilometers
    
    delta_ra_km_s = proper_motion_ra1_mas - proper_motion_ra2_mas
    delta_dec_km_s = proper_motion_dec1_mas - proper_motion_dec2_mas
    
    delta_v_ra_dec = np.sqrt(delta_ra_km_s**2 + delta_dec_km_s**2)  # Calculate velocity difference in km/s
    
    delta_radial_velocity = radial_velocity1_km_s - radial_velocity2_km_s
    
    total_velocity_difference = np.sqrt(delta_v_ra_dec**2 + delta_radial_velocity**2)
    
    return total_velocity_difference

def escape_velocity(distance_au, mass1_solar_mass, mass2_solar_mass):
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
    M_sun = 1.989e30  # Mass of the Sun in kg
    
    distance_m = distance_au * 1.496e11  # Convert AU to meters
    distance_pc = distance_m / 3.086e16  # Convert meters to parsecs
    
    mass1_kg = mass1_solar_mass * M_sun  # Convert solar masses to kg
    mass2_kg = mass2_solar_mass * M_sun
    
    total_mass = mass1_kg + mass2_kg
    
    escape_velocity = np.sqrt(2 * G * total_mass / distance_m)  # Calculate escape velocity
    
    return escape_velocity

def orbiting_status(distance_au, mass1_solar_mass, mass2_solar_mass, proper_motion_ra1_mas, proper_motion_dec1_mas, proper_motion_ra2_mas, proper_motion_dec2_mas, distance_pc, radial_velocity1_km_s, radial_velocity2_km_s):
    grav_force = gravitational_bond(distance_au, mass1_solar_mass, mass2_solar_mass)
    vel_difference = velocity_difference(proper_motion_ra1_mas, proper_motion_dec1_mas, proper_motion_ra2_mas, proper_motion_dec2_mas, distance_pc, radial_velocity1_km_s, radial_velocity2_km_s)
    esc_velocity = escape_velocity(distance_au, mass1_solar_mass, mass2_solar_mass)
    print(esc_velocity, vel_difference)
    
    if vel_difference <= esc_velocity:
        return "The stars are orbiting each other."
    else:
        return "The stars are not orbiting each other."

# Example usage:
distance_au = 3040  # Example distance in AU
mass1_solar_mass = 1.27  # Example mass of star 1 in solar masses
mass2_solar_mass = 0.7  # Example mass of star 2 in solar masses
proper_motion_ra1_mas = -66.71272791555867  # Example proper motion in RA for star 1 in milliarcseconds
proper_motion_dec1_mas = -2.8372787751705295  # Example proper motion in Dec for star 1 in milliarcseconds
proper_motion_ra2_mas = -66.18234160241312  # Example proper motion in RA for star 2 in milliarcseconds
proper_motion_dec2_mas = -3.970262523025885  # Example proper motion in Dec for star 2 in milliarcseconds
distance_pc = 121  # Example distance from Earth in parsecs
radial_velocity1_km_s = -25.264393  # Example radial velocity of star 1 in km/s
radial_velocity2_km_s = -24.513504  # Example radial velocity of star 2 in km/s

status = orbiting_status(distance_au, mass1_solar_mass, mass2_solar_mass, proper_motion_ra1_mas, proper_motion_dec1_mas, proper_motion_ra2_mas, proper_motion_dec2_mas, distance_pc, radial_velocity1_km_s, radial_velocity2_km_s)
print(status)

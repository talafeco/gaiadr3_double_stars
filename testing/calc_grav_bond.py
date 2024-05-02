'''
To calculate the probability of gravitational binding between two stars, we need to compare the escape velocity of the system with the relative velocity of the stars. The escape velocity represents the minimum velocity required for an object to escape the gravitational influence of another object. If the relative velocity of the stars is less than the escape velocity, then they are bound together by gravity.
This function takes as input the masses of the two stars (in solar masses), the distance between them (in AU), and their 3D motion vectors (in km/s). It returns the probability of gravitational binding between the stars. Remember to provide the velocity vectors of the stars in the same reference frame.
'''

import numpy as np

def calculate_gravitational_binding_probability(mass_star1, mass_star2, distance, velocity_vector_star1, velocity_vector_star2):
    # Constants
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
    solar_mass_to_kg = 1.9885e30  # Conversion factor from solar mass to kg
    au_to_m = 1.496e11  # Conversion factor from AU to meters

    # Convert masses to kg
    mass1_kg = mass_star1 * solar_mass_to_kg
    mass2_kg = mass_star2 * solar_mass_to_kg

    # Convert distance to meters
    distance_m = distance * au_to_m

    # Calculate the gravitational force between the stars
    gravitational_force = G * mass1_kg * mass2_kg / distance_m**2

    # Calculate the escape velocity of the system
    escape_velocity = np.sqrt(2 * gravitational_force / mass1_kg)

    # Calculate relative velocity between the stars
    relative_velocity = np.linalg.norm(velocity_vector_star1 - velocity_vector_star2)

    # Calculate the probability of gravitational binding
    probability_binding = relative_velocity / escape_velocity

    return probability_binding

# Example usage
mass_star1 = 1.5  # Mass of star 1 in solar masses
mass_star2 = 2.0  # Mass of star 2 in solar masses
distance = 10.0  # Distance between stars in AU
velocity_vector_star1 = np.array([10, 0, 0])  # 3D velocity vector of star 1 in km/s
velocity_vector_star2 = np.array([-5, 0, 0])  # 3D velocity vector of star 2 in km/s

probability = calculate_gravitational_binding_probability(mass_star1, mass_star2, distance, velocity_vector_star1, velocity_vector_star2)
print("Probability of gravitational binding:", probability)
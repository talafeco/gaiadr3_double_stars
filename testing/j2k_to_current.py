from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.table import QTable

def transform_coordinates_to_current_date(qtable, current_date):
    # Convert the current date to an Astropy Time object
    time = Time(current_date)

    # Create an empty list to store transformed coordinates
    transformed_coords = []

    # Iterate over each row in the QTable
    for row in qtable:
        # Extract RA and Dec from the row
        ra = row['ra']
        dec = row['dec']

        # Create a SkyCoord object for the J2000 coordinates
        coord_j2000 = SkyCoord(ra=ra, dec=dec, unit='deg', equinox='J2000')

        # Transform the coordinates to the current date
        coord_current = coord_j2000.transform_to(time)

        # Append the transformed coordinates to the list
        transformed_coords.append(coord_current)

    # Add a new column to the QTable with the transformed coordinates
    qtable['ra_current'] = [coord.ra.deg for coord in transformed_coords]
    qtable['dec_current'] = [coord.dec.deg for coord in transformed_coords]

    return qtable

# Example usage:
# Assuming qtable is your QTable with 'ra' and 'dec' columns
qtable = QTable()
qtable['ra'] = [10.0, 20.0, 30.0]  # Example RA values in degrees
qtable['dec'] = [45.0, 50.0, 55.0]  # Example Dec values in degrees

# Current date in ISO format (YYYY-MM-DD)
current_date = '2024-04-23'

# Transform coordinates to the current date
transformed_qtable = transform_coordinates_to_current_date(qtable, current_date)

# Print the transformed QTable
print(transformed_qtable)

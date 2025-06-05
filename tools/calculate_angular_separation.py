from astropy.coordinates import SkyCoord
import astropy.units as u

def get_coordinate_input(label):
    print(f"\nEnter coordinates for {label}:")
    ra = float(input("  RA (decimal degrees): "))
    dec = float(input("  Dec (decimal degrees): "))
    return SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')

def main():
    print("Angular Separation Calculator\n")

    # Get user input for two coordinates
    coord1 = get_coordinate_input("Object 1")
    coord2 = get_coordinate_input("Object 2")

    # Compute separation
    separation = coord1.separation(coord2)

    # Display results
    print(f"\nAngular separation:")
    print(f"  {separation.degree:.6f} degrees")
    print(f"  {separation.arcminute:.6f} arcminutes")
    print(f"  {separation.arcsecond:.6f} arcseconds")

if __name__ == "__main__":
    main()

import pandas as pd
import numpy as np
import sys
import datetime

print('### Script start! ###')
print('# Timestamp: ' + str(datetime.datetime.now()))

# Input argument: path to fixed-width WDS star catalog text file
wdss_file = sys.argv[1]

# Define column metadata
column_names = [
    'WDSS identifier', 'Component identifier', 'First/last observation',
    'Number of astrometric observations', 'Position angle', 'Sepatarion',
    'Flag for separation units', 'Magnitude', 'G Filter', 'Infrared magnitude',
    'Filter', 'Spectral type', 'Proper motion', 'Parallax', 'Alternate name',
    'Note flags', 'Coord (RA)', 'Coord (DEC)', 'Designation in main WDS',
    'Discoverer designation', 'Component designation'
]
col_starts = [1, 15, 24, 29, 33, 37, 43, 45, 50, 52, 57, 59, 65, 82, 90, 115, 118, 127, 137, 148, 155]
col_ends = [14, 23, 28, 32, 36, 42, 44, 50, 51, 57, 58, 64, 81, 89, 114, 117, 126, 136, 147, 155, 160]
col_widths = [end - start + 1 for start, end in zip(col_starts, col_ends)]
print('col widths: ', col_widths)

# Read data
wdss_data = pd.read_fwf(wdss_file, widths=col_widths, names=column_names, dtype=str)
print(wdss_data.head())

# Function to convert RA hhmmss.ss to degrees
def ra_to_deg(ra_str):
    try:
        if pd.isna(ra_str):
            return np.nan
        ra_str = ra_str.strip()
        h = int(ra_str[0:2])
        m = int(ra_str[2:4])
        s = float(ra_str[4:])
        return (h + m/60 + s/3600) * 15
    except:
        return np.nan

# Function to convert DEC ddmmss.s to degrees
def dec_to_deg(dec_str):
    try:
        if pd.isna(dec_str):
            return np.nan
        dec_str = dec_str.strip()
        sign = -1 if dec_str[0] == '-' else 1
        d = int(dec_str[0:2])
        m = int(dec_str[2:4])
        s = float(dec_str[4:])
        return sign * (abs(d) + m/60 + s/3600)
    except:
        return np.nan

# Apply conversion
wdss_data["RA_deg"] = wdss_data["Coord (RA)"].apply(ra_to_deg)
wdss_data["DEC_deg"] = wdss_data["Coord (DEC)"].apply(dec_to_deg)

# Save to CSV
output_file = 'wdss_output.csv'
wdss_data.to_csv(output_file, index=False)
print(f'Data successfully saved to {output_file}')

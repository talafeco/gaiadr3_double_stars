# Működő kód konverzióra:


import pandas as pd
import numpy as np
import sys
import datetime

print('### Script start! ###')
print('# Timestamp: ' + str(datetime.datetime.now()))

# Input: path to fixed-width WDS catalog file
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
#col_starts = [1, 15, 24, 29, 33, 37, 43, 45, 50, 52, 57, 59, 65, 82, 90, 115, 116, 125, 134, 148, 155]
#col_ends = [14, 23, 28, 32, 36, 42, 44, 50, 51, 57, 58, 64, 81, 89, 114, 115, 125, 133, 147, 155, 160]
#col_widths = [end - start + 1 for start, end in zip(col_starts, col_ends)]
#print('col widths: ', col_widths)

col_specs = [(1,14), (15, 18), (24, 28), (30, 32), (33, 36), (38, 43), (44, 44), (46, 50), (51, 51), (53, 57), (58, 58), (59, 64), (66, 81), (83, 89), (90, 114), (116, 117), (118, 127), (127, 136), (137, 147), (148, 155), (155, 160)]

# Read fixed-width file
#wdss_data = pd.read_fwf(wdss_file, widths=col_widths, names=column_names, dtype=str)
wdss_data = pd.read_fwf(wdss_file, colspecs=col_specs, names=column_names, dtype=str)

print(wdss_data.head())

# Function to convert RA hhmmss.ss to degrees
def ra_to_deg(ra_str):
    try:
        if pd.isna(ra_str):
            return np.nan
        ra_str = ra_str.strip()

        # Fix short formats
        if len(ra_str) < 6:
            return np.nan

        # Pad right side if no decimal
        if '.' not in ra_str:
            ra_str += '.00'

        # Parse into components
        h = int(ra_str[0:2])
        m = int(ra_str[2:4])
        s = float(ra_str[4:])

        return (h + m / 60 + s / 3600) * 15
    except Exception as e:
        print(f"RA parse error: '{ra_str}' -> {e}")
        return np.nan

# Function to convert DEC ddmmss.s to degrees
def dec_to_deg(dec_str):
    try:
        if pd.isna(dec_str):
            return np.nan
        dec_str = dec_str.strip().zfill(8)
        sign = -1 if dec_str[0] == '-' else 1
        if dec_str[0] in '+-':
            d = int(dec_str[1:3])
            m = int(dec_str[3:5])
            s = float(dec_str[5:])
        else:
            d = int(dec_str[0:2])
            m = int(dec_str[2:4])
            s = float(dec_str[4:])
        return sign * (d + m / 60 + s / 3600)
    except Exception:
        return np.nan

# Apply conversion
wdss_data["RA_deg"] = wdss_data["Coord (RA)"].apply(ra_to_deg)
wdss_data["DEC_deg"] = wdss_data["Coord (DEC)"].apply(dec_to_deg)

# Save to CSV
output_file = 'wdss_output.csv'
wdss_data.to_csv(output_file, index=False)
print(f'Data successfully saved to {output_file}')




# Multicore kód:


'''import pandas as pd
import numpy as np
import sys
import datetime
from joblib import Parallel, delayed

print('### Script start! ###')
print('# Timestamp: ' + str(datetime.datetime.now()))

# Input: path to fixed-width WDS catalog file
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
col_starts = [1, 15, 24, 29, 33, 37, 43, 45, 50, 52, 57, 59, 65, 82, 90, 115, 116, 125, 134, 148, 155]
col_ends = [14, 23, 28, 32, 36, 42, 44, 50, 51, 57, 58, 64, 81, 89, 114, 115, 125, 133, 147, 155, 160]
col_widths = [end - start + 1 for start, end in zip(col_starts, col_ends)]
print('col widths: ', col_widths)

# Read fixed-width file
wdss_data = pd.read_fwf(wdss_file, widths=col_widths, names=column_names, dtype=str)
print(wdss_data.head())

# Function to convert RA hhmmss.ss to degrees
def ra_to_deg(ra_str):
    try:
        if pd.isna(ra_str):
            return np.nan
        ra_str = ra_str.strip()
        if len(ra_str) < 6:
            return np.nan
        if '.' not in ra_str:
            ra_str += '.00'
        h = int(ra_str[0:2])
        m = int(ra_str[2:4])
        s = float(ra_str[4:])
        if not (0 <= h < 24 and 0 <= m < 60 and 0 <= s < 60):
            return np.nan
        return (h + m / 60 + s / 3600) * 15
    except Exception as e:
        print(f"[RA ERROR] Could not parse: '{ra_str}' — {e}")
        return np.nan

# Function to convert DEC ddmmss.s to degrees
def dec_to_deg(dec_str):
    try:
        if pd.isna(dec_str):
            return np.nan
        dec_str = dec_str.strip()
        sign = -1 if dec_str.startswith("-") else 1
        dec_str = dec_str.lstrip("+-")
        if len(dec_str) < 6:
            return np.nan
        d = int(dec_str[0:2])
        m = int(dec_str[2:4])
        s = float(dec_str[4:])
        if not (0 <= m < 60 and 0 <= s < 60):
            return np.nan
        return sign * (d + m / 60 + s / 3600)
    except Exception as e:
        print(f"[DEC ERROR] Could not parse: '{dec_str}' — {e}")
        return np.nan

# Parallel processing
n_jobs = -1  # Use all available cores
#n_jobs = 2  # Use all available cores

wdss_data["RA_deg"] = Parallel(n_jobs=n_jobs)(delayed(ra_to_deg)(v) for v in wdss_data["Coord (RA)"])
wdss_data["DEC_deg"] = Parallel(n_jobs=n_jobs)(delayed(dec_to_deg)(v) for v in wdss_data["Coord (DEC)"])

# Save to CSV
output_file = 'wdss_output.csv'
wdss_data.to_csv(output_file, index=False)
print(f'Data successfully saved to {output_file}')'''

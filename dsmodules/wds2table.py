from astropy.table import Table, vstack, hstack
import numpy as np
from dscalculation import calculate_wds_ra_hourangle, calculate_wds_dec_hourangle, create_unique_id, delete_invalid_lines_wds

def create_wds_table(wds_file):
    wds_converters = {  '2000 Coord': np.str_, 'Discov': np.str_, 'Comp': np.str_, 'Date (first)': np.str_, 'Date (last)': np.str_, 'Obs': np.str_, 'PA_f': np.float64, 'PA_l': np.float64, 'Sep_f': np.float64, 'Sep_l': np.float64, 'Mag_A': np.str_, 'Mag_B': np.str_, 'Spectral A/B': np.str_, 'PM_A_ra': np.str_, 'PM_A_dec': np.str_, 'PM_B_ra': np.str_, 'PM_B_dec': np.str_, 'D. Number': np.str_, 'Notes': np.str_, 'Coord (RA)': np.str_, 'Coord (DEC)': np.str_}

    wds_data = Table.read(wds_file,
                        names=('2000 Coord', 'Discov', 'Comp', 'Date (first)', 'Date (last)', 'Obs',
                                'PA_f', 'PA_l', 'Sep_f', 'Sep_l', 'Mag_A',
                                'Mag_B', 'Spectral A/B', 'PM_A_ra', 'PM_A_dec',
                                'PM_B_ra', 'PM_B_dec', 'D. Number', 'Notes', 'Coord (RA)',
                                'Coord (DEC)'
                            ),
                        converters=wds_converters,
                        format='ascii.fixed_width',
                        header_start=2, data_start=5,
                        col_starts=(0, 10, 17, 23, 28, 33, 38, 42, 46, 52, 58, 64, 70, 80, 84, 89, 93, 98, 107, 112, 121),
                        col_ends=(9, 16, 21, 26, 31, 36, 40, 44, 50, 56, 61, 68, 78, 83, 87, 92, 96, 105, 110, 120, 129),
                        )

    wds_table = hstack([wds_data, calculate_wds_ra_hourangle(wds_data['Coord (RA)'])])
    wds_table.rename_column('col0', 'Coord (RA) hms')
    wds_table = hstack([wds_table, calculate_wds_dec_hourangle(wds_data['Coord (DEC)'])])
    wds_table.rename_column('col0', 'Coord (DEC) dms')
    wds_table = hstack([wds_table, create_unique_id(wds_data['2000 Coord'], wds_data['Discov'])])
    wds_table.rename_column('col0', 'Unique ID')
    wds_table = delete_invalid_lines_wds(wds_table)
    
    return wds_table
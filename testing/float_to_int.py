import numpy as np

def convert_to_int16(arr_float64):
    # Convert to int16
    arr_int16 = arr_float64.astype(np.int16)
    return arr_int16

# Example usage:
arr_float64 = np.array([1.5, 2.7, 3.9, 4.1])
arr_int16 = convert_to_int16(arr_float64)
print(arr_int16)
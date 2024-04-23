'''
This function subtract_and_clip takes two arguments: arr, the NumPy array, and num, the number to subtract from each element. It subtracts num from each element of arr, and then any elements that become negative are set to zero. Finally, it returns the resulting array.
'''

import numpy as np

def subtract_and_clip(arr, num):
    result = arr - num
    result[result < 0] = 0
    return result

# Example usage:
arr = np.array([5, 8, 10, 3, 1])
num = 3
result = subtract_and_clip(arr, num)
print(result)
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
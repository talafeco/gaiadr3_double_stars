'''
Ensure that your CSV file has headers named "Absolute Magnitude" and "B-V Index" (or adjust the code accordingly if the headers are different).

This code will generate a HR diagram where stars are plotted with their B-V index on the x-axis, absolute magnitude on the y-axis, and colored based on their B-V index. Brighter stars will appear towards the top of the plot due to the inverted y-axis.
'''


import matplotlib.pyplot as plt
import pandas as pd

# Load data from CSV file
data = pd.read_csv("stars_data.csv")

# Extracting data
abs_magnitude = data['Absolute Magnitude']
b_v_index = data['B-V Index']

# Define colors based on B-V index
colors = b_v_index

# Plot HR diagram
plt.figure(figsize=(10, 8))
plt.scatter(b_v_index, abs_magnitude, c=colors, cmap='coolwarm', edgecolors='none', alpha=0.8)
plt.title('Hertzsprung-Russell Diagram')
plt.xlabel('B-V Index')
plt.ylabel('Absolute Magnitude')
plt.gca().invert_yaxis()  # Invert y-axis to display brighter stars at the top
plt.colorbar(label='B-V Index')
plt.grid(True)
plt.show()

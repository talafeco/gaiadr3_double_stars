"""
Description
This tool calculates the most probable common distance of a binary star system (double star) by combining the parallax measurements of the two components.

How it works:
    - Input: It takes the parallax (in mas) and parallax error for two stars (Star A and Star B).
    - Probability Modeling: It treats each star's measurement as a Gaussian (Normal) probability distribution.
    - Joint Estimation: It multiplies the two probability distributions together. The resulting "peak" represents the parallax value where the two measurements most strongly agree.
    - Scaling: Since multiplying small probabilities results in tiny numbers, the code automatically scales the joint curve so it is visible on the same plot as the input stars.
    - Distance Calculation: It converts the most probable parallax into distance using the formula d(pc)=1000/π(mas).
    - Visualization: It produces a plot showing the two stars' measurement curves and the "consensus" joint probability curve.

How to use in CLI (Bash):
python common_distance.py <par_a> <err_a> <par_b> <err_b> -m "Designation of the main star" (optional) -c "Designation of the companion star" (optional)

How to use in other Python program:
Example
# analytical_script.py
from common_distance import CalculateDistance

# 1. Define your data
parallax1, error1 = 3.45, 0.1
parallax2, error2 = 3.30, 0.4

# 2. Initialize the class
# Note: Position arguments match *args, named arguments match **kwargs
calculator = CalculateDistance(p1, e1, p2, e2, star_a="Alpha Cen A", star_b="Alpha Cen B")

# 3. Access the results
print(f"Most Probable Parallax: {calculator.most_probable_parallax:.4f} mas")
print(f"Calculated Distance: {calculator.distance_in_parsec:.2f} pc")

# 4. Generate the plot
calculator.plot()
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import argparse

class CalculateDistance():
    """
    Finds the most probable common distance and plots it with 
    SCALED joint probability for better visibility.
    """
    def __init__(self, *args, **kwargs):
        self.define_star_names(star_a=kwargs["star_a"], star_b=kwargs["star_b"])
        self.par_a = args[0]
        self.err_a = args[1]
        self.par_b = args[2]
        self.err_b = args[3]
        
        # 1. Define Range
        self.start = min(self.par_a, self.par_b) - 4 * max(self.err_a, self.err_b)
        self.end = max(self.par_a, self.par_b) + 4 * max(self.err_a, self.err_b)
        
        # 2. Calculate X Axis
        self.x_axis = np.linspace(self.start, self.end, 1000)
        
        # 3. Calculate PDFs
        self.calc_pdf() 
        
        # 4. Calculate Joint Probability
        self.joint_pdf_raw = self.pdf_a * self.pdf_b
        
        # 5. Calculate Scale Factors
        self.max_height_inputs = max(np.max(self.pdf_a), np.max(self.pdf_b))
        self.scale_factor = self.calc_max_height_joint(self.joint_pdf_raw)
        self.joint_pdf_plot = self.joint_pdf_raw * self.scale_factor
        
        # 6. Find the Peak
        self.peak_index = np.argmax(self.joint_pdf_raw) 
        
        # 7. Auto-calculate distance (optional, but ensures attributes exist)
        self.calc_distance()

    def calc_max_height_joint(self, joint_pdf_raw):
        self.max_height_joint = np.max(self.joint_pdf_raw)
        # Avoid division by zero safety check
        if self.max_height_joint > 0:
           result  = self.max_height_inputs / self.max_height_joint
        else:
           result = 1

        return result

    def define_star_names(self, **kwargs):
        if kwargs["star_a"]:
            self.star_designation_a = kwargs["star_a"]
        else:
            self.star_designation_a = "Main star"

        if kwargs["star_b"]:
            self.star_designation_b = kwargs["star_b"]
        else:
            self.star_designation_b = "Companion"


    def calc_pdf(self):
        # Use the stored self.x_axis instead of recalculating
        self.pdf_a = norm.pdf(self.x_axis, loc=self.par_a, scale=self.err_a)
        self.pdf_b = norm.pdf(self.x_axis, loc=self.par_b, scale=self.err_b)

    def calc_distance(self):
        self.most_probable_parallax = self.x_axis[self.peak_index]
        # Avoid division by zero
        if self.most_probable_parallax != 0:
            self.distance_in_parsec = 1000.0 / self.most_probable_parallax
        else:
            self.distance_in_parsec = np.inf

    def plot(self):
        # 5. Plotting
        plt.figure(figsize=(10, 6))
        
        # Plot individual stars
        plt.plot(self.x_axis, self.pdf_a, label=self.star_designation_a, linestyle='--', alpha=0.7)
        plt.plot(self.x_axis, self.pdf_b, label=self.star_designation_b, linestyle='--', alpha=0.7)
        
        # Plot the SCALED joint distribution
        plt.plot(self.x_axis, self.joint_pdf_plot, label='Joint Probability (Scaled)', color='red', linewidth=2.5)
        plt.fill_between(self.x_axis, self.joint_pdf_plot, color='red', alpha=0.05)
        
        # Mark the peak
        plt.axvline(self.most_probable_parallax, color='black', linestyle=':', label=f'Peak: {self.most_probable_parallax:.3f} mas')
        
        plt.xlabel('Parallax (mas)')
        plt.ylabel('Probability Density (Scaled)')
        plt.title(f'Common Distance Estimation\nMost Probable Distance: {self.distance_in_parsec:.1f} pc')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()

# Define input argumenst (given when using command line)
parser = argparse.ArgumentParser()
parser.add_argument("parallax_a", help="Parallax of star A (mas)", type=float)
parser.add_argument("parallax_err_a", help="Parallax error of star A (mas)", type=float)
parser.add_argument("parallax_b", help="Parallax of star B (mas)", type=float)
parser.add_argument("parallax_err_b", help="Parallax error of star B (mas)", type=float)
parser.add_argument("-m", "--main_star", help="Designation of the main star", type=str)
parser.add_argument("-c", "--companion_star", help="Designation of the companion star", type=str)
args = parser.parse_args()

double_star_parallaxes = CalculateDistance(args.parallax_a, args.parallax_err_a, args.parallax_b, args.parallax_err_b, star_a=args.main_star, star_b=args.companion_star)
double_star_parallaxes.plot()

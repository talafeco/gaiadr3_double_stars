import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

class CalculateDistance:
    """
    Finds the most probable common distance and plots it with 
    SCALED joint probability for better visibility.
    """
    def __init__(self, par_a, err_a, par_b, err_b, star_a, star_b):
        self.star_designation_a = star_a
        self.star_designation_b = star_b
        self.par_a = par_a
        self.err_a = err_a
        self.par_b = par_b
        self.err_b = err_b
        
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
        self.max_height_joint = np.max(self.joint_pdf_raw)
        
        # Avoid division by zero safety check
        if self.max_height_joint > 0:
            self.scale_factor = self.max_height_inputs / self.max_height_joint
        else:
            self.scale_factor = 1

        self.joint_pdf_plot = self.joint_pdf_raw * self.scale_factor
        
        # 6. Find the Peak
        self.peak_index = np.argmax(self.joint_pdf_raw) 
        
        # 7. Auto-calculate distance (optional, but ensures attributes exist)
        self.calc_distance()

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

 
# Example Run
double_star_parallaxes = CalculateDistance(19.200, 0.404, 19.808, 0.769, 'Gaia DR3 133984292234128896', 'Gaia DR3 133984360953300736') # args: star A parallax, star A parallax error, star B parallax, star B parallax error, star designation A, star designation B
print(double_star_parallaxes.most_probable_parallax, double_star_parallaxes.distance_in_parsec)
double_star_parallaxes.plot()


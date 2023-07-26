'''

Using Sci py

'''
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import minimize_scalar


'''
Calculates the total energy of a system per unit volume
 of ferromagnetic material in an externally applied magnetic field H.

'''
def calculate_energy_sweep(H, Phi_values, J_AF, M, Ku, K1, anisotropy_axis):
    # Convert angles to radians
    Phi_rad = np.radians(Phi_values)
    
    # Calculate the energy terms
    phi_1 = Phi_rad[:-1]
    phi_2 = Phi_rad[1:]
    
    energy_H = -0.5 * H * M * (np.cos(phi_1) + np.cos(phi_2))
    energy_J_AF = -J_AF * np.cos(phi_1 - phi_2)
    
    ''' Anisotropic energy term '''
    
    energy_anis = (1/8) * K1 * (np.sin(2 * (phi_1 - anisotropy_axis)) ** 2 + np.sin(2 * (phi_2 - anisotropy_axis)) ** 2) # Cubic Symmetry [001] growth (see Folkerts)

    #energy_anis = 0.5 * Ku * (np.sin(Phi_rad[:-1] - anisotropy_axis) ** 2 + np.sin(Phi_rad[1:] - anisotropy_axis) ** 2) # Uniaxial
    
    energy_total = energy_H + energy_J_AF + energy_anis
    
    return energy_total


'''

Application 

'''

# Derivative of the energy with respect to the applied field angle
def energy_derivative(H, M):
    energy_values = calculate_energy_sweep(H, Phi_values, J_AF, M, Ku, K1, anisotropy_axis)
    energy_derivative = np.gradient(energy_values, Phi_values[:-1]) # Derivative  of energy with respect to applied field angle
    return -np.max(np.abs(energy_derivative))  # Return the negative maximum absolute derivative

##'''
##
##Material Properties
##
##'''
##J_AF = 0  # Antiferromagnetic coupling strength (small value)
##Ku = 0  # Anisotropy constant (uniaxial anisotropy)
##K1 = 0  # First-order anisotropy parameter (small value)
##anisotropy_axis = np.radians(0)  # Anisotropy axis angle in radians (uniaxial along y-axis)
##
##
##
### Applied angle sweep 'range'
##Phi_values = np.linspace(0, 360, 100)  # Array of applied field angles in degrees with 100 points
##
##
##'''
##
##Editable 'sweep 'Parameters
##
##'''
##
##H_values = np.linspace(-2, 2, 100)  # Array of applied magnetic field strengths 
##
### Define the range and number of points for magnetization
min_magnetization = 0  # Minimum magnetization value (Tesla)
max_magnetization = 0.5   # Maximum magnetization value (Tesla)
num_points_magnetization = 2  # Number of points for magnetization (same as H_values)
'''

Soft Parameters?

'''

# Parameters for Soft Magnetic Material
H_values = np.linspace(-1, 1, 100)  # Applied magnetic field strengths from -1 to 1 Tesla with 100 points
Phi_values = np.linspace(0, 360, 100)  # Applied field angles in degrees with 100 points
J_AF = 0.2  # Antiferromagnetic coupling strength (small value)
Ku = 0.1  # Anisotropy constant (low value)
K1 = 0.1  # First-order anisotropy parameter (small value)
anisotropy_axis = np.radians(90)  # Anisotropy axis angle at 90 degrees (no strong preference)


'''

Optimization

'''
# Calculate the step size for magnetization
magnetization_step_size = (max_magnetization - min_magnetization) / (num_points_magnetization - 1)

# Array to store the optimized H values for each M (both increasing and decreasing)
H_optimized_increasing = []
H_optimized_decreasing = []
M_values = []

# Loop through the M values to find the hysteresis loop for both increasing and decreasing magnetization
for i in range(num_points_magnetization):
    M = min_magnetization + i * magnetization_step_size
    M_values.append(M)  # Append the M value to the list

    # Find the optimized H for the given M using minimize_scalar on energy_derivative
    result_increasing = minimize_scalar(energy_derivative, args=(M,), bounds=(H_values[0], H_values[-1]))
    H_optimized_increasing.append(result_increasing.x)

    # Find the optimized H for the given M using minimize_scalar on energy_derivative with negative magnetization
    result_decreasing = minimize_scalar(energy_derivative, args=(-M,), bounds=(H_values[0], H_values[-1]))
    H_optimized_decreasing.append(result_decreasing.x)

# Concatenate the H_optimized_increasing and H_optimized_decreasing to create the full hysteresis loop
H_optimized = np.concatenate([H_optimized_increasing, H_optimized_decreasing[::-1]])

# Create an array of corresponding M values for the full hysteresis loop
M_values = np.concatenate([M_values, M_values[::-1]])

'''

Plot Hysteresis Loop

'''
plt.plot(H_optimized, M_values, marker='o')
#plt.plot(M_values,H_optimized , marker='o')
plt.xlabel('Optimized Applied Magnetic Field (H)')
plt.ylabel('Magnetization (M)')
plt.title('Hysteresis Loop')
plt.grid(True)
plt.show()


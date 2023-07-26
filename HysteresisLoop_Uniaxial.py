import numpy as np
import math
from scipy.optimize import minimize

# ... (previous code remains unchanged)

def calculate_energy(H, phi, J_AF, M, Ku, K1, anisotropy_axis):
    # Convert angle to radians
    phi_rad = np.radians(phi)
    
    # Calculate the energy terms
    energy_H = -0.5 * H * M * (np.cos(phi_rad) + np.cos(phi_rad + anisotropy_axis))
    energy_J_AF = -J_AF * np.cos(phi_rad - (phi_rad + anisotropy_axis))
    energy_anis = 0.5 * Ku * (np.sin(phi_rad - anisotropy_axis) ** 2)
    
    energy_total = energy_H + energy_J_AF + energy_anis
    return energy_total

# Function to find the local minimum energy for a given H
def find_local_min_energy(H):
    # Use minimize function to find the local minimum of the energy function
    result = minimize(calculate_energy, x0=0.0, args=(H, J_AF, M, Ku, K1, anisotropy_axis))
    
    # The minimum energy is in 'result.fun'
    min_energy = result.fun
    
    # The corresponding angle is in 'result.x[0]'
    min_angle_radians = result.x[0]
    
    # Convert angle from radians to degrees
    min_angle_degrees = np.degrees(min_angle_radians)
    
    # Ensure the angle is within [0, 360] degrees
    min_angle_degrees = min_angle_degrees % 360
    
    return min_angle_degrees, min_energy

# Example usage
H_values = np.linspace(-1.5, 1.5, 100)  # Array of applied magnetic field strengths (both positive and negative)
J_AF = 1  # Antiferromagnetic coupling strength
M = 1  # Magnitude of magnetization vectors in layers A in Tesla
Ku = 1  # Uniaxial anisotropy constant
K1 = 1  # First-order anisotropy parameter
anisotropy_axis = np.radians(45)  # Anisotropy axis angle in radians

# Initialize lists to store local minimum angles and energies
local_minima_angles = []
local_minima_energies = []

# Find and store the local minimums for each value of H
for H in H_values:
    min_angle, min_energy = find_local_min_energy(H)
    local_minima_angles.append(min_angle)
    local_minima_energies.append(min_energy)



# Print the local minimums
print("Local Minima:")
for angle, energy in zip(local_minima_angles, local_minima_energies):
    angle_degrees = np.degrees(angle)  # Convert angle from radians to degrees for printing
    print(f"Angle: {angle_degrees:.2f} degrees, Energy: {energy:.2f}")


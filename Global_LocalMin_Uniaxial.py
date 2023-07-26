import numpy as np
import math
import matplotlib.pyplot as plt

def calculate_energy_sweep(H, Phi_values, J_AF, M, Ku, K1, anisotropy_axis):
    # Convert angles to radians
    Phi_rad = np.radians(Phi_values)
    
    # Calculate the energy terms
    phi_1 = Phi_rad[:-1]
    phi_2 = Phi_rad[1:]
    
    energy_H = -0.5 * H * M * (np.cos(phi_1) + np.cos(phi_2))
    energy_J_AF = -J_AF * np.cos(phi_1 - phi_2)
    
    # Uniaxial Term
    energy_anis = 0.5 * Ku * (np.sin(Phi_rad[:-1] - anisotropy_axis) ** 2 + np.sin(Phi_rad[1:] - anisotropy_axis) ** 2)
    energy_total = energy_H + energy_J_AF + energy_anis
    
    return energy_total

# Example usage
#H_values = [-1.5,-1.0,-0.75, -0.5, 0, 0.5, 0.75, 1.0, 1.5]  # Array of applied magnetic field strengths
H_values = [ 0, 0.5, 0.75, 1.0, 1.5]  # Array of applied magnetic field strengths
Phi_values = np.linspace(0, 360, 100)  # Array of applied field angles in degrees
J_AF = 0  # Antiferromagnetic coupling strength
M = 1  # Magnitude of magnetization vectors in layers A in Tesla
Ku = 1  # Uniaxial anisotropy constant
K1 = 1  # First-order anisotropy parameter
anisotropy_axis = np.radians(45)  # Anisotropy axis angle in radians

# Lists to store local and global minimum angles for each H value
local_min_angles_per_H = []
global_min_angles_per_H = []

plt.figure()  # Create a new figure

for H in H_values:
    energy_values = calculate_energy_sweep(H, Phi_values, J_AF, M, Ku, K1, anisotropy_axis)

    # Find and track local minimum angle for the current H
    min_energy = np.min(energy_values)
    min_angle = Phi_values[:-1][np.argmin(energy_values)]
    local_min_angles_per_H.append(min_angle)

    # Check if the current minimum energy is less than the global minimum energy for the current H
    if min_energy < global_min_energy:
        global_min_energy = min_energy
        global_min_angle = min_angle

    # Track the global minimum angle for the current H
    global_min_angles_per_H.append(global_min_angle)

# Calculate cosine of the local minimum angles and store them in arrays
for angle, energy in local_minima:
    cosine_angle = np.cos(np.radians(angle))
    local_minima_angles = np.append(local_minima_angles, angle)
    local_minima_cosines = np.append(local_minima_cosines, cosine_angle)

     # Find and track global minimum angle and its cosine
    min_energy = np.min(energy_values)
    min_angle = Phi_values[:-1][np.argmin(energy_values)]
    min_cosine = np.cos(np.radians(min_angle))
    
    global_min_angles = np.append(global_min_angles, min_angle)
    global_min_cosines = np.append(global_min_cosines, min_cosine)
    global_min_energy_values = np.append(global_min_energy_values, min_energy)

    print(global_min_angle)

# New plot for cosine of global minimum angles vs. applied field
plt.figure()
plt.plot(H_values, np.cos(np.radians(global_min_angles_per_H)), 'o-', label='Cosine of Global Minimum Angle')
plt.xlabel('Applied Magnetic Field Strength (H)')
plt.ylabel('Cosine of Global Minimum Angle')
plt.title('Cosine of Global Minimum Angle vs. Applied Magnetic Field')
plt.legend()
plt.grid(True)

plt.ylim(-1, 1)
plt.show()

##print("Global Minimum Energy:", global_min_energy)
##print("Global Minimum Angle:", global_min_angle)
##print("Local Minima:")
##for angle, energy in local_minima:
##    print(f"Angle: {angle:.2f} degrees, Energy: {energy:.2f}")

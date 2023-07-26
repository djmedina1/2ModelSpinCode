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
H_values = [0, 0.75, 1.0, 1.5]  # Array of applied magnetic field strengths
Phi_values = np.linspace(0, 360, 100)  # Array of applied field angles in degrees
J_AF = 1  # Antiferromagnetic coupling strength
M = 0.2  # Magnitude of magnetization vectors in layers A in Tesla
Ku = 1  # Uniaxial anisotropy constant
K1 = 1  # First-order anisotropy parameter
anisotropy_axis = np.radians(45)  # Anisotropy axis angle in radians

# Initialize variables to keep track of global and local minima
global_min_energy = float('inf')
global_min_angle = None
local_minima = []

plt.figure()  # Create a new figure

# Calculate and plot energy for each value of H
for H in H_values:
    energy_values = calculate_energy_sweep(H, Phi_values, J_AF, M, Ku, K1, anisotropy_axis)
    plt.plot(Phi_values[:-1], energy_values, label=f'H={H}')  # Plot energy values with a label
    
    # Find and track global minimum energy and angle
    min_energy = np.min(energy_values)
    min_angle = Phi_values[:-1][np.argmin(energy_values)]
    
    if min_energy < global_min_energy:
        global_min_energy = min_energy
        global_min_angle = min_angle
    
    # Find and track local minima
    local_minima_indices = np.where(energy_values == min_energy)[0]
    local_minima_angles = Phi_values[:-1][local_minima_indices]
    local_minima.extend(list(zip(local_minima_angles, [min_energy]*len(local_minima_indices))))

# Set custom x-axis tick locations and labels
x_ticks = np.arange(0, 361, 45) 
plt.xticks(x_ticks)

plt.xlabel('Applied Field Angle (degrees)')
plt.ylabel('Energy')
plt.title('Energy vs. Applied Field Angle for different H values')
plt.legend()  # Show the legend
plt.show()

print("Global Minimum Energy:", global_min_energy)
print("Global Minimum Angle:", global_min_angle)
print("Local Minima:")
for angle, energy in local_minima:
    print(f"Angle: {angle:.2f} degrees, Energy: {energy:.2f}")

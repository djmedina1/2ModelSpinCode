import numpy as np
import math
import matplotlib.pyplot as plt

def find_energy_minima(Phi_values, energy_values):
    # Find the indices where the derivative of energy with respect to Phi_values changes sign
    dE_dPhi_sign_changes = np.where(np.diff(np.sign(np.gradient(energy_values))) != 0)[0]

    # Indices of potential critical points (where dE_dPhi is zero)
    candidate_indices = [i for i in dE_dPhi_sign_changes if i > 0 and i < len(Phi_values) - 1]

    # Check each candidate index to determine if it's a local minimum
    minima_angles = []
    for index in candidate_indices:
        if energy_values[index] < energy_values[index - 1] and energy_values[index] < energy_values[index + 1]:
            # We found a local minimum
            minima_angles.append(Phi_values[index])

    return minima_angles

def calculate_energy_sweep(H, Phi_values, J_AF, M, Ku, K1, anisotropy_axis):
    # Convert angles to radians
    Phi_rad = np.radians(Phi_values)

    # Calculate the energy terms
    phi_1 = Phi_rad[:-1]
    phi_2 = Phi_rad[1:] # Using both phi's => two spin
    
    energy_H = -1 * H * M * (np.cos(phi_1))
    energy_J_AF = -J_AF * np.cos(phi_1 - phi_2)
    
    ''' Uniaxial Term ''' 
    #energy_anis = (1/8) * K1 * (np.sin(2 * (phi_1 - anisotropy_axis)) ** 2 + np.sin(2 * (phi_2 - anisotropy_axis)) ** 2)
    energy_anis = np.sin(Phi_rad[1:] - anisotropy_axis) ** 2
    energy_total = energy_H + energy_J_AF + energy_anis

    # Calculate the derivative of energy with respect to H (approximated using finite differences)
    dH = H_values[1] - H_values[0]
    dE_dH = np.gradient(energy_total, dH, edge_order=2)

    # Calculate magnetization M using the derivative of energy with respect to H
    magnetization = -dE_dH

    # Calculate the derivative of energy with respect to Phi_values (approximated using finite differences)
    dPhi = Phi_values[1] - Phi_values[0]
    dE_dPhi = np.gradient(energy_total, dPhi, edge_order=2)
    
    return energy_total, magnetization, dE_dPhi




# Example usage (unchanged)
H_values = [-1.5, -1.0, -0.75, 0, 0.75, 1.0, 1.5]
Phi_values = np.linspace(0, 360, 100)
J_AF = 0 # Set to zero for now.
M = 1
Ku = 1
K1 = 1
anisotropy_axis = np.radians(45)
'----'

# Lists to store local and global minimum angles for each H value
local_min_angles_per_H = []
global_min_angles_per_H = []

plt.figure()  # Create a new figure for the energy plot

# Calculate and plot energy for each value of H
for H in H_values:
    energy_values = calculate_energy_sweep(H, Phi_values, J_AF, M, Ku, K1, anisotropy_axis)
    plt.plot(Phi_values[:-1], energy_values, label=f'H={H}')  # Plot energy values with a label

    # Find and track local minimum angle for the current H
    min_energy = np.min(energy_values)
    min_angle = Phi_values[:-1][np.argmin(energy_values)]
    local_min_angles_per_H.append(min_angle)

    # Check if the current minimum energy is less than the global minimum energy for the current H
    if min_energy < global_min_energy:
        global_min_energy = min_energy
        global_min_angle = min_angle

    # Track the global minimum angle for the current H
    global_min_angles_per_H.append(global_min_an
'---'

# Set custom x-axis tick locations and labels
x_ticks = np.arange(0, 361, 45) 
plt.xticks(x_ticks)    

plt.xlabel('Applied Field Angle (degrees)')
plt.ylabel('Energy')
plt.title('Energy vs. Applied Field Angle for different H values')
plt.legend()  # Show the legend
plt.show()

# New plot for cosine of global minimum angles vs. applied field
plt.figure()
plt.plot(H_values, np.cos(np.radians(global_min_angles_per_H)), 'o-', label='Cosine of Global Minimum Angle')
plt.xlabel('Applied Magnetic Field Strength (H)')
plt.ylabel('Cosine of Global Minimum Angle')
plt.title('Cosine of Global Minimum Angle vs. Applied Magnetic Field')
plt.legend()
plt.grid(True)

# Set y-axis limits to -1 and 1
plt.ylim(-1, 1)

plt.show()
                


# Arrays to accumulate energy and magnetization values for each H
energy_values_list = []
magnetization_values_list = []
dE_dPhi_list = []

for H in H_values:
    energy_values, magnetization_values, dE_dPhi = calculate_energy_sweep(H, Phi_values, J_AF, M, Ku, K1, anisotropy_axis)
    energy_values_list.append(energy_values)
    magnetization_values_list.append(magnetization_values)
    dE_dPhi_list.append(dE_dPhi)

# Plot energy values for each H
plt.figure(figsize=(10, 5))
for i, H in enumerate(H_values):
    energy_values = energy_values_list[i]
    plt.plot(Phi_values[:-1], energy_values, label=f'H={H}')

plt.ylabel('Energy')
plt.xlabel('Applied Field Angle θ (deg)')
plt.grid(True)

# Plot derivative of energy with respect to applied field angle
plt.figure(figsize=(6, 4))
for i, H in enumerate(H_values):
    dE_dPhi = dE_dPhi_list[i]
    plt.plot(Phi_values[:-1], dE_dPhi, label=f'H={H}')
plt.ylabel('dE_dθ')
plt.xlabel('Applied Field Angle θ(deg)')
plt.grid(True)



plt.show()




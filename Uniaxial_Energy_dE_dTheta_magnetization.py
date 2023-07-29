import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.signal import find_peaks


def calculate_energy_sweep(H, Phi_values, J_AF, M, Ku, K1, anisotropy_axis):
    # Convert angles to radians
    Phi_rad = np.radians(Phi_values)
    
    # Calculate the energy terms
    phi_1 = Phi_rad[:-1]
    phi_2 = Phi_rad[1:]
    energy_H = -1 * H * M * (np.cos(Phi_rad[:len(Phi_rad)-1]))
    #energy_J_AF = -J_AF * np.cos(phi_1 - phi_2)
    
    ''' Uniaxial Term ''' 
    #energy_anis = (1/8) * K1 * (np.sin(2 * (phi_1 - anisotropy_axis)) ** 2 + np.sin(2 * (phi_2 - anisotropy_axis)) ** 2)
    energy_anis = Ku * np.sin(Phi_rad[:len(Phi_rad)-1] - anisotropy_axis) ** 2
    #debug
    energy_total = energy_H + 0 + energy_anis
    

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
H_values = [-1.5, -1.0, -0.75, -0.5, 0, 0.5, 0.75, 1.0, 1.5]
Phi_values = np.linspace(0, 361,1000)
J_AF = 0
M = 1
Ku = 1
K1 = 1
anisotropy_axis = np.radians(45)


# Arrays to accumulate energy and magnetization values for each H
energy_values_list = []
magnetization_values_list = []
dE_dPhi_list = []
easy_axis_min_angles_list = []   # To store the angles of the easy axis minimum for different H values

for H in H_values:
    energy_values, magnetization_values, dE_dPhi = calculate_energy_sweep(H, Phi_values, J_AF, M, Ku, K1, anisotropy_axis)
    energy_values_list.append(energy_values)
    magnetization_values_list.append(magnetization_values)
    dE_dPhi_list.append(dE_dPhi)

    # Find peaks (minima) in the energy values array for each H
    min_energy_indices, _ = find_peaks(-energy_values)
    if len(min_energy_indices) > 0:
        min_energy_index = min_energy_indices[0]
        easy_axis_min_angle = Phi_values[:-1][min_energy_index]
        easy_axis_min_angles_list.append(easy_axis_min_angle)
    else:
        easy_axis_min_angles_list.append(None)

# Plot energy values for each H
plt.figure(figsize=(10, 5))
for i, H in enumerate(H_values):
    energy_values = energy_values_list[i]
    label_text = f'H={H}'
    plt.plot(Phi_values[:-1], energy_values, label=label_text)

# Mark the easy axis minimum with a special marker for each H value
for i, H in enumerate(H_values):
    easy_axis_min_angle = easy_axis_min_angles_list[i]
    if easy_axis_min_angle is not None:
        plt.scatter(easy_axis_min_angle, energy_values_list[i][np.argmin(energy_values_list[i])], color='red', marker='o', label=f'Easy Axis Minimum (H={H})')

plt.ylabel('Energy')
plt.xlabel('Applied Field Angle θ (deg)')
plt.legend()
plt.grid(True)

# Print the angles of the easy axis minimum for different H values
print("Angles of the easy axis minimum for different H values:")
for i, H in enumerate(H_values):
    easy_axis_min_angle = easy_axis_min_angles_list[i]
    if easy_axis_min_angle is not None:
        print(f"H={H}: {easy_axis_min_angle:.2f} degrees")
    else:
        print(f"H={H}: No easy axis minimum found.")

# Convert angles to cosine of angles
cos_easy_axis_min_angles_list = [np.cos(np.radians(angle)) if angle is not None else None for angle in easy_axis_min_angles_list]

reflected_cos_easy_axis_min_angles_list = cos_easy_axis_min_angles_list[::-1]

# Plot cos(easy_axis_min_angles_list) vs H_values
plt.figure(figsize=(8, 5))
plt.scatter(H_values, cos_easy_axis_min_angles_list, color='red', marker='o')
plt.xlabel('Applied Field Strength (H)')
plt.ylabel('Cosine of Easy Axis Minimum Angle')
plt.title('Cosine of Easy Axis Minimum Angle vs Applied Field Strength')
plt.grid(True)
plt.show()

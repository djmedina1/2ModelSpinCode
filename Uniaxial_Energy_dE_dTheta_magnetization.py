import numpy as np
import math
import matplotlib.pyplot as plt

def calculate_energy_sweep(H, Phi_values, J_AF, M, Ku, K1, anisotropy_axis):
    # Convert angles to radians
    Phi_rad = np.radians(Phi_values)
    
    # Calculate the energy terms
    phi_1 = Phi_rad[:-1]
    phi_2 = Phi_rad[1:]
    
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
H_values = [0, 0.75, 1.0, 1.5]
Phi_values = np.linspace(0, 360, 100)
J_AF = 0
M = 1
Ku = 1
K1 = 1
anisotropy_axis = np.radians(45)

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
plt.legend()
plt.grid(True)

# Plot derivative of energy with respect to applied field angle
plt.figure(figsize=(6, 4))
for i, H in enumerate(H_values):
    dE_dPhi = dE_dPhi_list[i]
    plt.plot(Phi_values[:-1], dE_dPhi, label=f'H={H}')
plt.ylabel('dE_dθ')
plt.xlabel('Applied Field Angle θ(deg)')
plt.legend()
plt.grid(True)

    
# Plot magnetization vs applied field strength
magnetization_array = np.array(magnetization_values_list)
plt.figure(figsize=(6, 4))
plt.plot(H_values, magnetization_array[:, 0], marker='o')
plt.xlabel('Applied Field Strength (H)')
plt.ylabel('Magnetization (M)')
plt.title('Magnetization vs Applied Field Strength')
plt.grid(True)

plt.show()



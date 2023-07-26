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
H_values = [0, 0.5, 0.75, 1.0, 1.5]
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

# Arrays to store global and local minimum energy values and corresponding angles for each H
global_min_energy_values = []
global_min_angles = []
local_min_energy_values = []
local_min_angles = []

for H in H_values:
    energy_values, magnetization_values, dE_dPhi = calculate_energy_sweep(H, Phi_values, J_AF, M, Ku, K1, anisotropy_axis)
    energy_values_list.append(energy_values)
    magnetization_values_list.append(magnetization_values)
    dE_dPhi_list.append(dE_dPhi)

    # Find the global minimum for each H value
    global_min_energy = np.min(energy_values)
    global_min_angle = Phi_values[np.argmin(energy_values)]
    global_min_energy_values.append(global_min_energy)
    global_min_angles.append(global_min_angle)

    # Find the local minimums for each H value
    local_min_energies = []
    local_min_angles = []
    for i in range(1, len(energy_values) - 1):
        if energy_values[i] < energy_values[i - 1] and energy_values[i] < energy_values[i + 1]:
            local_min_energies.append(energy_values[i])
            local_min_angles.append(Phi_values[i])

    local_min_energy_values.append(local_min_energies)
    local_min_angles.append(local_min_angles)
    

    print(f"For H={H}:")
    print(f"Global Minimum: Energy = {global_min_energy}, Angle = {global_min_angle:.2f} degrees")
    print(f"Local Minimums: Energy = {local_min_energies}, Angles = {local_min_angles}")
    print()

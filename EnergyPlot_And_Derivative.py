import numpy as np
import math
import matplotlib.pyplot as plt

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


### Example usage ###

H_values = [0, 0.75, 1.0, 1.5]  # Array of applied magnetic field strengths
Phi_values = np.linspace(0, 360, 100)  # Array of applied field angles in degrees
J_AF = 1  # Antiferromagnetic coupling strength
M = 0.2  # Magnitude of magnetization vectors in layers A (Ferromagnetic layer) in Tesla
Ku = 1  # Anisotropy constant
K1 = 1  # First-order anisotropy parameter
anisotropy_axis = np.radians(45)  # Anisotropy axis angle in radians

fig, ax = plt.subplots(2)

# Calculate and plot energy for each value of H
for H in H_values:
    energy_values = calculate_energy_sweep(H, Phi_values, J_AF, M, Ku, K1, anisotropy_axis)
    energy_derivative = np.gradient(energy_values, Phi_values[:-1])
    ax[0].plot(Phi_values[:-1], energy_values, label=f'H={H}')  # Plot energy values with a label
    ax[1].plot(Phi_values[:-1], energy_derivative, label=f'H={H} T')  # Plot derivative values with a label


# Set custom x-axis tick locations and labels
x_ticks = np.arange(0, 361, 45) 
plt.xticks(x_ticks)

plt.xlabel('Applied Field Angle (degrees)')
plt.ylabel('Energy')
ax[0].set_title('Energy vs. Applied Field Angle')
ax[1].set_title('dE/dθ vs. Applied Field Angle')
plt.legend()  # Show the legend
plt.show()

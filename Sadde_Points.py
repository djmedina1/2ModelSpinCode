import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import minimize



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

def calculate_energy(H, Phi, J_AF, M, Ku, K1, anisotropy_axis):
    energy = K1 * (np.cos(Phi - anisotropy_axis) ** 2) - H * M
    return energy

# Example usage (unchanged)
H_values = [-1.5, -1.0, -0.75, -0.5, 0, 0.5, 0.75, 1.0, 1.5]
Phi_values = np.linspace(0, 361,1000)
J_AF = 0
M = 1
Ku = 1
K1 = 1
anisotropy_axis = np.radians(45)


# List to keep track of saddle points
saddle_points = []

for H in H_values:
    # Define the energy function for a given H, keeping other parameters constant
    energy_func = lambda Phi: calculate_energy(H, Phi, J_AF, M, Ku, K1, anisotropy_axis)

    # Find the critical points using numerical optimization
    result = minimize(energy_func, x0=1.57)  # Start the optimization from Phi=0.0 degrees

    # Check if the optimization was successful
    if result.success and result.fun is not None:
        Phi_critical = result.x[0]
        energy_at_critical = result.fun

        # Calculate the Hessian matrix to classify critical point
        hessian_matrix = np.gradient(np.gradient(energy_func(Phi_values)), Phi_values)
        hessian_matrix = np.diag(hessian_matrix)  # Convert to a diagonal matrix

        eigenvalues = np.linalg.eigvals(hessian_matrix)

        if np.any(np.real(eigenvalues) > 0) and np.any(np.real(eigenvalues) < 0):
            # Critical point has both positive and negative eigenvalues, indicating a saddle point
            saddle_points.append((H, np.degrees(Phi_critical), energy_at_critical))

# Print saddle points
print("Saddle points (H, Phi, Energy):")
for point in saddle_points:
    print(point)

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


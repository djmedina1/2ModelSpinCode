import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.signal import find_peaks


def calculate_energy_sweep(H, Phi_values, J_AF, M, Ku, K1, anisotropy_axis):
    # Convert angles to radians
    Phi_rad = np.radians(Phi_values)
    
    energy_H = 1 * H * M * (np.cos(Phi_rad[:len(Phi_rad)-1]))
    #energy_J_AF = -J_AF * np.cos(phi_1 - phi_2)
    
    ''' Uniaxial Term ''' 
    #energy_anis = (1/8) * K1 * (np.sin(2 * (phi_1 - anisotropy_axis)) ** 2 + np.sin(2 * (phi_2 - anisotropy_axis)) ** 2)
    energy_anis = Ku * np.sin(Phi_rad[:len(Phi_rad)-1] - anisotropy_axis) ** 2
    energy_total = energy_H + 0 + energy_anis

    # Calculate the derivative of energy with respect to Phi_values (approximated using finite differences)
    dPhi = Phi_values[1] - Phi_values[0]
    dE_dPhi = np.gradient(energy_total, dPhi, edge_order=2)
    d2E_dPhi2 = np.gradient(dE_dPhi, dPhi, edge_order=2)

    return energy_total, dE_dPhi, d2E_dPhi2  # Include the second derivative in the return statement

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
dE_dPhi_list = []
easy_axis_min_angles_list = []   # To store the angles of the easy axis minimum for different H values
d2E_dPhi2_list = []

# Arrays to accumulate points of inflection for each H and track minimum found status
inflection_points_list = []
minimum_found_list = [False] * len(H_values)
inflection_points_list = []
minimum_found_list = [False] * len(H_values)
inflection_included_list = [False] * len(H_values)  # To track if an inflection point is already included
inflection_found = False


for i, H in enumerate(H_values):
    energy_values, dE_dPhi, d2E_dPhi2 = calculate_energy_sweep(H, Phi_values, J_AF, M, Ku, K1, anisotropy_axis)
    energy_values_list.append(energy_values)

    dE_dPhi_list.append(dE_dPhi)
    d2E_dPhi2_list.append(d2E_dPhi2)

    # Find the points of inflection using the sign change in the second derivative
    inflection_indices = np.where(np.diff(np.sign(d2E_dPhi2)))[0]
    inflection_points = Phi_values[:-1][inflection_indices]

    # Append inflection points only if a minimum has not been found yet for this H value
    if not minimum_found_list[i]:
        energy_min_indices, _ = find_peaks(-energy_values)  # Negative energy values to find minima
        if len(energy_min_indices) > 0:
            min_energy_index = energy_min_indices[0]
            easy_axis_min_angle = Phi_values[:-1][min_energy_index]
            # Mark that a minimum has been found for this H value, so no more inflection points will be considered
            minimum_found_list[i] = True

            # Check if any inflection points were found before the minimum and at or past 90 degrees (pi/2 radians)
            valid_inflection_points = [angle for angle in inflection_points if round(angle) >= 90 and angle < easy_axis_min_angle and not inflection_found]
            if len(valid_inflection_points) > 0:
                # Include the first valid inflection point that occurs (if not already included)
                for inflection_point in valid_inflection_points:
                    if not inflection_included_list[i]:
                        inflection_found = True
                        easy_axis_min_angles_list.append([easy_axis_min_angle - 180, inflection_point])
                        inflection_included_list[i] = True

                # Include the minimum
                #easy_axis_min_angles_list.append(easy_axis_min_angle)
                
            else:
                # If no valid inflection points were found, add only the minimum to the list
                easy_axis_min_angles_list.append(easy_axis_min_angle)
        else:
            # If no minimum is found, check if any inflection point is available and include it
            if len(inflection_points) > 0:
                # Check if any inflection points were found at or past 90 degrees (pi/2 radians)
                valid_inflection_points = [angle for angle in inflection_points if round(angle) >= 90]
                if len(valid_inflection_points) > 0:
                    # Include the first valid inflection point (occurs at or past 90 degrees) if not already included
                    if not inflection_included_list[i]:
                        easy_axis_min_angles_list.append(valid_inflection_points[0])
                        inflection_included_list[i] = True
                    # No need to mark the minimum found here since we are including an inflection point
                else:
                    easy_axis_min_angles_list.append(None)
            else:
                easy_axis_min_angles_list.append(None)

    inflection_points_list.append(inflection_points)

### Print the angles of the minimum easy axis and the corresponding field strengths
##print("Angles of the easy axis minimum(s) and the corresponding fields:")
##angle_counter = 0
##for i, H in enumerate(H_values):
##    easy_axis_min_angles = easy_axis_min_angles_list[i]
##    if easy_axis_min_angles is not None:
##        if isinstance(easy_axis_min_angles, list):
##            # If it's an array of values, print each angle with the corresponding field strength
##            for angle in easy_axis_min_angles:
##                print(f"H={H}: {angle:.2f} degrees")
##        else:
##            # If it's a single value, print the angle along with the corresponding field strength
##            print(f"H={H}: {easy_axis_min_angles:.2f} degrees")
##    else:
##        print(f"H={H}: No minimum or valid inflection point found.")

# Accumulate the easy axis minimum angles and corresponding field strengths
min_angle_degrees_list = []
field_strengths_list = []
for i, H in enumerate(H_values):
    easy_axis_min_angles = easy_axis_min_angles_list[i]
    if easy_axis_min_angles is not None:
        # Convert the angle from radians to degrees
        easy_axis_min_angle_deg = easy_axis_min_angles
        # Append the easy axis minimum angles and the corresponding field strength to the lists
        if isinstance(easy_axis_min_angles, list):
            # Handle the case of multiple angles for the same field strength
            for angle in easy_axis_min_angles:
                min_angle_degrees_list.append(np.cos(np.radians(angle)))
                field_strengths_list.append(H)
        else:
            min_angle_degrees_list.append(-1 * np.cos(np.radians(easy_axis_min_angle_deg)))
            field_strengths_list.append(H)


# Convert the field_strengths_list to a NumPy array
field_strengths_array = np.array(field_strengths_list, dtype=float)

# Initialize an empty list to store the angles for which global minimums occur
global_min_angles = []

for i, H in enumerate(H_values):
    energy_values = energy_values_list[i]

    # Find the index of the minimum energy value for this H value
    min_energy_index = np.argmin(energy_values)

    # Get the corresponding angle for the global minimum
    global_min_angle = Phi_values[min_energy_index]

    # Since the energy is periodic, we need to check if there's a better minimum at the periodic boundary
    if min_energy_index > 0 and min_energy_index < len(Phi_values) - 1:
        if energy_values[min_energy_index] > energy_values[min_energy_index - 1] or energy_values[min_energy_index] > energy_values[min_energy_index + 1]:
            # If the minimum at the boundary is lower, update the global_min_angle
            if energy_values[min_energy_index - 1] < energy_values[min_energy_index + 1]:
                global_min_angle = Phi_values[min_energy_index - 1]
            else:
                global_min_angle = Phi_values[min_energy_index + 1]

    # Append the angle to the list of global minimum angles
    global_min_angles.append(global_min_angle)

# Convert the global_min_angles list to a NumPy array
global_min_angles_array = np.array(global_min_angles)


print("Angles for which global minimums occur for each applied field strength:")
print(global_min_angles_array)
# Create a scatter plot of easy axis minimum angles vs. applied field strength
plt.figure(figsize=(8, 5))
plt.scatter(field_strengths_array, min_angle_degrees_list, color='red', marker='o')

# Convert H_values to a NumPy array
H_values_array = np.array(H_values)

# Separate global minimum angles based on H values
global_min_angles_array_pos = global_min_angles_array[H_values_array >= 0]
global_min_angles_array_neg = global_min_angles_array[H_values_array < 0]

# Convert the field_strengths_list and min_angle_degrees_list to NumPy arrays
field_strengths_array = np.array(field_strengths_list, dtype=float)
min_angle_degrees_array = np.array(min_angle_degrees_list)

# Plot easy axis minimum angles vs. applied field strength
plt.figure(figsize=(8, 5))
plt.scatter(field_strengths_array, min_angle_degrees_array, color='red', marker='o', label='Easy Axis Minimums')

# Plot global minimum angles for H values >= 0
plt.scatter(H_values_array[H_values_array >= 0], np.cos(np.radians(global_min_angles_array_pos - 180)), color='blue', marker='x', label='Global Minimums (H >= 0)')

# Plot global minimum angles for H values < 0
#plt.scatter(H_values_array[H_values_array < 0], np.cos(np.radians(global_min_angles_array_neg)), color='green', marker='x', label='Global Minimums (H < 0)')

plt.xlabel('Applied Field Strength (H)')
plt.ylabel('Easy Axis Minimum Angle (degrees)')
plt.title('Easy Axis Minimum Angle vs Applied Field Strength')
plt.ylim(-1, 1)  # Set the y-axis limits to -1 and 1
plt.grid(True)
plt.show()

        


#Plot energy_values for each H
plt.figure(figsize=(10, 5))
for i, H in enumerate(H_values):
    energy_values = energy_values_list[i]
    label_text = f'H={H}'
    plt.plot(Phi_values[:len(Phi_values)-1],energy_values, label=label_text)
# Add the legend to the plot
plt.legend()
plt.show()


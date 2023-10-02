import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

def calculate_energy_sweep(H, Phi_values, J_AF, M, Ku, K1, anisotropy_axis):
    # Convert angles to radians
    Phi_rad = np.radians(Phi_values)

    phi_1 = Phi_rad[:]
    phi_2 = Phi_rad[::-1]

    
    energy_H = -1 * H * M * (np.cos(phi_1) + np.cos(phi_2))
    
    energy_J_AF = -J_AF * np.cos(phi_2 - phi_1) #coupling term
    
    ''' Uniaxial Term ''' 
    #energy_anis = (1/8) * K1 * (np.sin(2 * (phi_1 - anisotropy_axis)) ** 2 + np.sin(2 * (phi_2 - anisotropy_axis)) ** 2) 

    energy_anis = Ku * np.sin(Phi_rad[:] - anisotropy_axis) ** 2
    energy_total = energy_H + energy_J_AF + energy_anis

    # Calculate the derivative of energy with respect to Phi_values (approximated using finite differences)
    dPhi = Phi_values[1] - Phi_values[0]
    dE_dPhi = np.gradient(energy_total, dPhi, edge_order=2)
    d2E_dPhi2 = np.gradient(dE_dPhi, dPhi, edge_order=2)

    return energy_total, dE_dPhi, d2E_dPhi2  # Include the second derivative in the return statement

### Example usage ###

# Parameters
H_values = np.linspace(1.5, -1.5, 9)
Phi_values = np.linspace(0, 360,1000)
J_AF = 0
M = -1 
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

#Inflection Point and Minimum Calculations (Decreasing Field)
for i, H in enumerate(H_values):
    energy_values, dE_dPhi, d2E_dPhi2 = calculate_energy_sweep(H, Phi_values + 180, J_AF, M, Ku, K1, anisotropy_axis )
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
                        easy_axis_min_angles_list.append([easy_axis_min_angle , inflection_point])
                        inflection_included_list[i] = True
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

'---'



# Accumulate the easy axis minimum angles and corresponding field strengths
min_angle_degrees_list = []
field_strengths_list = []

# Pair the found easy axis minimums with their corrosponding field strengths. 
for i, H in enumerate(H_values[:-1]): #Excludes last point 
    easy_axis_min_angles = easy_axis_min_angles_list[i]
    if easy_axis_min_angles is not None:
        # Convert the angle from radians to degrees
        easy_axis_min_angle_deg = easy_axis_min_angles
        # Append the easy axis minimum angles and the corresponding field strength to the lists
        if isinstance(easy_axis_min_angles, list):
            # Handle the case of multiple angles for the same field strength
            for angle in easy_axis_min_angles:
                min_angle_degrees_list.append(round(np.cos(np.radians(angle)), 2))
                field_strengths_list.append(H)
        else:
            min_angle_degrees_list.append( round(np.cos(np.radians(easy_axis_min_angle_deg)), 2))
            field_strengths_list.append(H)


# Convert the field_strengths_list to a NumPy array
field_strengths_array = np.array(field_strengths_list, dtype=float)
plt.scatter(field_strengths_array, min_angle_degrees_list, color='green', marker='o')
plt.plot(field_strengths_array, min_angle_degrees_list)



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

H_values_Inc = H_values[1:]

for i, H in enumerate(H_values_Inc[::-1]):
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


min_angle_degrees_list = []
field_strengths_list = []
for i, H in enumerate(H_values_Inc[::-1]):
    easy_axis_min_angles = easy_axis_min_angles_list[i]
    if easy_axis_min_angles is not None:
        # Convert the angle from radians to degrees
        easy_axis_min_angle_deg = easy_axis_min_angles
        # Append the easy axis minimum angles and the corresponding field strength to the lists
        if isinstance(easy_axis_min_angles, list):
            # Handle the case of multiple angles for the same field strength
            for angle in easy_axis_min_angles:
                min_angle_degrees_list.append(round(np.cos(np.radians(angle)), 2))
                field_strengths_list.append(H)
        else:
            min_angle_degrees_list.append( round(- 1 * np.cos(np.radians(easy_axis_min_angle_deg)), 2))
            field_strengths_list.append(H)




# Convert the field_strengths_list to a NumPy array
field_strengths_array = np.array(field_strengths_list, dtype=float)
plt.scatter(field_strengths_array, min_angle_degrees_list, color='red', marker='x')
plt.plot(field_strengths_array, min_angle_degrees_list)



#PLOT

plt.xlabel('Applied Field Strength H (T)')
plt.ylabel('M/M_s')
plt.title('M/M_s vs Applied Field Strength')
plt.grid(True)
plt.show()

        

'''
Debugging
'''
###Plot energy_values for each H
##plt.figure(figsize=(10, 5))
##for i, H in enumerate(H_values):
##    energy_values = energy_values_list[i]
##    label_text = f'H={H}'
##    plt.plot(Phi_values[:],energy_values, label=label_text)
### Add the legend to the plot
##plt.legend()
##x_ticks = np.arange(0, 361, 45) 
##plt.xticks(x_ticks)
##
##plt.show()

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

'''
The pupose of this code is to simulate a 2 Dimensional 1 spin system and it's response to an applied field
'''

def calculate_energy_sweep(H, Phi_values, J_AF, M, Ku, K1, anisotropy_axis):
    # Convert angles to radians
    Phi_rad = np.radians(Phi_values)
    
    energy_H = -1 * H * M * (np.cos(Phi_rad[:]))
    #energy_J_AF = -J_AF * np.cos(phi_1 - phi_2)
    
    ''' Uniaxial Term ''' 
    #energy_anis = (1/8) * K1 * (np.sin(2 * (phi_1 - anisotropy_axis)) ** 2 + np.sin(2 * (phi_2 - anisotropy_axis)) ** 2)
    energy_anis = Ku * np.sin(Phi_rad[:] - anisotropy_axis) ** 2
    energy_total = energy_H + 0 + energy_anis

    # Calculate the derivative of energy with respect to Phi_values (approximated using finite differences)
    dPhi = Phi_values[1] - Phi_values[0]
    dE_dPhi = np.gradient(energy_total, dPhi, edge_order=2)
    d2E_dPhi2 = np.gradient(dE_dPhi, dPhi, edge_order=2)

    return energy_total, dE_dPhi, d2E_dPhi2  # Include the second derivative in the return statement

# Example usage 
H_values = np.linspace(1.5, -1.5, 9)
Phi_values = np.linspace(0, 360,1000)
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

#Inflection Point and Minimum Calculations (Decreasing Field)
for i, H in enumerate(H_values):
    energy_values, dE_dPhi, d2E_dPhi2 = calculate_energy_sweep(H, Phi_values, J_AF, M, Ku, K1, anisotropy_axis )
    energy_values_list.append(energy_values)

    dE_dPhi_list.append(dE_dPhi)
    d2E_dPhi2_list.append(d2E_dPhi2)

    # Find the points of inflection using the sign change in the second derivative
    inflection_indices = np.where(np.diff(np.sign(d2E_dPhi2)))[0]
    inflection_points = Phi_values[:-1][inflection_indices] 


    energy_min_indices, _ = find_peaks(-energy_values)  # Negative energy values to find minima
    if len(energy_min_indices) > 0:
        min_energy_index = energy_min_indices[0]
        easy_axis_min_angle = Phi_values[:][min_energy_index]
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
        #print(f'{H}:')
        #print(easy_axis_min_angle_deg)
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


#END OF DECREASING FIELD PLOTTING



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
    energy_values, dE_dPhi, d2E_dPhi2 = calculate_energy_sweep(H, Phi_values + 180, J_AF, M, Ku, K1, anisotropy_axis)
    energy_values_list.append(energy_values)

    dE_dPhi_list.append(dE_dPhi)
    d2E_dPhi2_list.append(d2E_dPhi2)

    # Find the points of inflection using the sign change in the second derivative
    inflection_indices = np.where(np.diff(np.sign(d2E_dPhi2)))[0]
    inflection_points = Phi_values[:-1][inflection_indices]

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

    inflection_points_list.append(inflection_points)



min_angle_degrees_list = []
field_strengths_list = []
for i, H in enumerate(H_values_Inc[::-1]):
    easy_axis_min_angles = easy_axis_min_angles_list[i]
    if easy_axis_min_angles is not None:
        # Convert the angle from radians to degrees
        easy_axis_min_angle_deg = easy_axis_min_angles
        #print(f'{H}:')
        #print(easy_axis_min_angle_deg)
        # Append the easy axis minimum angles and the corresponding field strength to the lists
        if isinstance(easy_axis_min_angles, list):
            # Handle the case of multiple angles for the same field strength
            for angle in easy_axis_min_angles:
                min_angle_degrees_list.append(round(np.cos(np.radians(angle)), 2))
                field_strengths_list.append(H)
        else:
            min_angle_degrees_list.append( round(- 1 * np.cos(np.radians(easy_axis_min_angle_deg)), 2))
            field_strengths_list.append(H)


# Hysteresis Chart

### Convert the field_strengths_list to a NumPy array
field_strengths_array = np.array(field_strengths_list, dtype=float)
plt.scatter(field_strengths_array, min_angle_degrees_list, color='red', marker='x')

plt.xlabel('Applied Field Strength H (T)')
plt.ylabel('M/M_s')
plt.title('M/M_s vs Applied Field Strength')
plt.grid(True)
plt.show()

        

'''
Debugging (Energy Plots)
'''
###Plot energy_values for each H
##plt.figure(figsize=(10, 5))
##for i, H in enumerate(H_values):
##    energy_values = energy_values_list[i]
##    print(energy_values_list)
##    label_text = f'H={H}'
##    plt.plot(Phi_values[:],energy_values, label=label_text)
### Add the legend to the plot
##plt.legend()
##x_ticks = np.arange(0, 361, 45) 
##plt.xticks(x_ticks)
##
##plt.show()

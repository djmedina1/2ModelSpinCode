import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

def calculate_energy(H, Phi_values, J_AF, M, Ku, K1, anisotropy_axis):
    # Your energy calculation code here
    # ...

def find_easy_axis_min_angles(H_values, Phi_values, J_AF, M, Ku, K1, anisotropy_axis):
    easy_axis_min_angles_list = []
    # Rest of your calculation and inflection point code here
    # ...

    return easy_axis_min_angles_list

def main():
    # Constants and variables
    # ...

    # Calculate easy axis minimum angles for decreasing field
    easy_axis_min_angles_decreasing = find_easy_axis_min_angles(
        H_values, Phi_values + 180, J_AF, M, Ku, K1, anisotropy_axis
    )

    # Calculate easy axis minimum angles for increasing field (reverse H_values)
    easy_axis_min_angles_increasing = find_easy_axis_min_angles(
        H_values[::-1], Phi_values, J_AF, M, Ku, K1, anisotropy_axis
    )

    # Extract data for plotting
    min_angle_degrees_list = []
    field_strengths_list = []

    for H, angles_list in zip(H_values, easy_axis_min_angles_decreasing):
        if angles_list is not None:
            for angle in angles_list:
                min_angle_degrees_list.append(round(np.cos(np.radians(angle)), 2))
                field_strengths_list.append(H)

    for H, angles_list in zip(H_values[1:], easy_axis_min_angles_increasing):
        if angles_list is not None:
            for angle in angles_list:
                min_angle_degrees_list.append(round(-1 * np.cos(np.radians(angle)), 2))
                field_strengths_list.append(H)

    # Convert the field_strengths_list to a NumPy array
    field_strengths_array = np.array(field_strengths_list, dtype=float)

    # Plotting
    plt.scatter(field_strengths_array, min_angle_degrees_list, color='green', marker='o')
    # ...

    # Show the plot
    plt.xlabel('Applied Field Strength H (T)')
    plt.ylabel('M/M_s')
    plt.title('M/M_s vs Applied Field Strength')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()

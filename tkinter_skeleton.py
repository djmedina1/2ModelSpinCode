import tkinter as tk
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# The rest of your simulation functions here...

def plot_graph():
    # Call your simulation functions and perform calculations
    # ...

    # Plotting
    plt.scatter(field_strengths_array, min_angle_degrees_list, color='green', marker='o')
    # ...

    # Show the plot
    plt.xlabel('Applied Field Strength H (T)')
    plt.ylabel('M/M_s')
    plt.title('M/M_s vs Applied Field Strength')
    plt.grid(True)
    plt.show()

def run_simulation():
    # Call your simulation functions and perform calculations
    # ...

    # Generate the plot
    plot_graph()

def on_click_run():
    run_simulation()

# Create the main window
root = tk.Tk()
root.title("2D 1-Spin System Simulation")

# Add GUI elements
# ...

# Create a button to run the simulation
run_button = tk.Button(root, text="Run Simulation", command=on_click_run)
run_button.pack()

# Start the main event loop
root.mainloop()

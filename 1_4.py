""" fd3d_1_4.py: 1D FDTD
Simulation of a sinusoidal wave hitting a dielectric medium
"""
import numpy as np
from math import pi, sin
from matplotlib import pyplot as plt
ke = 200
ex = np.zeros(ke)
hy = np.zeros(ke)
ddx = 0.01 # Cell size
dt = ddx / 6e8 # Time step size
freq_in = 700e6
boundary_low = [0, 0]
boundary_high = [0, 0]
# Create Dielectric Profile
cb = np.ones(ke)
cb = 0.5 * cb
cb_start = 100
epsilon = 4
cb[cb_start:] = 0.5 / epsilon
nsteps = 425
# Dictionary to keep track of desired points for plotting
plotting_points = [
    {'num_steps': 150, 'data_to_plot': None, 'label': ''},
    {'num_steps': 425, 'data_to_plot': None, 'label': 'FDTD cells'}
]
# Main FDTD Loop
for time_step in range(1, nsteps + 1):
    # Calculate the Ex field
    for k in range(1, ke):
        ex[k] = ex[k] + cb[k] * (hy[k - 1] - hy[k])
    # Put a sinusoidal at the low end
    pulse = sin(2 * pi * freq_in * dt * time_step)
    ex[5] = pulse + ex[5]
    # Absorbing Boundary Conditions
    ex[0] = boundary_low.pop(0)
    boundary_low.append(ex[1])
    ex[ke - 1] = boundary_high.pop(0)
    boundary_high.append(ex[ke - 2])
    # Calculate the Hy field
    for k in range(ke - 1):
        hy[k] = hy[k] + 0.5 * (ex[k] - ex[k + 1])
    # Save data at certain points for later plotting
    for plotting_point in plotting_points:
        if time_step == plotting_point['num_steps']:
            plotting_point['data_to_plot'] = np.copy(ex)
# Plot the outputs in Fig. 1.5
plt.rcParams['font.size'] = 12
fig = plt.figure(figsize=(8, 3.5))
def plot_e_field(data, timestep, epsilon, cb, label):
    """Plot of E field at a single time step"""
    plt.plot(data, color='k', linewidth=1)
    plt.ylabel('E$_x$', fontsize='14')
    plt.xticks(np.arange(0, 199, step=20))
    plt.xlim(0, 199)
    plt.yticks(np.arange(-1, 1.2, step=1))
    plt.ylim(-1.2, 1.2)
    plt.text(50, 0.5, 'T = {}'.format(timestep),
    horizontalalignment='center')
    plt.plot((0.5 / cb - 1) / 3, 'k--', linewidth=0.75)
    # The math on cb above is just for scaling
    plt.text(170, 0.5, 'Eps = {}'.format(epsilon),
    horizontalalignment='center')
    plt.xlabel('{}'.format(label))
# Plot the E field at each of the time steps saved earlier
for subplot_num, plotting_point in enumerate(plotting_points):
    ax = fig.add_subplot(2, 1, subplot_num + 1)
    plot_e_field(plotting_point['data_to_plot'],
        plotting_point['num_steps'], epsilon, cb,
        plotting_point['label'])
plt.subplots_adjust(bottom=0.2, hspace=0.45)
plt.show()

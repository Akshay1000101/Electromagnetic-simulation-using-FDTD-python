""" fd3d_1_5.py: 1D FDTD
Simulation of a sinusoid wave hitting a lossy dielectric
"""


import numpy as np
from math import pi, sin
from matplotlib import pyplot as plt


ke = 200
ex = np.zeros(ke)
hy = np.zeros(ke)


22
ONE-DIMENSIONAL SIMULATION WITH THE FDTD METHOD
ddx = 0.01
# Cell size
dt = ddx / 6e8
# Time step size
freq_in = 700e6
boundary_low = [0, 0]
boundary_high = [0, 0]
# Create Dielectric Profile
epsz = 8.854e-12
epsilon = 4
sigma = 0.04
ca = np.ones(ke)
cb = np.ones(ke) * 0.5
cb_start = 100
eaf = dt * sigma / (2 * epsz * epsilon)
ca[cb_start:] = (1 - eaf ) / (1 + eaf )
cb[cb_start:] = 0.5 / (epsilon * (1 + eaf ))
nsteps = 500
# Main FDTD Loop
for time_step in range(1, nsteps + 1):
# Calculate the Ex field
for k in range(1, ke):
ex[k] = ca[k] * ex[k] + cb[k] * (hy[k - 1] - hy[k])
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
# Plot the outputs in Fig. 1.6
plt.rcParams['font.size'] = 12
plt.figure(figsize=(8, 2.25))
ONE-DIMENSIONAL SIMULATION WITH THE FDTD METHOD
23
plt.plot(ex, color='k', linewidth=1)
plt.ylabel('E$_x$', fontsize='14')
plt.xticks(np.arange(0, 199, step=20))
plt.xlim(0, 199)
plt.yticks(np.arange(-1, 1.2, step=1))
plt.ylim(-1.2, 1.2)
plt.text(50, 0.5, 'T = {}'.format(time_step),
horizontalalignment='center')
plt.plot((0.5 / cb - 1) / 3, 'k--',
linewidth=0.75)
# The math on cb is just for scaling
plt.text(170, 0.5, 'Eps = {}'.format(epsilon),
horizontalalignment='center')
plt.text(170, -0.5, 'Cond = {}'.format(sigma),
horizontalalignment='center')
plt.xlabel('FDTD cells')
plt.subplots_adjust(bottom=0.25, hspace=0.45)
plt.show()
24
ONE-DIMENSIONAL SIMULATION WITH THE FDTD METHOD

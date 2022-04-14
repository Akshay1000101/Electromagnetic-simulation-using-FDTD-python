""" fd2d_3_2.py: 2D FDTD
TM program with the PML added
"""
import numpy as np
from math import sin, pi
from matplotlib import pyplot as plt
import mpl_toolkits.mplot3d.axes3d

ie = 60
je = 60
ic = int(ie / 2 - 5)
jc = int(je / 2 - 5)


ez = np.zeros((ie, je))
dz = np.zeros((ie, je))
hx = np.zeros((ie, je))
hy = np.zeros((ie, je))
ihx = np.zeros((ie, je))
ihy = np.zeros((ie, je))

ddx = 0.01
# Cell size
dt = ddx / 6e8  # Time step size


# Create Dielectric Profile
epsz = 8.854e-12

# Pulse Parameters
t0 = 40
spread = 12

gaz = np.ones((ie, je))

# Calculate the PML parameters
gi2 = np.ones(ie)
gi3 = np.ones(ie)
fi1 = np.zeros(ie)
fi2 = np.ones(ie)
fi3 = np.ones(ie)

gj2 = np.ones(ie)
gj3 = np.ones(ie)
fj1 = np.zeros(ie)
fj2 = np.ones(ie)
fj3 = np.ones(ie)


# Create the PML as described in Section 3.2
npml = 8
for n in range(npml):
    xnum = npml - n
    xd = npml
    xxn = xnum / xd
    xn = 0.33 * xxn ** 3
    gi2[n] = 1 / (1 + xn)
    gi2[ie - 1 - n] = 1 / (1 + xn)
    gi3[n] = (1 - xn) / (1 + xn)
    gi3[ie - 1 - n] = (1 - xn) / (1 + xn)
    gj2[n] = 1 / (1 + xn)
    gj2[je - 1 - n] = 1 / (1 + xn)
    gj3[n] = (1 - xn) / (1 + xn)
    gj3[je - 1 - n] = (1 - xn) / (1 + xn)
    xxn = (xnum - 0.5) / xd
    xn = 0.33 * xxn ** 3
    fi1[n] = xn
    fi1[ie - 2 - n] = xn
    fi2[n] = 1 / (1 + xn)
    fi2[ie - 2 - n] = 1 / (1 + xn)
    fi3[n] = (1 - xn) / (1 + xn)
    fi3[ie - 2 - n] = (1 - xn) / (1 + xn)
    fj1[n] = xn
    fj1[je - 2 - n] = xn
    fj2[n] = 1 / (1 + xn)
    fj2[je - 2 - n] = 1 / (1 + xn)
    fj3[n] = (1 - xn) / (1 + xn)
    fj3[je - 2 - n] = (1 - xn) / (1 + xn)

nsteps = 100

# Dictionary to keep track of desired points for plotting
plotting_points = [
    {'num_steps': 40, 'data_to_plot': None},
    {'num_steps': nsteps, 'data_to_plot': None},
]
# Main FDTD Loop
for time_step in range(1, nsteps + 1):
   
    # Calculate Dz
    for j in range(1, je):
         for i in range(1, ie):
            dz[i, j] = gi3[i] * gj3[j] * dz[i, j] + \
                        gi2[i] * gj2[j] * 0.5 * \
                        (hy[i, j] - hy[i - 1, j] -
                        hx[i, j] + hx[i, j - 1])
    # Put a Gaussian pulse in the middle
    pulse = sin(2 * pi * 1500 * 1e6 * dt * time_step)
    dz[ic, jc] = pulse
    
    ez = gaz * dz
    
    # Calculate the Ez field from Dz
    # Calculate the Hx field
    for j in range(je - 1):
        for i in range(ie - 1):
            curl_e = ez[i, j] - ez[i, j + 1]
            ihx[i, j] = ihx[i, j] + curl_e
            hx[i, j] = fj3[j] * hx[i, j] + fj2[j] * \
                        (0.5 * curl_e + fi1[i] * ihx[i, j])

# Calculate the Hy field
for j in range(0, je - 1):
    for i in range(0, ie - 1):
        curl_e = ez[i, j] - ez[i + 1, j]
        ihy[i, j] = ihy[i, j] + curl_e
        hy[i, j] = fi3[i] * hy[i, j] - fi2[i] * \
                (0.5 * curl_e + fj1[j] * ihy[i, j])
    
    # Save data at certain points for later plotting
    for plotting_point in plotting_points:
        if time_step == plotting_point['num_steps']:
            plotting_point['data_to_plot'] = np.copy(ez)
    
# Plot Fig. 3.4
plt.rcParams['font.size'] = 12
plt.rcParams['grid.color'] = 'gray'
plt.rcParams['grid.linestyle'] = 'dotted'
fig = plt.figure(figsize=(8, 8))

X, Y = np.meshgrid(range(je), range(ie))
def plot_e_field(ax, data, timestep):
    """3d Plot of E field at a single time step"""
    ax.set_zlim(-0.5, 0.5)
    ax.view_init(elev=30., azim=-135)
    ax.plot_surface(X, Y, data, rstride=1, cstride=1, color='white',edgecolor='black', linewidth=.25)
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel(r' $E_{Z}$', rotation=90, labelpad=10,fontsize=14)
    ax.set_xlabel('cm')
    ax.set_ylabel('cm')
    ax.set_xticks(np.arange(0, 60, step=20))
    ax.set_yticks(np.arange(0, 60, step=20))
    ax.set_zticks([-0.5, 0, 0.5])
    ax.text2D(0.6, 0.7, "T = {}".format(timestep),transform=ax.transAxes)
    ax.xaxis.pane.fill = ax.yaxis.pane.fill = \
        ax.zaxis.pane.fill = False
    plt.gca().patch.set_facecolor('white')
    ax.dist = 11


def plot_e_field_contour(ax, data):
    """Contour Plot of E field at a single time step"""
    CP = plt.contour(X, Y, data, colors='black',linestyles='solid')

    CP.collections[4].remove()
    # above removes extraneous outer contour display
    ax.set_xticks(np.arange(0, 60, step=20))
    ax.set_yticks(np.arange(0, 60, step=20))
    plt.xlabel('cm')
    plt.ylabel('cm')
   
    # Plot the E field at each of the time steps saved earlier
    for subplot_num, plotting_point in enumerate(plotting_points):
        ax = fig.add_subplot(2, 2, subplot_num * 2 + 1,projection='3d')
    plot_e_field(ax, plotting_point['data_to_plot'],plotting_point['num_steps'])
    ax = fig.add_subplot(2, 2, subplot_num * 2 + 2)
    plot_e_field_contour(ax, plotting_point['data_to_plot'])
plt.tight_layout()
plt.subplots_adjust(left=0.05)
plt.show()

import numpy as np
from math import exp
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D, get_test_data
ie = 60
je = 60
ic = int(ie / 2)
jc = int(je / 2)
ia = 7
ib = ie - ia - 1
ja = 7
jb = je - ja - 1
ez = np.zeros((ie, je))
dz = np.zeros((ie, je))
hx = np.zeros((ie, je))
hy = np.zeros((ie, je))
ihx = np.zeros((ie, je))
ihy = np.zeros((ie, je))
ez_inc = np.zeros(je)
hx_inc = np.zeros(je)
ddx = 0.01 # Cell size
dt = ddx / 6e8 # Time step size
# Absorbing Boundary Conditions
boundary_low = [0, 0]
boundary_high = [0, 0]
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
# Create Dielectric Profile
epsz = 8.854e-12
# Pulse Parameters
t0 = 20
spread = 8
gaz = np.ones((ie, je))
nsteps = 115
# Dictionary to keep track of desired points for plotting
plotting_points = [
    {'label': 'a', 'num_steps': 30, 'data_to_plot': None},
    {'label': 'b', 'num_steps': 60, 'data_to_plot': None},
    {'label': 'c', 'num_steps': 90, 'data_to_plot': None},
    {'label': 'd', 'num_steps': 115, 'data_to_plot': None},
]
# Main FDTD Loop
for time_step in range(1, nsteps + 1):
    for j in range(1, je):
        ez_inc[j] = ez_inc[j] + 0.5 * (hx_inc[j - 1] - hx_inc[j])
    # Absorbing Boundary Conditions
    ez_inc[0] = boundary_low.pop(0)
    boundary_low.append(ez_inc[1])
    ez_inc[je - 1] = boundary_high.pop(0)
    boundary_high.append(ez_inc[je - 2])
    # Calculate the Dz field
    for j in range(1, je):
        for i in range(1, ie):
            dz[i, j] = gi3[i] * gj3[j] * dz[i, j] + \
                    gi2[i] * gj2[j] * 0.5 * \
                    (hy[i, j] - hy[i - 1, j] -
                    hx[i, j] + hx[i, j - 1])
    # Source
    pulse = exp(-0.5 * ((t0 - time_step) / spread) ** 2)
    ez_inc[3] = pulse
    # Incident Dz values
    for i in range(ia, ib):
        dz[i, ja] = dz[i, ja] + 0.5 * hx_inc[ja - 1]
        dz[i, jb] = dz[i, jb] - 0.5 * hx_inc[jb - 1]
    ez = gaz * dz # Calculate the Ez field from Dz
    for j in range(0, je - 1):
        hx_inc[j] = hx_inc[j] + 0.5 * (ez_inc[j] - ez_inc[j + 1])
    for j in range(0, je - 1):
        for i in range(0, ie - 1):
            curl_e = ez[i, j] - ez[i, j + 1]
            ihx[i, j] = ihx[i, j] + curl_e
            hx[i, j] = fj3[j] * hx[i, j] + fj2[j] * \
            (0.5 * curl_e + fi1[i] * ihx[i, j])
    # Incident Hx values
    for i in range(ia, ib):
        hx[i, ja - 1] = hx[i, ja - 1] + 0.5 * ez_inc[ja]
        hx[i, jb] = hx[i, jb] - 0.5 * ez_inc[jb]
    # Calculate the Hy field
    for j in range(0, je - 1):
        for i in range(0, ie - 1):
            curl_e = ez[i, j] - ez[i + 1, j]
            ihy[i, j] = ihy[i, j] + curl_e
            hy[i, j] = fi3[i] * hy[i, j] - fi2[i] * \
                    (0.5 * curl_e + fj1[j] * ihy[i, j])
    # Incident Hy values
    for j in range(ja, jb):
        hy[ia - 1, j] = hy[ia - 1, j] - 0.5 * ez_inc[j]
        hy[ib - 1, j] = hy[ib - 1, j] + 0.5 * ez_inc[j]
    # Save data at certain points for later plotting
    for plotting_point in plotting_points:
        if time_step == plotting_point['num_steps']:
            plotting_point['data_to_plot'] = np.copy(ez)
# Plot Fig. 3.7
plt.rcParams['font.size'] = 12
plt.rcParams['grid.color'] = 'gray'
plt.rcParams['grid.linestyle'] = 'dotted'
fig = plt.figure(figsize=(8, 8))
X, Y = np.meshgrid(range(je), range(ie))
def plot_e_field(ax, data, timestep, label):
    """3d Plot of E field at a single timestep"""
    ax.set_zlim(-0.5, 1)
    ax.view_init(elev=15., azim=25)
    ax.plot_surface(Y, X, data, rstride=1, cstride=1,
            color='white',
            edgecolor='black', linewidth=.25)
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel(r' $E_{Z}$', rotation=90,
                labelpad=10, fontsize=14)
    ax.set_zticks([-0.5, 0, 0.5, 1])
    ax.set_xlabel('cm')
    ax.set_ylabel('cm')
    ax.set_xticks(np.arange(0, 61, step=20))
    ax.set_yticks(np.arange(0, 61, step=20))
    ax.text2D(0.25, 0.3, "T = {}".format(timestep),
                transform=ax.transAxes)
    ax.xaxis.pane.fill = ax.yaxis.pane.fill = \
            ax.zaxis.pane.fill = False
    plt.gca().patch.set_facecolor('white')
    ax.text2D(-0.05, 0.8, "({})".format(label),
            transform=ax.transAxes)
    ax.dist = 11
# Plot the E field at each of the four time steps saved earlier
for subplot_num, plotting_point in enumerate(plotting_points):
    ax = fig.add_subplot(2, 2, subplot_num + 1, projection='3d')
    plot_e_field(ax, plotting_point['data_to_plot'],
                plotting_point['num_steps'], plotting_point
                ['label'])
fig.tight_layout()
plt.show()

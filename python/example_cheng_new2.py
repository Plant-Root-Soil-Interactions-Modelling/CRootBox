# from pyevtk.hl import gridToVTK
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors
import py_rootbox as rb
from rb_tools import *


def v2a(vd):  # rb.std_vector_double_ to numpy array
    l = np.zeros(len(vd))
    for i in range(0, len(vd)):
        l[i] = vd[i]
    return l


#
# Root system
#
rootsystem = rb.RootSystem()
name = "Anagallis_straight_simple"  # name = "Zea"
rootsystem.openFile(name)

for i in range(0, 10):
    p = rootsystem.getRootTypeParameter(i + 1)
    p.gf = 2  # linear growth function

rootsystem.initialize()
simtime = 2
rootsystem.simulate(simtime, True);
rootsystem.write(name + ".vtp")

#
# Grid parameter
#
nodes = vv2a(rootsystem.getNodes())
boxmin = nodes.min(axis = 0); boxmax = nodes.max(axis = 0);  # cm
width = abs(max(boxmax[0], boxmax[1]) - min(boxmin[0], boxmin[1])) + 6;  # cm
depth = abs(boxmin[2]) + 3
xres = 0.3;
yres = 0.3;
zres = 0.3;
nx = int(width / xres);
ny = int(width / yres);
nz = int(depth / zres);
print("Width", width, "Depth", depth, " at a Resolution", nx, ny, nz)

#
# Model parameter
#
model = rb.ExudationModel(width, width, depth, nx, ny, nz, rootsystem)
model.Q = 4  # µg/d/tip
model.Dl = 2.43e-6 * 3600 * 24  # cm2/d
model.theta = 0.3
model.R = 16.7  # 16.7  # -
model.k = 2.60e-6 * 3600 * 24  # d-1
model.l = 4  # cm (for line source only)

#
# Numerical parameter
#
model.type = rb.IntegrationType.mls;  # mps, mps_straight, mls
model.n0 = 10  # integration points per cm
model.thresh13 = 1.e-15;  # threshold to neglect diffusing g (eqn 13)
model.calc13 = True;  # turns Eqn 13  on (True) and off (False)
model.observationRadius = 5;  # limits computational domain around roots [cm]

t = time.time()
C = model.calculate(0, simtime)
elapsed = time.time() - t
print("Computation took", elapsed, "s")

#
# post processing...
#
C = np.reshape(v2a(C), (nx, ny, nz))  # hope that works, it does not :-(, or does it?
X = np.linspace(-width / 2, width / 2, nx)
Y = np.linspace(-width / 2, width / 2, ny)
Z = np.linspace(-depth, 0, nz)
X_, Y_, Z_ = np.meshgrid(X, Y, Z, indexing = "ij")

num_th = (C > 0).sum()  # number of points for which concentration is larger than threshold
print("volume of concentration above threshold: ", num_th * width / nx * width / ny * depth / nz)
print("this is", num_th / (nx * ny * nz) * 100, "% of the overall volume")

# gridToVTK("./Exudates", X, Y, Z, pointData = {"Exudates":C})

fig1 = plt.figure()
ax = plt.axes()
C_ = C[:, int(ny / 2), :]
levels = np.logspace(np.log10(np.max(C_)) - 5, np.log10(np.max(C_)), 100)
cs = ax.contourf(X_[:, int(ny / 2), :], Z_[:, int(ny / 2), :], C_, levels = levels, locator = ticker.LogLocator(), cmap = 'jet')
ax.set_xlabel('x')
ax.set_ylabel('z')
plt.axis('equal')
cbar = fig1.colorbar(cs)

fig2 = plt.figure()
nodes = rootsystem.getNodes()
node_tip = nodes[-10]
idy = (np.abs(Y - node_tip.y)).argmin()  # y-index of soil element closest to the point on the root axis
idz = (np.abs(Z - node_tip.z)).argmin()  # z-index of soil element closest to the point on the root axis
C_ = C[:, idy, idz]  # 1d array
plt.plot(X, C_, 'k-')
plt.xlabel('x')
plt.ylabel('c')

fig3 = plt.figure()
nodes = rootsystem.getNodes()
node_tip = nodes[-1]
idy = (np.abs(Y - node_tip.y)).argmin()  # y-index of soil element closest to the point on the root axis
idx = (np.abs(X - node_tip.x)).argmin()  # x-index of soil element closest to the point on the root axis
C_ = C[idx, idy, :]  # 1d array
plt.plot(Z, C_, 'k-')
plt.xlabel('z')
plt.ylabel('c')

plt.show()

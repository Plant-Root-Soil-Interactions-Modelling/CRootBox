import math

import numpy as np
import scipy
from scipy import integrate

import matplotlib.pyplot as plt

#
# Model parameter
#
Q = 4  # µg/d/tip
Dl = 2.43e-6 * 3600 * 24  # cm2/d
theta = 0.3
R = 16.7  # 16.7  # -
k = 2.60e-6 * 3600 * 24  # d-1
l = 4  # cm (for line source only)

r = 2  # initial growth speed (linear growth)
simtime = 10  # days
depth = 26  # cm


def to32(x):
    return np.sqrt(x * x * x)


def fun11a(x, t):
    c = -R / (4 * Dl * t)
    d = 8 * theta * to32(math.pi * Dl * t)
    xtip = -r * t - 3
    z = x - xtip
    return (Q * math.sqrt(R)) / d * np.exp(c * z * z - k / R * t)


def fun11b(x, t, l):
    c = -R / (4 * Dl * t)
    d = 8 * theta * to32(math.pi * Dl * t);
    xtip = -r * t - 3
    z = x - xtip - l;
    return (Q * math.sqrt(R)) / d * math.exp(c * z * z - k / R * t)


N = 200
z_ = np.linspace(0, -depth, N)
c = np.zeros((N,))
tend = simtime

for i, z in enumerate(z_):
    f = lambda t: fun11a(z, t)
    res = integrate.quad(f, 0, tend)
    c[i] = res[-1]

plt.plot(z_, c)
plt.show()

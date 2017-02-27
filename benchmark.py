import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#
# Exact lengths of a root system with 0 and 1 order without standard deviation
#

def maxRootLength(la,lb,ln,nob): # maximal length the root will reach
    return la+lb+(nob-1)*ln


def rootLength(t,r,k): # root length at a certain age
    return k*(1-np.exp(-r*t/k))


def rootAge(l,r,k): # root age at a certain length
    return -np.log(1-l/k)*k/r


def rootLateralLength(t,et,r,k): # length of first order laterals (without second order laterals)
    i = 0
    l = 0
    while et[i] <t:
        age = t-et[i]
        l += rootLength(age,r,k)
        i += 1
    return l

#
# Root parameter
#

# 0 order root
la = 10.
nob = 20.
lb = 1.
ln = 89./19.
r = 1.
k = maxRootLength(la,lb,ln,nob)

# 1st order root
k1 = 25.
r1 = 2.

#
# Plant parameter (neglecting shoot borne)
#
maxB = 100
firstB = 10
delayB = 3

times = np.array([7.,15.,30.,60.])

#
# Benchmark 1: single root, no lateral
#
print("Benchmark 1")
l = rootLength(times,r,k)
print("times   ", times)
print("lenghts ", l, "\n")

#
# Benchmark 2: single root
#
print("Benchmark 2")
i=0
et = np.zeros(nob)
while i<nob:
    et[i] = rootAge(la+lb+ln*i,r,k+0.01)
    i += 1
# print("lateral emergence times", et)

l = rootLength(times,r,k) # zero order lengths
l1 = np.zeros(times.size)
j=0
for t in times :
    l1[j] = rootLateralLength(t,et,r1,k1)
    j=j+1

print("times               ", times)
print("zero order length   ", l)
print("first order length ", l1)
print("total length        ", l+l1, "\n")

#
# Benchmark 3: plant, no laterals
#
print("Benchmark 3")
etB = np.array(range(maxB))*delayB + np.ones(maxB)*firstB # basal root emergence times

bl = np.zeros(times.size)
j = 0 # time counter
for t in times:
    i = 0 # basal root counter
    while t-etB[i]>0:
        bl[j] += rootLength(t-etB[i],r,k)
        i += 1
    j += 1

print("times               ", times)
print("tap root lenght     ", l)
print("summed basal length ", bl)
print("total length        ", l+bl, "\n")

#
# Benchmark 4: plant, with laterals
#
print("Benchmark 4")
# etB as berfore
# et a before

bl = np.zeros(times.size)
j = 0 # time counter
for t in times:
    i = 0 # basal root counter
    while t-etB[i]>0:
        bl[j] += ( rootLateralLength(t-etB[i],et,r1,k1) + rootLength(t-etB[i],r,k) )
        i += 1
    j += 1

print("times               ", times)
print("tap root lenght     ", l+l1)
print("summed basal length ", bl)
print("total length        ", l+l1+bl, "\n")

#
# Matlab like plotting
#
t_ = np.linspace(0,100,100)
l_ = rootLength(t_,r,k)

l1_ = np.zeros(t_.size)
i = 0
for t in t_:
    l1_[i] = rootLateralLength(t, et, r1, k1)
    i += 1

plt.plot(t_,l_)
plt.plot(t_,l1_,'g')
plt.plot(times,l,'ro')
h0 = mpatches.Patch(color='blue', label='zero order')
h1 = mpatches.Patch(color='green', label='first order')

plt.xlabel("Age (days)")
plt.ylabel("Length (cm)")
plt.legend([h0, h1],['zero order','first order'])
plt.show()


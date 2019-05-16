import py_rootbox as rb
from rb_tools import *
import math

rs = rb.RootSystem()

# Root type parameter
p0 = rb.RootTypeParameter() # with default values, 
p1 = rb.RootTypeParameter() # all standard deviations are 0

p0.name = "taproot"
p0.type = 1
p0.lb = 1
p0.la = 10
p0.nob = 20
p0.ln = 89./19.
p0.theta = 30./180.*math.pi
p0.r = 1
p0.dx = 0.5
p0.successor = a2i([2]) # add successors
p0.successorP = a2v([1])
p0.tropismT = rb.TropismType.gravi
p0.tropismN = 1.
p0.tropismS = 0.2

p1.name = "lateral"
p1.type = 2
p1.la = 25
p1.las = 10 # add standard deviation
p1.ln = 0
p1.r = 2
p1.dx = 0.1
p1.tropismS = 0.3

rs.setRootTypeParameter(p0)
rs.setRootTypeParameter(p1)

# Root system parameter (neglecting shoot borne)
maxB = 100
firstB = 10.
delayB = 3.
rsp = rb.RootSystemParameter()
rsp.set(-3., firstB, delayB, maxB, 0, 1.e9, 1.e9,  1.e9, 0., 0.)

rs.setRootSystemParameter(rsp)

rs.initialize()
rs.simulate(40, True)

rs.write("results/example_3c.vtp")

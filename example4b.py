import py_rootbox as rb
import math

rs = rb.RootSystem()
name = "Anagallis_femina_Leitner_2010" 
rs.openFile(name)

# box with a left and a right compartment for analysis
sideBox = rb.SDF_PlantBox(10,20,50)
left = rb.SDF_RotateTranslate(sideBox, rb.Vector3d(-4.99,0,0))
right = rb.SDF_RotateTranslate(sideBox, rb.Vector3d(4.99,0,0))
leftright = rb.SDF_Union(left,right)
rs.setGeometry(leftright)  

# left compartment has a minimum of 0, 1 elsewhere
maxS = 1. # maximal 
minS = 0. # minimal 
slope = 1. # [cm] linear gradient between min and ma
leftC = rb.SDF_Complement(left)
soilprop = rb.SoilPropertySDF(leftC, maxS, minS, slope) # for root elongation 
soilprop2 = rb.SoilPropertySDF(left, 1., 0.002, slope) # for branching

# class MySoil(rb.SoilPropertySDF):    # todo 
#     def getValue(p, r):
#         print("haha")
#         return 1.         
# mysoil = MySoil(leftC, maxS, minS, slope)

# Manually set scaling function and tropism parameters
sigma = [0.4, 1., 1., 1., 1. ] * 2
for i in range(0,10):  
    p = rs.getRootTypeParameter(i+1)
    p.dx = 0.25 # adjust resolution
    p.tropismS = sigma[i] 
    
#     # 1. Scale elongation
#     p.se = soilprop
    
#     # 2. Scale insertion angle
#     p.sa = soilprop
    
# 3. Scale branching probability
p = rs.getRootTypeParameter(2)
p.ln = p.ln/5
p.nob = p.nob*5
p = rs.getRootTypeParameter(3)
p.sbp = soilprop2

# simulation
simtime = 120.
dt = 1.
N = 120/dt

rs.initialize()
for i in range(0,round(N)):
    rs.simulate(dt,True)

# analyse
print()
print("Left compartment: ")
al = rb.SegmentAnalyser(rs)
al.crop(left)
ll = al.getSummed(rb.ScalarType.length)
print('Total root length', ll, 'cm')
lmct = al.getSummed(rb.ScalarType.time)/al.getSummed(rb.ScalarType.one)
print('Mean age', simtime-lmct, 'days')
lroots = al.getRoots()
lm_theta = 0
for r in lroots: 
    lm_theta += r.param.theta
lm_theta /= len(lroots)
print('Mean insertion angle is ', lm_theta/math.pi*180, 'degrees')
print()

print("Right compartment: ")
ar = rb.SegmentAnalyser(rs)
ar.crop(right)
lr = ar.getSummed(rb.ScalarType.length)
print('Total root length', lr, 'cm')
rmct = ar.getSummed(rb.ScalarType.time)/ar.getSummed(rb.ScalarType.one)
print('Mean age', simtime-rmct, 'days')
rroots = ar.getRoots()
rm_theta = 0
for r in lroots: 
    rm_theta += r.param.theta
rm_theta /= len(rroots)
print('Mean insertion angle is ', rm_theta/math.pi*180, 'degrees')
print()

# write results
rs.write("results/example_4b.py") # compartment geometry 
rs.write("results/example_4b.vtp") # root system

#
# Example 3
# 
# Creates N*N root systems at different locations, no confining geometry
# 
import py_rootbox as rb
import math
import random

name = "Maize_Pheno1_Leitner_et_al_2014"

N = 3
dist = 40 # cm
allRS = [ ]

#
# Creates N*N root systems
#
for i in range(0,N):
    for j in range(0,N):
         rs = rb.RootSystem()
         rs.openFile(name) 
         rs.getRootSystemParameter().seedPos = rb.Vector3d(dist*i,dist*j,-3) # set position of seed [cm]
         allRS.append(rs)

#
# Simulate
#
simtime = 120
for rs in allRS:
    rs.setSeed(random.random() )
    rs.initialize()
    rs.simulate(simtime)
    print(rs.getNumberOfNodes())

#
# Export results as single vtp files
#
c = 0
for rs in allRS:
      c += 1 # root system number
      vtpname = "results/"+name+"_"+str(c)+".vtp";
      rs.write(vtpname, rb.OutputType.polylines);
        
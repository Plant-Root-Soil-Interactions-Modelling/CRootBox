import py_rootbox as rb
from multiprocessing import Pool

name = "Zea_mays_4_Leitner_2014"
N = 3 # number of columns and rows
dist = 40 # distance between the root systems [cm]

# Creates and initializes N*N root systems
allRS = []
for i in range(0,N):
    for j in range(0,N):
         rs = rb.RootSystem()
         rs.openFile(name) 
         rs.getRootSystemParameter().seedPos = rb.Vector3d(dist*i,dist*j,-3) 
         allRS.append(rs)
         rs.initialize()

# Simulate parallel
simtime = 120
def simulate(i):
    allRS[i].simulate(simtime, True)
     
pool = Pool()
param_space = range(0,len(allRS))   
d = [1 for res in pool.imap(simulate,param_space)] 

# Export results as single vtp files (as polylines)
c = 0
ana = rb.SegmentAnalyser() # see example 3b
for rs in allRS:
      c += 1 # root system number
      vtpname = "results/example_2b_"+str(c)+".vtp"
      rs.write(vtpname)
      ana.addSegments(rs) # collect all
       
# Write all into single file (segments)
ana.write("results/example_2b_all.vtp") 
       
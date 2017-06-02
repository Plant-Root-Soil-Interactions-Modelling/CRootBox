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
         rs.getRootSystemParameter().seedPos = rb.Vector3d(dist*i,dist*j,-3) # set position of seed [cm]
         allRS.append(rs)
         rs.initialize()

# Simulate parallel
simtime = 120
def simulate(i):
    print("simualte ",i)
    allRS[i]. simulate(simtime)
     
pool = Pool()
param_space = range(0,len(allRS))    
for res in pool.imap(simulate,param_space):
    # do nothing
    c = 0
            
# Export results as single vtp files
c = 0
for rs in allRS:
      c += 1 # root system number
      vtpname = "results/"+name+"_"+str(c)+".vtp"
      rs.write(vtpname)
      
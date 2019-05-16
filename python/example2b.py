import py_rootbox as rb

name = "Zea_mays_4_Leitner_2014"
simtime = 120
N = 3 # number of columns and rows
dist = 40 # distance between the root systems [cm]

# Initializes N*N root systems
allRS = []
for i in range(0,N):
    for j in range(0,N):
        rs = rb.RootSystem()
        rs.openFile(name)      
        rs.getRootSystemParameter().seedPos = rb.Vector3d(dist*i,dist*j,-3.) # cm
        rs.initialize()
        allRS.append(rs)

# Simulate
for rs in allRS:    
    rs.simulate(simtime, True) 
    
# Export results as single vtp files (as polylines)
ana = rb.SegmentAnalyser() # see example 3b
for i,rs in enumerate(allRS):      
      vtpname = "results/example_2b_"+str(i)+".vtp"
      rs.write(vtpname)
      ana.addSegments(rs) # collect all
       
# Write all into single file (segments)
ana.write("results/example_2b_all.vtp") 
       
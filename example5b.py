import py_rootbox as rb

rs = rb.RootSystem()

# Open plant and root parameter from a file
name = "Cauliflower_Vansteenkiste_et_al_2014" 
rs.openFile(name) # as a second argument you could add the path (default = "modelparameter/"

# Initialize
rs.initialize() 

# simulation
simtime = 120.
dt = 1.
N = 120/dt

rs.initialize()
for i in range(0,round(N)):
    rs.simulate(dt,True)





# Export final result (as vtp)
rs.write("results/example_5b.vtp") # roots are exported as polylines

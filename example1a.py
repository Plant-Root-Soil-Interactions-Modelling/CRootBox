import py_rootbox as rb

rootsystem = rb.RootSystem()

# Open plant and root parameter from a file
name = "Anagallis_femina_Leitner_2010" 
rootsystem.openFile(name) # as a second argument you could add the path (default = "modelparameter/"

# Initialize
rootsystem.initialize() 

# Simulate
rootsystem.simulate(30)  # or 20, 40, 60 days

# Export final result (as vtp)
rootsystem.write("results/"+name+".vtp") # roots are exported as polylines

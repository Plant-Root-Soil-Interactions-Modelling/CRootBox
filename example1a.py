import py_rootbox as rb

rootsystem = rb.RootSystem()

# Open plant and root parameter from a file
name = "Anagallis_femina_Leitner_2010" 
rootsystem.openFile(name) 

# Initialize
rootsystem.initialize() 

# Simulate
rootsystem.simulate(30) 

# Export final result (as vtp)
rootsystem.write("results/example_1a.vtp") 
import py_rootbox as rb
from rb_tools import *

rootsystem = rb.RootSystem()

# Open plant and root parameter from a file
name = "Zea_mays_1_Leitner_2010"  # "Anagallis_femina_Leitner_2010"
rootsystem.openFile(name)

# Initialize
rootsystem.initialize()

# Simulate
rootsystem.simulate(30, True)

# Export final result (as vtp)
rootsystem.write("../results/example_1a.vtp")

print("done!")
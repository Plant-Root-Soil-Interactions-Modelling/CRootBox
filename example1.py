#
# The Python version of example1.h:
#
#
# 1) Opens plant and root parameters from a file
# 2) Simulates root growth
# 3) Outputs a VTP (for vizualisation in ParaView)
#    In Paraview: use tubePLot.py script for fancy visualisation (Macro/Add new macro...), apply after opening file
#
#  Additionally, exports the line segments as .txt file to import into Matlab for postprocessing
#
import py_rootbox as rb

rootsystem = rb.RootSystem()
name = "Anagallis_femina_Leitner_2010" 

#
# Open plant and root parameter from a file
#
rootsystem.openFile(name)
# rootsystem.writeParameters() # not exposed to python yet

#
# Initialize
#
rootsystem.initialize() 

#
# Simulate
#
simtime = 30  # or 20, 40, 60 days
dt = 1 # try other values here
N = round(simtime/dt) # steps

for i in range(0,int(N)):
    rootsystem.simulate(dt);

#
# Export final result (as vtp)
#
rootsystem.write("results/"+name+".vtp",rb.OutputType.segments) # use ot_polylines for nicer visualization, ot_segments for animations

#
# Export segments for Matlab analysis
#    
analysis = rb.SegmentAnalyser(rootsystem)
analysis.write("results/"+name+".txt")

#
# Total length and surface
#     
l = analysis.getSummed(rb.ScalarType.length)  
print("Total root system length is "+str(l)+" cm")    

print("Finished with a total of "+str(rootsystem.getNumberOfNodes())+ " nodes")

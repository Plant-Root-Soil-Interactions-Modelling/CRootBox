import py_rootbox as rb

rootsystem = rb.RootSystem()
name = "Anagallis_femina_Leitner_2010" 
rootsystem.openFile(name)

# Manually set tropism to hydrotropism for the first ten root types
sigma = [0.4, 1., 1., 1., 1. ] * 2
for i in range(0,10):  
    rootsystem.getRootTypeParameter(i+1).dx = 0.25 # adjust resolution
    rootsystem.getRootTypeParameter(i+1).tropismT = rb.TropismType.hydro
    rootsystem.getRootTypeParameter(i+1).tropismN = 2 # strength of tropism
    rootsystem.getRootTypeParameter(i+1).tropismS = sigma[i] # sigma
     
# Static soil property
maxS = 0.7 # maximal 
minS = 0.1 # minimal 
slope = 5 # linear gradient between min and max (cm)
box = rb.SDF_PlantBox(30,30,2) # cm
layer = rb.SDF_RotateTranslate(box, rb.Vector3d(0,0,-16))
soil_prop = rb.SoilPropertySDF(layer, maxS, minS, slope)

# Set the soil properties before calling initialize
rootsystem.setSoil(soil_prop)

# Initialize
rootsystem.initialize()

# Simulate
simtime = 100 # e.g. 30 or 60 days
dt = 1
N = round(simtime/dt)
for _ in range(0,N):
    # in a dynamic soil setting you would need to update soil_prop        
    rootsystem.simulate(dt)   

# Export results (as vtp)    
rootsystem.write("results/example_4a.vtp")

# Export geometry of static soil  
rootsystem.setGeometry(layer)  # just for vizualisation
rootsystem.write("results/example_4a.py")


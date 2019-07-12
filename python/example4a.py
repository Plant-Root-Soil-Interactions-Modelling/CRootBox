import py_rootbox as rb

rs = rb.RootSystem()
name = "Anagallis_femina_Leitner_2010"
rs.openFile(name)

# Manually set tropism to hydrotropism for the first ten root types
sigma = [0.4, 1., 1., 1., 1. ] * 2
for p in rs.getRootTypeParameter():
        p.dx = 0.25  # adjust resolution
        p.tropismT = rb.TropismType.hydro
        p.tropismN = 2  # strength of tropism
        p.tropismS = sigma[p.subType - 1]

# Static soil property
maxS = 0.7  # maximal
minS = 0.1  # minimal
slope = 5  # linear gradient between min and max (cm)
box = rb.SDF_PlantBox(30, 30, 2)  # cm
layer = rb.SDF_RotateTranslate(box, rb.Vector3d(0, 0, -16))
soil_prop = rb.SoilLookUpSDF(layer, maxS, minS, slope)

# Set the soil properties before calling initialize
rs.setSoil(soil_prop)

# Initialize
rs.initialize()

# Simulate
simtime = 100  # e.g. 30 or 60 days
dt = 1
N = round(simtime / dt)
for _ in range(0, N):
    # in a dynamic soil setting you would need to update soil_prop
    rs.simulate(dt)

# Export results (as vtp)
rs.write("../results/example_4a.vtp")

# Export geometry of static soil
rs.setGeometry(layer)  # just for vizualisation
rs.write("../results/example_4a.py")

print("done.")

import py_rootbox as rb

rs = rb.RootSystem()

# Open plant and root parameter from a file
name = "Anagallis_femina_Leitner_2010"
rs.openFile(name)

# Set Geometry
soilcore = rb.SDF_PlantContainer(5, 5, 40, False)
rs.setGeometry(soilcore)

# Modify axial resolution
for i in range(0, 10):
    p = rs.getRootTypeParameter(i + 1)
    p.dx = 0.1  # adjust resolution

# Simulate
rs.initialize()
rs.simulate(60, True)  # days

# Export results as segments
rb.SegmentAnalyser(rs).write("../results/example_3e.vtp")

# Export container geometry as Paraview Python script
rs.write("../results/example_3e.py")

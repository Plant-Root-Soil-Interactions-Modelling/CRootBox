import vtk
from vtk import *
from IPython.display import Image

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

# Set the background color
colors = vtk.vtkNamedColors()
bkg = map(lambda x: x / 255.0, [26, 51, 102, 255])
colors.SetColor("BkgColor", *bkg)

#
# Xiaoran VTK
#
reader = vtk.vtkXMLPolyDataReader()
path = "../results/example_1a.vtp"  # path or name of the vtp output
reader.SetFileName(path)
reader.Update()
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputData(reader.GetOutput())
plantActor = vtk.vtkActor()
plantActor.SetMapper(mapper)

plantActor.RotateX(-90.0)
# plantActor.RotateY(-45.0)

# Create the graphics structure. The renderer renders into the render
# window. The render window interactor captures mouse events and will
# perform appropriate camera or actor manipulation depending on the
# nature of the events.
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

# Add the actors to the renderer, set the background and size
ren.AddActor(plantActor)
ren.SetBackground(colors.GetColor3d("BkgColor"))
renWin.SetSize(800, 800)
renWin.SetWindowName('test')

# This allows the interactor to initalize itself. It has to be
# called before an event loop.
iren.Initialize()

# We'll zoom in a little by accessing the camera and invoking a "Zoom"
# method on it.
ren.ResetCamera()
# ren.GetActiveCamera().Zoom(1.5)
renWin.Render()

# Start the event loop.
iren.Start()

# renWin = vtk.vtkRenderWindow()
# renWin.AddRenderer(ren1)
# renWin.SetSize(600, 600)
#
# iren = vtk.vtkRenderWindowInteractor()
# iren.SetRenderWindow(renWin)

# ren1 = vtk.vtkRenderer()
# ren1.AddActor(plantActor)
# ren1.SetActiveCamera(camera)
# ren1.SetBackground(0.1, 0.2, 0.4)


from vtk_tools import *
import vtkInterface
import distmesh as dm
import numpy as np
from vtkInterface.grid import VTK_TRIANGLE
from vtk.util.numpy_support import numpy_to_vtkIdTypeArray

fd = lambda p: np.sqrt((p**2).sum(1))-1.0

p, t = dm.distmesh2d(fd, dm.huniform, 0.2, (-1,-1,1,1))

# make the unstructured grid
points = vtkPoints(p)
cells = vtkCells(t)
grid = vtk.vtkUnstructuredGrid()
grid.SetPoints(points)
grid.SetCells(vtk.VTK_TRIANGLE,cells)
# print(cells)
# print(points)

# plot the grid
vtkInterface.Plot(grid)

# mapper = vtk.vtkDataSetMapper()
# mapper.SetInputData(grid)
# actor = vtk.vtkActor()
# actor.SetMapper(mapper)
# 
# # Create a renderer, render window, and interactor
# renderer = vtk.vtkRenderer()
# renderWindow = vtk.vtkRenderWindow()
# renderWindow.AddRenderer(renderer)
# renderWindowInteractor = vtk.vtkRenderWindowInteractor()
# renderWindowInteractor.SetRenderWindow(renderWindow)
#  
# # Add the actor to the scene
# renderer.AddActor(actor)
# renderer.SetBackground(.3, .6, .3) # Background color green
#  
# # Render and interact
# renderWindow.Render()
# renderWindowInteractor.Start()


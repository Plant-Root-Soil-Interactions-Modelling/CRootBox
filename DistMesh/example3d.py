from vtk_tools import *
import vtk
import vtkInterface
import distmesh as dm
import numpy as np

def fd10(p):
    r, z = np.sqrt(p[:,0]**2 + p[:,1]**2), p[:,2]
    d1, d2, d3 = r-1.0, z-1.0, -z-1.0
    d4, d5 = np.sqrt(d1**2+d2**2), np.sqrt(d1**2+d3**2)
    d = dm.dintersect(dm.dintersect(d1, d2), d3)
    ix = (d1>0)*(d2>0); d[ix] = d4[ix]
    ix = (d1>0)*(d3>0); d[ix] = d5[ix]
    return dm.ddiff(d, dm.dsphere(p, 0,0,0, 0.5))

def fh10(p):
    h1 = 4*np.sqrt((p**2).sum(1))-1.0
    return np.minimum(h1, 2.0)

p, t = dm.distmeshnd(fd10, fh10, 0.1, (-1,-1,-1, 1,1,1))
# p = np.zeros((10,3))
# t = np.ones((2,4))

# make the unstructured grid
points = vtkPoints(p)
cells = vtkCells(t)
grid = vtk.vtkUnstructuredGrid()
grid.SetPoints(points)
grid.SetCells(vtk.VTK_TETRA,cells)
# print(cells)
# print(points)

# plot the grid
vtkInterface.Plot(grid)

write_msh("test.msh",grid)
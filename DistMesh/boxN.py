from vtk_tools import *
import vtk
import vtkInterface
import pydistmesh.distmesh as dm
import numpy as np

def fd_box(p,box):
    b = np.zeros(3,)
    mid = np.zeros(3,)
    for i in range(0,3):
        b[i] = 0.5*(box[i+3]-box[i])               
        mid[i] = 0.5*(box[i+3]+box[i])               
    m = np.minimum(b[0]+(p[:,0]-mid[0]),b[0]-(p[:,0]-mid[0]))
    for i in range(1,3):
        m = np.minimum(m,b[i]+(p[:,i]-mid[i]))
        m = np.minimum(m,b[i]-(p[:,i]-mid[i]))
    return -m
 
def fd_unit_ball(p):
    return np.sqrt((p**2).sum(1))-1.0

def fd_plant_container(p,rt=0.05, rb=0.05, h=1, square = False):
    z = p[:,2]/h # scale 0..1
    r =  (1+z)*rb - z*rt
    if square: # rectangular container        
        d = np.maximum(np.abs(p[:,0]),np.abs(p[:,1]))-r
    else: # round pot
        d = np.sqrt(np.square(p[:,0])+np.square(p[:,1]))-r
    return np.maximum(d,-np.minimum(h-p[:,2],0.+p[:,2]))
 
def fh_soil_layer(p, layerZ, hmin = 0.5, hmax = 2., slope = 1):
    h = 1e9
    for l in layerZ:
        h1 = np.maximum(slope*(np.abs(p[:,2]-l)),hmin) 
        h = np.minimum(h, np.minimum(h1, hmax))        
    return h

def fix_bbox(box):
    pfix = np.ones((8,3))
    for i in range(0,8):
        a = i % 2
        b = int(i/2) % 2
        c = int(i/4) % 2
        pfix[i,0] = box[0+a*3] 
        pfix[i,1] = box[1+b*3]   
        pfix[i,2]=  box[2+c*3] 
    return pfix

box = (0, 0, 0, 0.1, 0.2, 2) # b1
fun = lambda p : fd_box(p,box)




#
# BENCHMARK 1
#
cbox = (-0.1,-0.1,0,0.1,0.1,2)
cylinder = lambda p : fd_plant_container(p,0.1,0.1,2)

# p, t = dm.distmeshnd(cylinder, dm.huniform, 0.02, np.array(cbox))
# p, t = rebuild_grid(p,t)
# p = snap_to_box(p, cbox, 1e-5)
# name = "b1_ug.msh" # 9k dof

# fh1 = lambda p: fh_soil_layer(p, [2,1.5], 0.3, 1, 3)
# p, t = dm.distmeshnd(cylinder, fh1, 0.01, np.array(cbox), pfix = np.array([ [0.1,0,0],[-0.1,0,0],[0,0.1,0],[0,-0.1,0], [0.1,0,2],[-0.1,0,2],[0,0.1,2],[0,-0.1,2] ]))
# p, t = rebuild_grid(p,t)
# p = snap_to_box(p, cbox, 1e-5)
# name = "b1_ug2.msh" # 6K dof

#
# BENCHMARK 2
#
cbox2 = (-0.1,-0.1,0,0.1,0.1,.54)
cylinder2 = lambda p : fd_plant_container(p,0.1,0.1,.54)

# p, t = dm.distmeshnd(cylinder2, dm.huniform, 0.01, np.array(cbox2))
# p, t = rebuild_grid(p,t)
# p = snap_to_box(p, cbox2, 1e-5)
# name = "b2_ug.msh" # 17k

fh2 = lambda p: fh_soil_layer(p, [0.54], 0.3, 1, 3)
p, t = dm.distmeshnd(cylinder2, fh2, 0.0075, np.array(cbox2), pfix = np.array([ [0.1,0,0],[-0.1,0,0],[0,0.1,0],[0,-0.1,0], [0.1,0,.54],[-0.1,0,.54],[0,0.1,.54],[0,-0.1,.54] ]))
p, t = rebuild_grid(p,t)
p = snap_to_box(p, cbox2, 1e-5)
name = "b2_ug2.msh" # 11K dof

points = vtkPoints(p) # make the unstructured grid
cells = vtkCells(t)
grid = vtk.vtkUnstructuredGrid()
grid.SetPoints(points)
grid.SetCells(vtk.VTK_TETRA,cells)
celldata=np.hstack((np.ones((grid.GetNumberOfCells(),1)), 10*np.ones((grid.GetNumberOfCells(),1))))

write_msh(name, grid, celldata) 

print()
print("Points: ", p.shape)
print("Triangles", t.shape)
print()
grid_quality(p,t)

vtkInterface.Plot(grid)

print("done")

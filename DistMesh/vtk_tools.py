import vtk
import numpy as np
import matplotlib.pyplot as plt

#
# vtk_tools
#
# D. Leitner, 2018
#



def vtkPoints(p):
    da = vtk.vtkDataArray.CreateDataArray(vtk.VTK_DOUBLE)
    da.SetNumberOfComponents(3) # vtk point dimension is always 3
    da.SetNumberOfTuples(p.shape[0])   
    for i in range(0, p.shape[0]):
        if p.shape[1] == 2:
            da.InsertTuple3(i,p[i,0],p[i,1],0.)
        elif p.shape[1] == 3:        
            da.InsertTuple3(i,p[i,0],p[i,1],p[i,2])
    points = vtk.vtkPoints()      
    points.SetData(da)
    return points



def vtkCells(t):
    cellArray = vtk.vtkCellArray()
    for vert in t:        
        if t.shape[1] == 3:
            tetra = vtk.vtkTriangle()
        elif t.shape[1] == 4:
            tetra = vtk.vtkTetra() 
        for i,v in enumerate(vert):
            tetra.GetPointIds().SetId(i, int(v))        
        cellArray.InsertNextCell(tetra)
    return cellArray



def rebuild_grid(p, t):

    pp = np.zeros(p.shape[0], dtype = "bool") # initially all are unused
    for t_ in t: # find unused points 
        for n in t_:
            pp[n] = 1 # used

    upi = np.argwhere(pp==0) # unused point indices
      
    for k in upi[::-1]: # reverse            
        for i,t_ in enumerate(t): # update triangle indices  
            for j,n in enumerate(t_):
                if n>k:
                    t[i][j] -= 1

    p = np.delete(p, upi, axis=0) # delete unused points
                
    return p,t


def snap_to_box(p, box, eps):    
    for i,p_ in enumerate(p):
        for j in range(0,3):
            if p_[j]<box[j]+eps:
                p[i,j] = box[j]             
        for j in range(3,6):            
            if p_[j-3]>box[j]-eps:
                p[i,j-3] = box[j]             
    return p

def grid_quality(p,t):
    q = np.zeros(t.shape[0])
    for k,t_ in enumerate(t):
        d = p[t_]
        for i in range(0,3):
            d[i] -= d[3]
        q[k] = abs(np.linalg.det(d[0:3,:])/2.)
    print("Quality: min: ",np.min(q), "max:",np.max(q), "mean:",np.mean(q))
    return q


def write_msh(name, ug, celldata = np.zeros((0,0))):
    with open(name, "w") as f:
        # Init
        f.write("$MeshFormat\n")
        f.write("2.2 0 8\n") # version, file-type (ascii=0), data-size
        f.write("$EndMeshFormat\n")
        # Nodes        
        np_ = ug.GetNumberOfPoints()        
        f.write("$Nodes\n")
        f.write("{:d}\n".format(np_)) # number of nodes       
        for i in range(0,np_):
            p = ug.GetPoint(i)
            f.write('{:d} {:08.6f} {:08.6f} {:08.6f}\n'.format(i+1,p[0],p[1],p[2]))    # node number, x, y, z              
        f.write("$EndNodes\n")
        # Cells
        ind = np.zeros(4, dtype=int)
        nc = ug.GetNumberOfCells()
        dn = celldata.shape[1]
        f.write("$Elements\n")
        f.write("{:d}\n".format(nc+8)) # number of cells
        for i in range(0,8):
            f.write("{:d} 15 1 1 {:d}\n".format(i+1,i+1))
        for i in range(0,nc):
            tetra = ug.GetCell(i)
            c = tetra.GetPointIds()
            f.write("{:d} 4 {:d} ".format(i+1+8,dn)) # id, 4 = tetra
            if dn > 0:
                for j in range(0,dn):
                    f.write("{:g} ".format(celldata[i,j]))             
            for j in range(0,4):
                ind[j] = c.GetId(j)+1
#                 # desperation checks
#                 if ind[j] <= 0 :  
#                     print("less euqal Zero !!!!!")
#                 if ind[j] > np_ :
#                     print("value to big")
#                 if ind[j] == 1: 
#                     print("value for 1")
#                 if ind[j] == np_-1: 
#                     print("value for end node")
                
            f.write("{:d} {:d} {:d} {:d}\n".format(ind[0],ind[1],ind[2],ind[3]))
                            
        f.write("$EndElements\n")

#
# opens a vtp and returns the vtk polydata class
#
def read_polydata(name):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(name)
    reader.Update()   
    polydata = reader.GetOutput() 
    return polydata

#
# converts a vtp to dgf
#
def vtp2dgf(name):    
    pd = read_polydata(name+".vtp") # read vtp 

    file = open(name+".dgf","w")   # write dgf
    file.write("DGF\n") 
    # vertex
    file.write('Vertex\n')     
    Np = pd.GetNumberOfPoints()
    points = pd.GetPoints()    
    pdata = pd.GetPointData()
    Npd = pdata.GetNumberOfArrays()
    file.write('parameters {:g}\n'.format(Npd)) 
    for i in range(0,Np):
        p = np.zeros(3,)
        points.GetPoint(i,p)
        file.write('{:g} {:g} {:g} '.format(p[0],p[1],p[2]))
        for j in range(0,Npd): # write point data - todo lets pick ids
            pdataj = pdata.GetArray(j)
            d = pdataj.GetTuple(i)            
            file.write('{:g} '.format(d[0]))
        file.write('\n')

    file.write('#\n');
    file.write('Simplex\n'); 
    # cells
    Nc = pd.GetNumberOfCells() 
    cdata = pd.GetCellData()
    Ncd = cdata.GetNumberOfArrays()    
    file.write('parameters 2\n'.format(Ncd))
    for i in range(0,Nc-1):
        cpi = vtk.vtkIdList()
        pd.GetCellPoints(i,cpi)
        for j in range(0,cpi.GetNumberOfIds()): # write cell ids
            file.write('{:g} '.format(cpi.GetId(j)))
        for j in range(0,Ncd): # write cell data - todo lets pick ids
            cdataj = cdata.GetArray(j)
            d = cdataj.GetTuple(i)            
            file.write('{:g} '.format(d[0]))      
        file.write('\n')
    # i dont know how to use these in dumux        
    file.write('#\n')
    file.write('BOUNDARYSEGMENTS\n') # how do i get the boundary segments into DUMUX ?
    file.write('2 0\n')            
    file.write('3 {:g}\n'.format(Np-1)) # vertex id, but index starts with 0
    file.write('#\n');
    file.write('BOUNDARYDOMAIN\n')
    file.write('default 1\n');
    file.write('#\n')
 
    print("finished writing "+name+".dgf")
    file.close() 
        
#
# returns the cell or vertex data (index 0 and 2 hard coded)) of vtp file 
#
def read1D_vtp_data(name, cell = True):
    polydata = read_polydata(name)     
    
    if cell: 
        data = polydata.GetCellData()
    else:
        data = polydata.GetPointData()        
    
    nocd = data.GetNumberOfArrays()
    p = data.GetArray(0) # pressure   
    noa = p.GetNumberOfTuples()    
    
    p_ = np.ones(noa,)
    for i in range(0,noa):    
        d = p.GetTuple(i)
        p_[i] = d[0]            
    
    return p_

#
# returns the cell or vertex data (index 0 and 2 hard coded)) of vtp file 
#
def read3D_vtp_data(name, cell = True):
    polydata = read_polydata(name)     
    
    if cell: 
        data = polydata.GetCellData()
    else:
        data = polydata.GetPointData()        
    
    nocd = data.GetNumberOfArrays() 
    p = data.GetArray(0) # pressure
    noa = p.GetNumberOfTuples()
    
    p_ = np.ones(noa,)
    for i in range(0,noa):    
        d = p.GetTuple(i)
        p_[i] = d[0]
    
    Np = polydata.GetNumberOfPoints()
    z_ = np.zeros((Np,3))   
    points = polydata.GetPoints()    
    for i in range(0,Np):
        p = np.zeros(3,)
        points.GetPoint(i,p)
        z_[i,:] = p                    
    
    return p_, z_

if __name__ == "__main__":
    # manually set absolute path
    path = "/home/daniel/workspace/DUMUX/dumux-rosi/build-cmake/rosi_benchmarking/richards1d/"
    vtp2dgf(path+"jan1a-00000")
    
    
     

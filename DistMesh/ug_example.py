def meshToUnstructeredGrid(mesh):
    
    """Converts a FiPy mesh structure to a vtkUnstructuredGrid.
    
    Works for 2D and 3D meshes.
    
    Args:
        mesh (fipy.GmshImporter3D): Some Fipy mesh object.
        
    Returns:
        vtk.vtkUnstructuredGrid    
    """
    
    
    # Get vertex coordinates
    coords=mesh.vertexCoords
    
    if len(coords)==2:
        x,y=coords
        dim=2
    else:
        x,y,z=coords
        dim=3
        
    # Insert them as points
    points = vtk.vtkPoints()
    for i in range(len(x)):
        if dim==2:
            points.InsertNextPoint(x[i], y[i],0)
        else:    
            points.InsertNextPoint(x[i], y[i],z[i])

    # Insert tetrahedrons
    verts=mesh._getOrderedCellVertexIDs().T
        
    cellArray = vtk.vtkCellArray()
    for j,vert in enumerate(verts):
        
        if dim==3:
            tetra = vtk.vtkTetra()
        else:
            tetra = vtk.vtkTriangle()
            
        for i,v in enumerate(vert):
            tetra.GetPointIds().SetId(i, v)
        cellArray.InsertNextCell(tetra)

    # Grid
    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(points)
    
    if dim==3:
        grid.SetCells(vtk.VTK_TETRA, cellArray)
    else:
        grid.SetCells(vtk.VTK_TRIANGLE, cellArray)
    
    return grid 
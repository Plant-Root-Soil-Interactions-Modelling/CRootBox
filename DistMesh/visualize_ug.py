#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# by Artem Plotnikov, email: imartemplotnikov <> gmail

import vtk


def read_unstructured_grid(filepath):
    """
    Reads unstructured grid from a VTK file
    :param filepath: Path to *.vtk file, containing unstructured grid
    :return: instanse of vtkUnstructuredGridReader, containing pre-processed grid
    """
#    reader = vtk.vtkUnstructuredGridReader()
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(filepath)
    reader.Update()
    return reader


def apply_quality(grid):
    """
    Applies vtkMeshQuality filter to an input grid, for a more detailed description
    visit http://www.vtk.org/doc/nightly/html/classvtkMeshQuality.html#details
    :param grid: input grid, instance of vtkUnstructuredGrid
    :return: instance of vtkMeshQuality and string, describing used mesh quality metric
    """
    quality = vtk.vtkMeshQuality()
    quality.SetInputData(grid)
    quality.SetTetQualityMeasureToAspectRatio()
    quality.Update()
    return quality, "Aspect Ratio"


def filter_quality(grid, qmin=0.0, qmax=float("inf"), array="Quality"):
    """
    Filter all the cells in the input unstructured grid, which quality measure
    q < qmin or q > qmax.
    :param grid: input grid, instance of vtkUnstructuredGrid
    :param qmin: quality lower bound
    :param qmax: quality upper bound
    :param array: name of mesh quality cellular scalar field. By default, it is the
           same, as vtkMeshQuality produce
    :return:
    """
    threshold = vtk.vtkThreshold()
    threshold.SetInputData(grid)
    threshold.ThresholdBetween(qmin, qmax)
    threshold.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, array)
    threshold.Update()
    return threshold.GetOutput()


def color_unstructured_grid(grid, r, g, b):
    """
    Colors all the cells of input unstructured grid with a desired color
    :param grid: input grid, instance of vtkUnstructuredGrid
    :param r: R part of RGB
    :param g: G part of RGB
    :param b: B part of RGB
    :return:
    """
    colors = vtk.vtkUnsignedCharArray()
    colors.SetName("Colors")
    colors.SetNumberOfComponents(3)
    colors.SetNumberOfTuples(grid.GetNumberOfCells())

    for i in range(0, grid.GetNumberOfCells()):
        colors.InsertTuple3(i, r, g, b)

    grid.GetCellData().SetScalars(colors)


def yield_background_grid_actor(grid):
    color_unstructured_grid(grid, 224, 224, 224)

    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(grid)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetOpacity(0.1)
    actor.GetProperty().LightingOff()
    actor.GetProperty().EdgeVisibilityOn()
    actor.GetProperty().SetEdgeColor(0.36, 0.36, 0.36)
    actor.GetProperty().SetLineWidth(0.3)
    return actor


# see http://www.kennethmoreland.com/color-maps/
def yield_lookup_table():
    lookup_table = vtk.vtkLogLookupTable()
    num_colors = 256
    lookup_table.SetNumberOfTableValues(num_colors)

    transfer_func = vtk.vtkColorTransferFunction()
    transfer_func.SetColorSpaceToDiverging()
    transfer_func.AddRGBPoint(0, 0.230, 0.299,  0.754)
    transfer_func.AddRGBPoint(1, 0.706, 0.016, 0.150)

    for ii,ss in enumerate([float(xx)/float(num_colors) for xx in range(num_colors)]):
        cc = transfer_func.GetColor(ss)
        lookup_table.SetTableValue(ii, cc[0], cc[1], cc[2], 1.0)

    return lookup_table


def yield_filtered_grid_actor(grid, quality_min, quality_max):
    lookup_table = yield_lookup_table()

    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(grid)
    mapper.SetScalarRange(quality_min, quality_max)
    mapper.SetLookupTable(lookup_table)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().EdgeVisibilityOn()

    legend_actor = vtk.vtkScalarBarActor()
    legend_actor.SetLookupTable(lookup_table)

    return actor, legend_actor


def yield_corner_annotation(metric, qmin, qmax, num_cells_showed):
    metric_info = "Quality metric is: %s" % metric
    filter_info = "Filtering criteria is: %s <= quality <= %s. Number of cells found is %s" %\
        (qmin, qmax, num_cells_showed)

    corner_annotation = vtk.vtkCornerAnnotation()
    corner_annotation.SetLinearFontScaleFactor(2)
    corner_annotation.SetNonlinearFontScaleFactor(1)
    corner_annotation.SetMaximumFontSize(20)
    corner_annotation.SetText(0, metric_info)
    corner_annotation.SetText(2, filter_info)
    color = vtk.vtkNamedColors().GetColor3d('warm_grey')
    corner_annotation.GetTextProperty().SetColor(color)
    return corner_annotation


def main():
    
#    filepath = "/daniel/home/workspace/DUMUX/dumux-rosi/build-cmake/rosi_benchmarking/rootsystem/rootsystem-00000.vtp"
    filepath = "/daniel/home/workspace/DUMUX/dumux-rosi/build-cmake/rosi_benchmarking/richards1d/benchmark1a-00000.vtp"

    reader = read_unstructured_grid(filepath)
    grid = reader.GetOutput()
    quality, metric = apply_quality(grid)

    quality_min = quality.GetOutput().GetFieldData()\
        .GetArray('Mesh Tetrahedron Quality').GetComponent(0, 0)
    quality_max = quality.GetOutput().GetFieldData()\
        .GetArray('Mesh Tetrahedron Quality').GetComponent(0, 2)

    qmin = 1.5
    qmax = quality_max

    filtered_grid = filter_quality(quality.GetOutput(), qmin, qmax)
    filtered_grid_actor, legend_actor = yield_filtered_grid_actor(filtered_grid, quality_min, quality_max)

    renderer = vtk.vtkRenderer()
    renderer.AddActor(yield_background_grid_actor(grid))
    renderer.AddActor(filtered_grid_actor)
    renderer.AddActor(legend_actor)
    renderer.SetBackground(1, 1, 1)
    renderer.AddViewProp(yield_corner_annotation(metric, qmin, qmax, filtered_grid.GetNumberOfCells()))

    renderer_window = vtk.vtkRenderWindow()
    renderer_window.AddRenderer(renderer)

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renderer_window)

    interactor_style = vtk.vtkInteractorStyleTrackballCamera()  # ParaView-like interaction
    interactor.SetInteractorStyle(interactor_style)

    interactor.Initialize()
    interactor.Start()


if __name__ == "__main__":
    main()
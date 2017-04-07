#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, MultiLineString
from descartes.patch import PolygonPatch
import vtk

def centrelineToThinPolygon(centreline_polydata, width=0.25e-6):
    segment_collection = []
    for cell_id in range(centreline_polydata.GetNumberOfCells()):
        segment_collection.append([np.array(centreline_polydata.GetCell(cell_id).GetPoints().GetPoint(0)) / 1e6, 
                                   np.array(centreline_polydata.GetCell(cell_id).GetPoints().GetPoint(1)) / 1e6])

    centreline = MultiLineString(segment_collection)
    return centreline.buffer(width)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pixels_per_micron', type=float, help='Pixels per micron in segmented MA image')
    parser.add_argument('ma_body_roi_file', type=str, help='File containing the vertices of polygon delineating the MA body (.csv)')
    parser.add_argument('ma_centreline_file', type=str, help='File defining the MA skeleton or centreline (.vtp)')

    args = parser.parse_args()

    ma_body_roi_file = args.ma_body_roi_file
    ma_centreline_file = args.ma_centreline_file
    pixel_per_metre = 1e6 * args.pixels_per_micron

    # Read MA body delineation from .csv file
    whole_body_coords = np.genfromtxt (ma_body_roi_file, delimiter=",")
    # Create a closed polygon from the delineation by repeating the first vertex at the end
    whole_body_coords = np.append(whole_body_coords, [whole_body_coords[0]], axis=0)
    # Swap axes for consistent XY convention
    whole_body_coords[:, [0, 1]] = whole_body_coords[:, [1, 0]]
    # Turn pixel number into spatial coordinates
    whole_body_coords /= pixel_per_metre

    # Create polygon defining MA body
    ma_polygon = Polygon(whole_body_coords)

    # Create polygon by dilating MA centreline
    vtp_reader = vtk.vtkXMLPolyDataReader()
    vtp_reader.SetFileName(ma_centreline_file)
    vtp_reader.Update()
    vtk_centreline = vtp_reader.GetOutput()
    dilated_centreline = centrelineToThinPolygon(vtk_centreline)

    # Intersect previous two polygons to split MA body in two parts
    ma_hemespheres = ma_polygon.difference(dilated_centreline)

    # Draw MA split polygon 
    fig = plt.figure()
    ax = fig.add_subplot(111)    
    ax.add_patch(PolygonPatch(ma_hemespheres))
    ax.autoscale()
    plt.axis('equal')
    plt.savefig('asymmetry.png')

    # Compute asymmetry ratio
    assert len(ma_hemespheres) == 2, 'Centreline must split MA body in two parts'
    areas = [ma_hemespheres[0].area, ma_hemespheres[1].area]
    print 'Areas: ', areas
    print 'Asymmetry ratio: ', max(areas) / min(areas)

    # Compute body-to-neck ratio
    assert (vtk_centreline.GetPointData().GetNumberOfArrays() == 1) and (vtk_centreline.GetPointData().GetArrayName(0) == 'Radius'), 'Centreline file must contain a point data field named Radius'
    radii_range = vtk_centreline.GetPointData().GetScalars().GetValueRange()
    print 'Radii range: ', radii_range
    print 'Body-to-neck ratio: ', radii_range[1] / radii_range[0]

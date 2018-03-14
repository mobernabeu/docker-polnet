#!/usr/bin/env python

import argparse
import os.path
import numpy as np
from hemeTools.parsers.extraction import ExtractedProperty
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import pandas as pd
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from shapely.geometry import LineString, MultiLineString

flowVariables = {'tangentialprojectiontraction' : 'surface-tractions.xtr',
                 'shearrate' : 'wholegeometry-shearrate.xtr'}

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('results_folder', type=str, help='HemeLb results/ folder location')
    parser.add_argument('pixels_per_micron', type=float, help='Pixels per micron in segmented MA image')
    parser.add_argument('ma_body_roi_file', type=str, help='File containing the vertices of polygon delineating the MA body (.csv)')
    parser.add_argument('ma_centreline_file', type=str, help='File defining the MA skeleton or centreline (.vtp)')
    parser.add_argument('number_subregions', type=int, help='Number of subregions to divide the MA (positive even number)')
    parser.add_argument('subregion_of_interest_file', nargs='?', type=str, help='File containing the vertices of polygon delineating the subregion of interest in the MA (.csv)')

    args = parser.parse_args()

    resultsFolder = args.results_folder
    maBodyRoiFile = args.ma_body_roi_file
    subregionRoiFile = args.subregion_of_interest_file
    pixelPerMetre = 1e6 * args.pixels_per_micron
    ma_centreline_file = args.ma_centreline_file    
    assert (args.number_subregions > 0) and (args.number_subregions % 2 == 0), 'number_subregions argument should be a positive even number' 
    number_subregions_each_side_centreline = args.number_subregions / 2
    
    # Read MA body delineation from .csv file
    wholeBodyCoords = np.genfromtxt (maBodyRoiFile, delimiter=",")
    # Create a closed polygon from the delineation by repeating the first vertex at the end
    wholeBodyCoords = np.append(wholeBodyCoords, [wholeBodyCoords[0]], axis=0)
    # Swap axes for consistent XY convention
    wholeBodyCoords[:, [0, 1]] = wholeBodyCoords[:, [1, 0]]
    # Turn pixel number into spatial coordinates
    wholeBodyCoords /= pixelPerMetre
    # Path object allows inside-the-polygon lookups
    wholeBodyPolygon = Path(wholeBodyCoords)

    if subregionRoiFile != None:
        # Read subregion delineation from .csv file
        subregionCoords = np.genfromtxt (subregionRoiFile, delimiter=",")
        # Create a closed polygon from the delineation by repeating the first vertex at the end
        subregionCoords = np.append(subregionCoords, [subregionCoords[0]], axis=0)
        # Swap axes for consistent XY convention
        subregionCoords[:, [0, 1]] = subregionCoords[:, [1, 0]]
        # Turn pixel number into spatial coordinates
        subregionCoords /= pixelPerMetre
        # Path object allows inside-the-polygon lookups
        subregionPolygon = Path(subregionCoords)
    
    for (variableName, fileName) in flowVariables.iteritems():
        resultsFileName = os.path.join(resultsFolder, 'Extracted', fileName)

        resultsProperties = ExtractedProperty(resultsFileName)
        lastTimestep = resultsProperties.times[-1]
        resultsLastTimeStep = resultsProperties.GetByTimeStep(lastTimestep)

        # Flatten flow domain in 2D
        flowDomainCoords = resultsLastTimeStep.position[:,:2]

        # Visual check for geometrical consistency
        fig = plt.figure()
        ax = fig.add_subplot(111)
        patch = patches.PathPatch(wholeBodyPolygon, facecolor='orange', lw=2)
        ax.add_patch(patch)
        ax.plot(flowDomainCoords[:,0], flowDomainCoords[:,1], 'b+')
        plt.savefig('mean_drop_' + variableName + '.png')

        # Compute mask telling us which results are within the MA body
        ma_mask = wholeBodyPolygon.contains_points(flowDomainCoords)
        assert np.any(ma_mask), 'No results found within the ROI. Check intermediate .png files for overlap with the flow domain'
        
        # If no subregion of interest was provided we compute drops on whole MA
        if subregionRoiFile == None:
            roi_mask = ma_mask
        else:
            # Visual check for geometrical consistency
            fig = plt.figure()
            ax = fig.add_subplot(111)
            patch = patches.PathPatch(subregionPolygon, facecolor='orange', lw=2)
            ax.add_patch(patch)
            ax.plot(flowDomainCoords[:,0], flowDomainCoords[:,1], 'b+')
            plt.savefig('mean_drop_' + variableName + '_subregion.png')

            # Compute mask telling us which results are within the MA body
            roi_mask = subregionPolygon.contains_points(flowDomainCoords)
            assert np.any(roi_mask), 'No results found within the ROI. Check intermediate .png files for overlap with the flow domain'

        results = getattr(resultsLastTimeStep, variableName)
        if variableName in ['tangentialprojectiontraction', 'velocity']:
            results = np.linalg.norm(results, axis=1)

        # Load MA centreline and convert centreline coords to metres
        vtp_reader = vtk.vtkXMLPolyDataReader()
        vtp_reader.SetFileName(ma_centreline_file)
        transform = vtk.vtkTransform()
        transform.Scale(1e-6, 1e-6, 1e-6)
        transform_filter = vtk.vtkTransformFilter()
        transform_filter.SetInputConnection(vtp_reader.GetOutputPort());
        transform_filter.SetTransform(transform);
        transform_filter.Update()
        centreline = transform_filter.GetOutput()

        class CentrelineSideChecker():
            def __init__(self, centreline):
                self.centreline = centreline
                self.crossing_checker = self._centreline_to_MultiLineString(centreline)
                self.xy_inside_bend = np.mean([centreline.GetPoints().GetPoint(centreline_point_id) for centreline_point_id in range(centreline.GetNumberOfPoints())], axis=0)[:2]

            def _centreline_to_MultiLineString(self, centreline):
                segment_collection = []
                for cell_id in range(centreline.GetNumberOfCells()):
                    segment_collection.append([centreline.GetCell(cell_id).GetPoints().GetPoint(0)[:2], 
                                               centreline.GetCell(cell_id).GetPoints().GetPoint(1)[:2]])
                return MultiLineString(segment_collection)

            def is_point_in_concave_side(self, coord):
                xy_coord = coord[:2]
                coord_to_inside_point = LineString([xy_coord, self.xy_inside_bend])
                return not self.crossing_checker.crosses(coord_to_inside_point)

        centreline_side_checker = CentrelineSideChecker(centreline)
 
        # Centreline radii vector in metres
        radii = 1e-6 * vtk_to_numpy(centreline.GetPointData().GetArray('Radius'))

        # Loop over 3D coordinates of points within the roi and perform region classification
        point_coords_within_roi = resultsLastTimeStep.position[roi_mask]
        vtk_points_within_roi = vtk.vtkPoints()
        classification = np.empty(len(point_coords_within_roi))
        for coord_id, coord in enumerate(point_coords_within_roi):

            def unit_distance_to_centreline_point(centreline_point_id):
                centreline_point = centreline.GetPoints().GetPoint(centreline_point_id)
                radius = radii[centreline_point_id]

                distance = np.sqrt(vtk.vtkMath.Distance2BetweenPoints(coord, centreline_point))
                return distance / radius

            centreline_point_id = np.argmin([unit_distance_to_centreline_point(centreline_point_id) for centreline_point_id in range(centreline.GetNumberOfPoints())])
            unit_distance = unit_distance_to_centreline_point(centreline_point_id)

            region_unit_width = 1.0 / number_subregions_each_side_centreline
            region_id = int(np.floor(unit_distance / region_unit_width))
            assert unit_distance < 1.005, 'Unit distance exceeds 1: {}'.format(unit_distance)
            if unit_distance > 1:
                region_id -= 1

            inside_bend = centreline_side_checker.is_point_in_concave_side(coord)
            region_mapping = {True : range(number_subregions_each_side_centreline-1, -1, -1),
                              False : range(number_subregions_each_side_centreline, 2*number_subregions_each_side_centreline)}
            region_id = region_mapping[inside_bend][region_id]

            vtk_points_within_roi.InsertNextPoint(coord);
            classification[coord_id] = region_id

        # Write out .vtp with coordinate classification for debugging purposes
        region_classifier = vtk.vtkPolyData()
        region_classifier.SetPoints(vtk_points_within_roi);
        classification_vtk = numpy_to_vtk(classification)
        classification_vtk.SetName("point_classification");
        region_classifier.GetPointData().AddArray(classification_vtk);
        results_vtk_array = numpy_to_vtk(results[roi_mask].astype(np.float64)) # The cast seems to be required for numpy_to_vtk to work properly with results, a np.float32 array
        results_vtk_array.SetName(variableName)
        region_classifier.GetPointData().AddArray(results_vtk_array)
        polydata_writer = vtk.vtkXMLPolyDataWriter()
        polydata_writer.SetInputData(region_classifier)
        polydata_writer.SetFileName('{}_classification.vtp'.format(variableName))
        polydata_writer.Update()

        # Extract the subset of results contained inside/outside the ROI and calculate drop indices
        mean_roi = np.mean(results[roi_mask])
        mean_non_ma = np.mean(results[np.logical_not(ma_mask)])
        print "{} mean drop ratio: {}".format(variableName, mean_non_ma / mean_roi)
        std_roi = np.std(results[roi_mask])
        std_non_ma = np.std(results[np.logical_not(ma_mask)])
        print "{} std drop ratio: {}".format(variableName, std_non_ma / std_roi)

        # Calculate drop indices on a per region basis
        for region in range(args.number_subregions):
            region_mask = (classification == region)
            if np.sum(region_mask) > 0:
                mean_roi_region = np.mean(results[roi_mask][region_mask])
                print "{} mean drop ratio, region {}: {}".format(variableName, region, mean_non_ma / mean_roi_region)
                std_roi_region = np.std(results[roi_mask][region_mask])
                print "{} std drop ratio, region {}: {}".format(variableName, region, std_non_ma / std_roi_region)

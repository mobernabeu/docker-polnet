#!/usr/bin/env python

import argparse
import os.path
import numpy as np
from hemeTools.parsers.extraction import ExtractedProperty
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import pandas as pd

flowVariables = {'tangentialprojectiontraction' : 'surface-tractions.xtr',
                 'shearrate' : 'wholegeometry-shearrate.xtr'}

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('results_folder', type=str, help='HemeLb results/ folder location')
    parser.add_argument('pixels_per_micron', type=float, help='Pixels per micron in segmented MA image')
    parser.add_argument('region_of_interest_file', type=str, help='File containing the vertices of polygon delineating the region of interest (.csv)')

    args = parser.parse_args()

    resultsFolder = args.results_folder
    roiFile = args.region_of_interest_file
    pixelPerMetre = 1e6 * args.pixels_per_micron
    
    for (variableName, fileName) in flowVariables.iteritems():
        resultsFileName = os.path.join(resultsFolder, 'Extracted', fileName)

        resultsProperties = ExtractedProperty(resultsFileName)
        lastTimestep = resultsProperties.times[-1]
        resultsLastTimeStep = resultsProperties.GetByTimeStep(lastTimestep)

        # Read MA delineation from .csv file
        wholeBodyCoords = np.genfromtxt (roiFile, delimiter=",")
        # Create a closed polygon from the delineation by repeating the first vertex at the end
        wholeBodyCoords = np.append(wholeBodyCoords, [wholeBodyCoords[0]], axis=0)
        # Swap axes for consistent XY convention
        wholeBodyCoords[:, [0, 1]] = wholeBodyCoords[:, [1, 0]]
        # Turn pixel number into spatial coordinates
        wholeBodyCoords /= pixelPerMetre

        # Path object allows inside-the-polygon lookups
        wholeBodyPolygon = Path(wholeBodyCoords)

        # Flatten flow domain in 2D
        flowDomainCoords = resultsLastTimeStep.position[:,:2]

        # Visual check for geometrical consistency
        fig = plt.figure()
        ax = fig.add_subplot(111)
        patch = patches.PathPatch(wholeBodyPolygon, facecolor='orange', lw=2)
        ax.add_patch(patch)
        ax.plot(flowDomainCoords[:,0], flowDomainCoords[:,1], 'b+')
        plt.savefig('mean_drop_' + variableName + '.png')

        # Compute mask telling us which results are with the ROI
        mask = wholeBodyPolygon.contains_points(flowDomainCoords)
        assert np.any(mask), 'No results found within the ROI. Check intermediate .png files for overlap with the flow domain'
        
        results = getattr(resultsLastTimeStep, variableName)
        if variableName in ['tangentialprojectiontraction', 'velocity']:
            results = np.linalg.norm(results, axis=1)

        # Extract the subset of results contained inside/outside the ROI
        mean_roi = np.mean(results[mask])
        mean_non_roi = np.mean(results[np.logical_not(mask)])
        print "{} mean drop ratio: {}".format(variableName, mean_non_roi / mean_roi)

        std_roi = np.std(results[mask])
        std_non_roi = np.std(results[np.logical_not(mask)])
        print "{} std drop ratio: {}".format(variableName, std_non_roi / std_roi)

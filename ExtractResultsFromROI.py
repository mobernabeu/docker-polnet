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
                 'shearrate' : 'wholegeometry-shearrate.xtr',
                 'velocity' : 'wholegeometry-velocity.xtr'}

def SubsampleArray(array, targetLength=1000):
    n = np.size(array) / targetLength
    end =  n * int(len(array)/n)
    average_every_n = np.mean(array[:end].reshape(-1, n), 1)
    return average_every_n[:targetLength]
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('results_folder', type=str, help='HemeLb results folder')
    parser.add_argument('pixel_size', type=float, help='Pixel size in segmented MA image')
    parser.add_argument('region_of_interest_file', type=str, help='File containing the vertices of polygon delineating the region of interest (.csv)')
    parser.add_argument('output_prefix', type=str, help='Prefix for the .xlsx file to be generated as output')

    args = parser.parse_args()

    resultsFolder = args.results_folder
    roiFile = args.region_of_interest_file
    pixelSize = args.pixel_size
    outputPrefix = args.output_prefix

    output_data = {}
    
    for (variableName, fileName) in flowVariables.iteritems():
        resultsFileName = os.path.join(resultsFolder, 'results', 'Extracted', 
                                       fileName)

        # Catch command line arguments
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
        wholeBodyCoords *= pixelSize

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
        plt.savefig(outputPrefix + '_' + variableName + '.png')

        # Compute mask telling us which results are with the ROI
        mask = wholeBodyPolygon.contains_points(flowDomainCoords)
        
        results = getattr(resultsLastTimeStep, variableName)
        if variableName in ['tangentialprojectiontraction', 'velocity']:
            results = np.linalg.norm(results, axis=1)

        # Extract the subset of results contained in the ROI and subsample it
        output_data[variableName] = SubsampleArray(results[mask])

    # Write out data as spreadsheet
    df = pd.DataFrame(output_data)
    df.to_excel(outputPrefix + '.xlsx', index=False)

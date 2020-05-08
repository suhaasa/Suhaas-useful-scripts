# --- VTK-PYTHON SCRIPT FOR READING VTP, INTEGRATING VARIABLES OF INTEREST USING TRAPEZOIDAL RULE
# --- NOTE: THIS SCRIPT ASSUMES A TRIANGULAR SURFACE MESH
# --- BASED ON A SCRIPT BY AEKAANSH VERMA
# --- MODIFIED BY SUHAAS FOR SPLITING BY BRANCH

import sys
import vtk
import numpy as np

# Model that has QOIs data
vtk_file_prefix = 'all_results_mapped_interpolated_threshold_gradients'
# The QOIs to integrate
quantities_to_integrate = ['Threshold Dach1 (20)','Threshold Jag1 (20)','vTAWSS (dynes/cm^2)','Gradients','pressure_avg (mmHg)']
# vtp files of the branches to calculate
branches = ['wall_br_main','wall_br_1','wall_br_2','wall_br_3','wall_br_4','wall_br_1_1','wall_br_3_1', \
'wall_br_3_2','wall_br_4_1']

# First, read in the .vtp file containing your quantities of interest
datareader=vtk.vtkXMLPolyDataReader()
datareader.SetFileName(vtk_file_prefix + '.vtp')
datareader.Update()

# Read your data into another polydata variable for manipulation
model=vtk.vtkPolyData()
model=datareader.GetOutput()
  

output_filename = 'integrated_quantities.dat'
# Integrated quantities
output_collection = []
# Point Averaged quantities
output2_collection = []

for index in xrange(0, len(branches)):

  datareader2=vtk.vtkXMLPolyDataReader()
  datareader2.SetFileName(branches[index] + '.vtp')
  datareader2.Update()
  
  # Read your data into another polydata variable for manipulation
  branch=vtk.vtkPolyData()
  branch=datareader2.GetOutput()
  
  numPts=branch.GetNumberOfPoints()
  numCells=branch.GetNumberOfCells()
  
  # Read in your quantities of interest from the .vtp file
  QOI = []
  for i in xrange(0, len(quantities_to_integrate)):
    QOI.append(vtk.vtkDoubleArray())
    QOI[i] = model.GetPointData().GetArray(quantities_to_integrate[i])

  # Initialize data structures to store area, point total, and point counts
  total_area = 0.0
  integrated_variables = []
  for i in xrange(0, len(quantities_to_integrate)):
    integrated_variables.append(0.0)
  
  point_totals = []
  point_counts = []
  point_averages = []
  for iq in xrange(0, len(quantities_to_integrate)):
    point_totals.append(0.0)
    point_counts.append(0)
    point_averages.append(0.0)

  for pointID in xrange(0, numPts):
    # Get the branch coordinate
    pointCoordinate = branch.GetPoint(pointID) 
    for iq in xrange(0, len(quantities_to_integrate)):
        # Find the branch coordinate in the model
        vector = QOI[iq].GetTuple(model.FindPoint(pointCoordinate))
        mag = 0
        for i in xrange(0,len(vector)):
          mag = mag + float(vector[i])**2
        mag = mag**(.5)
        # Ignore 0 value points
        if(mag>0):
          point_totals[iq] = point_totals[iq] + mag
          point_counts[iq] = point_counts[iq] + 1

  # Calculate point averaged QOIs
  for iq in xrange(0, len(quantities_to_integrate)):
    if(point_counts[iq]>0):
      point_averages[iq] = point_totals[iq]/point_counts[iq]

  # Now loop over the cells and add in contribution to the integrated_variable
  # that is equal to the average of all the cell nodal values multiplied by
  # the area of the cell
  for icell in xrange(0,numCells):
    
    temp_cell = branch.GetCell(icell)
    pts_cell = temp_cell.GetPointIds()
    
    # First, get the area of this cell
    vtkpt = temp_cell.GetPoints()
    p0 = vtkpt.GetPoint(0)
    p1 = vtkpt.GetPoint(1)
    p2 = vtkpt.GetPoint(2)
    temp_area = temp_cell.TriangleArea(p0,p1,p2)
    total_area = total_area + temp_area
    
    # Now, sum up the values of the quantities of interest at the cell
    # vertices
    averages = []
    for iq in xrange(0, len(quantities_to_integrate)):
      averages.append(0.0)
    
    for ipt in xrange(0,pts_cell.GetNumberOfIds()):
      iid = pts_cell.GetId(ipt)
      pointCoordinate = branch.GetPoint(iid) 
      for iq in xrange(0, len(quantities_to_integrate)):
        vector = QOI[iq].GetTuple(model.FindPoint(pointCoordinate))
        mag = 0
        for i in xrange(0,len(vector)):
          mag = mag + float(vector[i])**2
        mag = mag**(.5)
        averages[iq] = averages[iq] + mag
        
    # To complete the trapezoidal rule integration, multiply each summed quantity
    # by the area of the cell, then divide by the number of vertices
    for iq in xrange(0,len(quantities_to_integrate)):
      integrated_variables[iq] = integrated_variables[iq] + \
        averages[iq] * temp_area / float(pts_cell.GetNumberOfIds())
        
  # Now that we have integrated the variables of interest, it is time to save
  # the results to a list for outputting later
  temp_collection = [total_area]
  for iq in xrange(0, len(quantities_to_integrate)):
    temp_collection.append(integrated_variables[iq])

  temp2_collection = []
  for iq in xrange(0, len(quantities_to_integrate)):
    temp2_collection.append(point_averages[iq])
  # Stores the integrated values
  output_collection.append(temp_collection)
  # Stores the point averaged values
  output2_collection.append(temp2_collection)
  
# Now that we have looped over all our .vtp files of interest and integrated
# the variables, it is time to save them to the output file. We also print
# out the integrated quantities divided by the surface area for post-processing
# convenience
  
outfile = open(output_filename, 'w')

# First print a header that tells what each integrated quantity of interest is
out_string = 'Branch, Surface area, '
for iq in xrange(0, len(quantities_to_integrate)):
  out_string = out_string + quantities_to_integrate[iq] + ', '
  out_string = out_string + quantities_to_integrate[iq] + "/Area, "
out_string = out_string + '\n'
outfile.write(out_string)
# Print data
for i in xrange(0, len(branches)):
  out_string = branches[i] + ',' + str(output_collection[i][0])
  for iq in xrange(1, len(quantities_to_integrate)+1):
    out_string = out_string + ', ' + str(output_collection[i][iq])
    out_string = out_string + ', ' + str(output_collection[i][iq]/output_collection[i][0])
    
  out_string = out_string + '\n'
  outfile.write(out_string)

outfile.write('\n')

# Second print a header point averaged quantity of interests
out_string = 'Branch, '
for iq in xrange(0, len(quantities_to_integrate)):
  out_string = out_string + quantities_to_integrate[iq] + ', '
out_string = out_string + '\n'
outfile.write(out_string)
# Print data
for i in xrange(0, len(branches)):
  out_string = branches[i]
  for iq in xrange(0, len(quantities_to_integrate)):
    out_string = out_string + ', ' + str(output2_collection[i][iq])
    
  out_string = out_string + '\n'
  outfile.write(out_string)

outfile.close()
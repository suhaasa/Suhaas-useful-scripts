# --- VTK-PYTHON SCRIPT FOR READING VTP, INTEGRATING VARIABLES OF INTEREST USING TRAPEZOIDAL RULE
# --- NOTE: THIS SCRIPT ASSUMES A TRIANGULAR SURFACE MESH
# --- BASED ON A SCRIPT BY AEKAANSH VERMA
# --- MODIFIED BY SUHAAS FOR SPLITING BY BRANCH

import sys
import vtk
import numpy as np

# Model that has QOIs data
vtk_file_prefix = 'all_results_mapped_interpolated_threshold_gradients'
Mfactors = ['Threshold Dach1 (20)']
Hfactors = ['vTAWSS (dynes/cm^2)','Gradients','pressure_avg (mmHg)']
# The QOIs to integrate
quantities_to_integrate = ['Threshold Dach1 (20)','vTAWSS (dynes/cm^2)','Gradients','pressure_avg (mmHg)']
# vtp files of the branches to calculate
branches = ['walls_combined','wall_br_main','wall_br_1','wall_br_2','wall_br_3','wall_br_4','wall_br_5','wall_br_6', \
'wall_br_7','wall_br_1_1','wall_br_1_2']
num_bins = 5
forced_max = [200, 10, 100]
MF_THRESHOLD = 20
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

output3_collection = np.zeros((len(branches),len(Hfactors),num_bins))

output4_collection  = np.zeros((len(branches),len(Mfactors),len(Hfactors),num_bins))

HF_bins = np.zeros((len(branches),len(Hfactors),num_bins))

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
  QOI_HF = []
  for i in xrange(0, len(Hfactors)):
    QOI_HF.append(vtk.vtkDoubleArray())
    QOI_HF[i] = model.GetPointData().GetArray(Hfactors[i])

  QOI_MF = []
  for i in xrange(0, len(Mfactors)):
    QOI_MF.append(vtk.vtkDoubleArray())
    QOI_MF[i] = model.GetPointData().GetArray(Mfactors[i])

  # Initialize data structures to store area, point total, and point counts
  total_area = 0.0
  integrated_variables = []
  for i in xrange(0, len(Hfactors)):
    integrated_variables.append(0.0)
  
  point_totals = []
  point_counts = []
  point_averages = []
  point_maxs = []
  for iq in xrange(0, len(Hfactors)):
    point_totals.append(0.0)
    point_counts.append(0)
    point_averages.append(0.0)
    point_maxs.append(0.0)

  for pointID in xrange(0, numPts):
    # Get the branch coordinate
    pointCoordinate = branch.GetPoint(pointID) 
    for iq in xrange(0, len(Hfactors)):
        # Find the branch coordinate in the model
        vector = QOI_HF[iq].GetTuple(model.FindPoint(pointCoordinate))
        mag = 0
        for i in xrange(0,len(vector)):
          mag = mag + float(vector[i])**2
        mag = mag**(.5)
        # Ignore 0 value points
        if(mag>0):
          point_totals[iq] = point_totals[iq] + mag
          point_counts[iq] = point_counts[iq] + 1
        if(mag>point_maxs[iq]):
          point_maxs[iq]=mag

  
  for HF in xrange(0,len(Hfactors)):
    for bin_i in xrange(0,num_bins):
      max = point_maxs[HF]
      if (max<forced_max[HF]):
        HF_bins[index][HF][bin_i] = bin_i * max/num_bins
      else:
        HF_bins[index][HF][bin_i] = bin_i * float(forced_max[HF]/num_bins)

  HF_areas = np.zeros((len(Hfactors),num_bins))
  MF_pos_areas = np.zeros((len(Mfactors),len(Hfactors),num_bins))
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
    cell_area = temp_cell.TriangleArea(p0,p1,p2)
    total_area = total_area + cell_area
    
    
    # Now, sum up the values of the quantities of interest at the cell
    # vertices
    HF_averages = []
    for iq in xrange(0, len(Hfactors)):
      HF_averages.append(0.0)
    MF_averages = []
    for MF in xrange(0,len(Mfactors)):
      MF_averages.append(0.0)

    for ipt in xrange(0,pts_cell.GetNumberOfIds()):
      iid = pts_cell.GetId(ipt)
      pointCoordinate = branch.GetPoint(iid) 
      for iq in xrange(0, len(Hfactors)):
        vector = QOI_HF[iq].GetTuple(model.FindPoint(pointCoordinate))
        mag = 0
        for i in xrange(0,len(vector)):
          mag = mag + float(vector[i])**2
        mag = mag**(.5)
        HF_averages[iq] = HF_averages[iq] + mag
      for MF in xrange(0,len(Mfactors)):
        val = QOI_MF[MF].GetValue(model.FindPoint(pointCoordinate))
        MF_averages[MF] = MF_averages[MF] + val

    for iq in xrange(0,len(Hfactors)):
      HF_averages[iq] = HF_averages[iq] / float(pts_cell.GetNumberOfIds())
      for bin_i in xrange(0,num_bins-1):
        if (HF_averages[iq]>HF_bins[index][iq][bin_i] and HF_averages[iq]<HF_bins[index][iq][bin_i+1]):
          HF_areas[iq][bin_i] = HF_areas[iq][bin_i] + cell_area
      if(HF_averages[iq]>HF_bins[index][iq][num_bins-1]):
        HF_areas[iq][num_bins-1] = HF_areas[iq][num_bins-1] + cell_area
      for MF in xrange(0,len(Mfactors)):
        if(MF_averages[MF]>MF_THRESHOLD):
          for bin_i in xrange(0,num_bins-1):
            if (HF_averages[iq]>HF_bins[index][iq][bin_i] and HF_averages[iq]<HF_bins[index][iq][bin_i+1]):
              MF_pos_areas[MF][iq][bin_i] = MF_pos_areas[MF][iq][bin_i] + cell_area
          if(HF_averages[iq]>HF_bins[index][iq][num_bins-1]):
            MF_pos_areas[MF][iq][num_bins-1] = MF_pos_areas[MF][iq][num_bins-1] + cell_area


  # Now that we have integrated the variables of interest, it is time to save
  # the results to a list for outputting later
  temp3_collection = []
  temp3_bins= [total_area]
  for HF in xrange(0,len(Hfactors)):
    for iq in xrange(0, num_bins):
      output3_collection[index][HF][iq] = HF_areas[HF][iq]

  for MF in xrange(0,len(Mfactors)):
    for HF in xrange(0,len(Hfactors)):
      for iq in xrange(0, num_bins):
        output4_collection[index][MF][HF][iq] = MF_pos_areas[MF][HF][iq]
# Now that we have looped over all our .vtp files of interest and integrated
# the variables, it is time to save them to the output file. We also print
# out the integrated quantities divided by the surface area for post-processing
# convenience
print(HF_bins)
outfile = open(output_filename, 'w')

out_string = 'Branch, Total Area, '
for iq in xrange(0, num_bins):
  out_string = out_string + 'bin ' + str(iq) + ', '
out_string = out_string + '\n'
outfile.write(out_string)
# Print data
for i in xrange(0, len(branches)):
  out_string = branches[i]
  for HF in xrange(0,len(Hfactors)):
    out_string = out_string + ', ' + Hfactors[HF]
    for iq in xrange(0, num_bins):
      out_string = out_string + ', ' + str(HF_bins[i][HF][iq])
    out_string = out_string + '\n'
    out_string = out_string + ' , Areas'
    for iq in xrange(0, num_bins):
      out_string = out_string + ', ' + str(output3_collection[i][HF][iq])
    out_string = out_string + '\n'

  for MF in xrange(0, len(Mfactors)):
    out_string = out_string + ', ' + Mfactors[MF]
    out_string = out_string + '\n'
    for HF in xrange(0,len(Hfactors)):
      out_string = out_string + ', ' + Hfactors[HF]
      for iq in xrange(0, num_bins):
        out_string = out_string + ', ' + str(output4_collection[i][MF][HF][iq])
      out_string = out_string + '\n'
  out_string = out_string + '\n'
  outfile.write(out_string)

outfile.close()
# --- VTK-PYTHON SCRIPT FOR READING VTP, INTEGRATING VARIABLES OF INTEREST USING TRAPEZOIDAL RULE
# --- NOTE: THIS SCRIPT ASSUMES A TRIANGULAR SURFACE MESH
# --- BASED ON A SCRIPT BY AEKAANSH VERMA

import sys
import vtk
import numpy as np
from vtk.util import numpy_support as VN

vtk_file_prefix = 'all_results'
quantities_to_integrate = ['vTAWSS']
CONV_FACTOR = 10000
output_filename = 'TAWSS_converted_results.dat'
output_collection = []
# First, read in the .vtp file containing your quantities of interest
datareader=vtk.vtkXMLPolyDataReader()
datareader.SetFileName(vtk_file_prefix + '.vtp')
datareader.Update()

# Read your data into another polydata variable for manipulation
model=vtk.vtkPolyData()
model=datareader.GetOutput()

numPts=model.GetNumberOfPoints()
numCells=model.GetNumberOfCells()

# Read in your quantities of interest from the .vtp file
TA_WSS = VN.vtk_to_numpy(model.GetPointData().GetArray('vTAWSS'))
QOI = []
for i in xrange(0, len(quantities_to_integrate)):
  QOI.append(vtk.vtkDoubleArray())
  QOI[i] = model.GetPointData().GetArray(quantities_to_integrate[i])


for icell in xrange(0,numCells):
    
    temp_cell = model.GetCell(icell)
    pts_cell = temp_cell.GetPointIds()
    
    
    for ipt in xrange(0,pts_cell.GetNumberOfIds()):
      iid = pts_cell.GetId(ipt)
#      for iq in xrange(0, len(quantities_to_integrate)):
#        output_collection = CONV_FACTOR * float(QOI[iq].GetTuple(iid)[0])
TA_WSS_numpy = 10000*TA_WSS
TA_WSS_vtk = VN.numpy_to_vtk(TA_WSS_numpy)
TA_WSS_vtk.SetName('vTAWSS (dynes/cm^2)') #rememebr to give an unique name
model.GetPointData().AddArray(TA_WSS_vtk) 


  
# Now that we have looped over all our .vtp files of interest and integrated
# the variables, it is time to save them to the output file. We also print
# out the integrated quantities divided by the surface area for post-processing
# convenience

w = vtk.vtkPolyDataWriter()
#w.SetInput(tpd1.GetOutput())
w.SetInputData(model)
w.SetFileName("test.vtk")
w.Write()
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      

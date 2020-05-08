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
TA_WSS_numpy = 10000*TA_WSS
TA_WSS_vtk = VN.numpy_to_vtk(TA_WSS_numpy)
TA_WSS_vtk.SetName('vTAWSS (dynes/cm^2)') 
model.GetPointData().AddArray(TA_WSS_vtk) 

pressure_avg = VN.vtk_to_numpy(model.GetPointData().GetArray('pressure_avg'))
pressure_numpy = pressure_avg/.13333333333333
pressure_vtk = VN.numpy_to_vtk(pressure_numpy)
pressure_vtk.SetName('pressure_avg (mmHg)') 
model.GetPointData().AddArray(pressure_vtk) 

  
# Now that we have looped over all our .vtp files of interest and integrated
# the variables, it is time to save them to the output file. We also print
# out the integrated quantities divided by the surface area for post-processing
# convenience

w = vtk.vtkPolyDataWriter()
#w.SetInput(tpd1.GetOutput())
w.SetInputData(model)
w.SetFileName('all_results_converted.vtk')
w.Write()
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      

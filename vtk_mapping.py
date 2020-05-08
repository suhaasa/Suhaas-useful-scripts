# --- VTK-PYTHON SCRIPT FOR READING MESH VTP AND DACH1 VTI TO MAP DACH1 VALUES ONTO MESH SURFACE
# --- WRITTEN BY SUHAAS ANBAZHAKAN
# --- BASED ON A SCRIPT BY AEKAANSH VERMA

import sys
import vtk
import numpy as np
from vtk.util import numpy_support as VN

vtk_file_prefix = 'all_results'
WSS_CONV_FACTOR = 10000
P_CONV_FACTOR = .133333333333333
INTERPOLATING_RADIUS = 5
# First, read in the .vtp file for mesh and Dach1
datareader=vtk.vtkXMLPolyDataReader()
datareader.SetFileName(vtk_file_prefix + '.vtp')
datareader.Update()

ref = vtk.vtkXMLImageDataReader()
ref.SetFileName('Dach1.vti')
ref.Update()

# Read your data into another polydata variable for reading
mesh=vtk.vtkPolyData()
mesh=datareader.GetOutput()

dach=vtk.vtkPolyData()
dach=ref.GetOutput()

numPts=mesh.GetNumberOfPoints()

# Create new Dach1 Array to fill
dach1Array = vtk.vtkDoubleArray()

# Loop through mesh points, find closest point in Dach1.vti, and add to Dach1 Array
for meshPointID in xrange(0, numPts):
	meshPointCoordinate = mesh.GetPoint(meshPointID)
	dachPointID = dach.FindPoint(meshPointCoordinate)
	dach1Value = dach.GetPointData().GetArray('ImageFile').GetValue(dachPointID)
	dach1Array.InsertNextValue(dach1Value)

# Add Dach1 Array to point data
dach1Array.SetName('Dach1') 
mesh.GetPointData().AddArray(dach1Array)

# Convert TA_WSS and Pressure to micron units and new arrays to point data
TA_WSS = VN.vtk_to_numpy(model.GetPointData().GetArray('vTAWSS'))
TA_WSS_numpy = WSS_CONV_FACTOR*TA_WSS
TA_WSS_vtk = VN.numpy_to_vtk(TA_WSS_numpy)
TA_WSS_vtk.SetName('vTAWSS (dynes/cm^2)') 
mesh.GetPointData().AddArray(TA_WSS_vtk) 

pressure_avg = VN.vtk_to_numpy(model.GetPointData().GetArray('pressure_avg'))
pressure_numpy = pressure_avg/P_CONV_FACTOR
pressure_vtk = VN.numpy_to_vtk(pressure_numpy)
pressure_vtk.SetName('pressure_avg (mmHg)') 
mesh.GetPointData().AddArray(pressure_vtk) 


# Write a new .vtp file that has all converted values and mapped values.
w = vtk.vtkXMLPolyDataWriter()
w.SetInputData(model)
w.SetFileName('all_results_mapped.vtp')
w.Write()
      
      
      
      
      
      
      
      
      
     
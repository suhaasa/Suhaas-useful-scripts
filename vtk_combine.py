# --- VTK-PYTHON SCRIPT FOR READING MESH VTP AND DACH1 VTI TO MAP DACH1 VALUES ONTO MESH SURFACE
# --- WRITTEN BY SUHAAS ANBAZHAKAN
# --- BASED ON A SCRIPT BY AEKAANSH VERMA

import sys
import vtk
import numpy as np
from vtk.util import numpy_support as VN

vtk_file_prefix = 'all_results'
WSS_CONV_FACTOR = 1
P_CONV_FACTOR = .133333333333333
INTERPOLATING_RADIUS = 5
factors = ['Dach1']
# First, read in the .vtp file for mesh and Dach1
datareader=vtk.vtkXMLImageDataReader()
datareader.SetFileName('test.vti')
datareader.Update()

mesh=vtk.vtkPolyData()
mesh=datareader.GetOutput()

for MF in xrange(0,len(factors)):
	ref = vtk.vtkXMLImageDataReader()
	ref.SetFileName(factors[MF] + '.vti')
	ref.Update()

	# Read your data into another polydata variable for reading

	dach=vtk.vtkPolyData()
	dach=ref.GetOutput()

	# Convert Dach1 and Pressure to micron units and new arrays to point data
	Dach1 = VN.vtk_to_numpy(dach.GetPointData().GetArray('ImageFile'))
	Dach1_numpy = WSS_CONV_FACTOR*Dach1
	Dach1_vtk = VN.numpy_to_vtk(Dach1_numpy)
	Dach1_vtk.SetName(factors[MF]) 
	mesh.GetPointData().AddArray(Dach1_vtk) 


# Write a new .vtp file that has all converted values and mapped values.
w = vtk.vtkXMLImageDataWriter()
w.SetInputData(mesh)
w.SetFileName('E175_H2_ImageData_Combined.vti')
w.Write()
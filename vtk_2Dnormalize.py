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
GLOBAL_AVERAGE = 160
a = 2
b = 3
# First, read in the .vtp file for mesh and Dach1
datareader=vtk.vtkXMLImageDataReader()
datareader.SetFileName('test.vti')
datareader.Update()

ref = vtk.vtkXMLImageDataReader()
ref.SetFileName('E175_H2_Dach1_ImageData_Combined.vti')
ref.Update()

# Read your data into another polydata variable for reading
mesh=vtk.vtkPolyData()
mesh=datareader.GetOutput()

imageData=vtk.vtkPolyData()
imageData=ref.GetOutput()

numPts=imageData.GetNumberOfPoints()
sum = 0
for imageDataPointID in xrange(0, numPts):
	imageDataPointCoordinate = imageData.GetPoint(imageDataPointID)
	x,y,z = imageData.GetPoint(imageDataPointID)
	kernel_size = np.array([a,b]);
	sum = sum + imageData.GetPointData().GetArray('Dach1').GetValue(imageDataPointID)
GLOBAL_AVERAGE = sum/numPts

normalizedDach1Array = vtk.vtkDoubleArray()

for imageDataPointID in xrange(0, numPts):
	imageDataPointCoordinate = imageData.GetPoint(imageDataPointID)
	x,y,z = imageData.GetPoint(imageDataPointID)
	kernel_size = np.array([a,b]);
	initialDach1Value = imageData.GetPointData().GetArray('Dach1').GetValue(imageDataPointID)
	sum = 0
	count = 1
	for kernel_x in xrange(-a, a):
		for kernel_y in xrange(-b, b):
			x = x - a
			y = y - b
			kernelPointCoordinate = x,y,z
			imageDataPointID = imageData.FindPoint(kernelPointCoordinate)
			sum = sum +  imageData.GetPointData().GetArray('Dach1').GetValue(kernelPointID)
			count = count + 1

	local_average = sum/count
	normalizedDach1 = GLOBAL_AVERAGE/local_average*initialDach1Value
	normalizedDach1Array.InsertNextValue(normalizedDach1)



# Convert Dach1 and Pressure to micron units and new arrays to point data
maxDach1Array.SetName('Normalized Dach1 (' + str(a) + ',' + str(b) + ')') 
mesh.GetPointData().AddArray(normalizedDach1Array)


# Write a new .vtp file that has all converted values and mapped values.
w = vtk.vtkXMLImageDataWriter()
w.SetInputData(mesh)
w.SetFileName('E175_H2_Dach1_ImageData_Normalized.vti')
w.Write()
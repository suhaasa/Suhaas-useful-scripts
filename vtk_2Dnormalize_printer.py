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
GLOBAL_AVERAGE = 80
THRESHOLD = 20
a = 50
b = 50
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

normalizedDach1Array = vtk.vtkDoubleArray()

ImageData_numpy = VN.vtk_to_numpy(imageData.GetPointData().GetArray('Dach1'))
GLOBAL_AVERAGE = np.mean(ImageData_numpy)
GLOBAL_max = np.max(ImageData_numpy)

sum = 0
count = 0
for imageDataPointID in xrange(0, numPts):
	imageDataPointCoordinate = imageData.GetPoint(imageDataPointID)
	x,y,z = imageData.GetPoint(imageDataPointID)
	kernel_size = np.array([a,b]);
	if(imageData.GetPointData().GetArray('Dach1').GetValue(imageDataPointID)>THRESHOLD):
		sum = sum + imageData.GetPointData().GetArray('Dach1').GetValue(imageDataPointID)
		count = count + 1
GLOBAL_AVERAGE = sum/numPts

for imageDataPointID in xrange(0, numPts):
	imageDataPointCoordinate = imageData.GetPoint(imageDataPointID)
	x,y,z = imageData.GetPoint(imageDataPointID)
	kernel_size = np.array([a,b]);
	initialDach1Value = imageData.GetPointData().GetArray('Dach1').GetValue(imageDataPointID)
	if(imageDataPointID % 100 ==0):
		print(initialDach1Value)
	sum = 0
	count = 1
	percentDone = float(imageDataPointID/numPts*100)
	for kernel_x in xrange(-a, a):
		for kernel_y in xrange(-b, b):
			if(x-a>0 and x+a<1000 and y-b>0 and y+b<1000):
				x = x - a
				y = y - b
				kernelPointCoordinate = x,y,z
				kernelPointID = imageData.FindPoint(kernelPointCoordinate)
				sum = sum +  imageData.GetPointData().GetArray('Dach1').GetValue(kernelPointID)
				count = count + 1
				

	local_average = sum/count
	if(local_average>0):
		normalizedDach1 = GLOBAL_AVERAGE/local_average*initialDach1Value
		normalizedDach1Array.InsertNextValue(normalizedDach1)



# Convert Dach1 and Pressure to micron units and new arrays to point data
normalizedDach1Array.SetName('Normalized Dach1 (' + str(a) + ',' + str(b) + ')') 
imageData.GetPointData().AddArray(normalizedDach1Array)


# Write a new .vtp file that has all converted values and mapped values.
w = vtk.vtkXMLImageDataWriter()
w.SetInputData(imageData)
w.SetFileName('E175_H2_Dach1_ImageData_Normalized.vti')
w.Write()
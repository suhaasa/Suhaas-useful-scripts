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
INTERPOLATING_RADIUS = 10
# First, read in the .vtp file for mesh and Dach1
datareader=vtk.vtkXMLPolyDataReader()
datareader.SetFileName(vtk_file_prefix + '.vtp')
datareader.Update()

ref = vtk.vtkXMLImageDataReader()
ref.SetFileName('E175_H2_Dach1_ImageData_Combined.vti')
ref.Update()

# Read your data into another polydata variable for reading
mesh=vtk.vtkPolyData()
mesh=datareader.GetOutput()

dach=vtk.vtkPolyData()
dach=ref.GetOutput()

numPts=mesh.GetNumberOfPoints()

# Create new Dach1 Array to fill
maxDach1Array = vtk.vtkDoubleArray()
meanDach1Array = vtk.vtkDoubleArray()
thresholdDach1Array = vtk.vtkDoubleArray()

# Loop through mesh points, find closest point in Dach1.vti, and add to Dach1 Array
for meshPointID in xrange(0, numPts):
	meshPointCoordinate = mesh.GetPoint(meshPointID)
	x,y,z = mesh.GetPoint(meshPointID)
	pointsPool = np.array([[x,y,z]]);
	pointsPool = np.append(pointsPool,[[x+INTERPOLATING_RADIUS,y,z],[x-INTERPOLATING_RADIUS,y,z]],axis=0)
	pointsPool = np.append(pointsPool,[[x,y+INTERPOLATING_RADIUS,z],[x,y-INTERPOLATING_RADIUS,z]],axis=0)
	pointsPool = np.append(pointsPool,[[x,y,z+INTERPOLATING_RADIUS],[x,y,z-INTERPOLATING_RADIUS]],axis=0)
	maxDach1Value=0
	meanDach1Value=0
	threshold = 1
	x = pointsPool[0,0]
	y = pointsPool[0,1]
	z = pointsPool[0,2]
	meshPointCoordinate = x,y,z
	thresholdCx40Value = dach.GetPointData().GetArray('Scalars_').GetValue(dach.FindPoint(meshPointCoordinate))
	for i in xrange(0,len(pointsPool)):
		x = pointsPool[i,0]
		y = pointsPool[i,1]
		z = pointsPool[i,2]
		meshPointCoordinate = x,y,z
		dachPointID = dach.FindPoint(meshPointCoordinate)
		dach1Value = dach.GetPointData().GetArray('Dach1').GetValue(dachPointID)
		meanDach1Value +=dach1Value
		Cx40Value = dach.GetPointData().GetArray('Scalars_').GetValue(dachPointID)
		if maxDach1Value<dach1Value:
			maxDach1Value = dach1Value
		if (Cx40Value>thresholdCx40Value):
			threshold = 0

	if threshold==1:
		x = pointsPool[0,0]
		y = pointsPool[0,1]
		z = pointsPool[0,2]
		meshPointCoordinate = x,y,z
		thresholdDach1Array.InsertNextValue(dach.GetPointData().GetArray('Dach1').GetValue(dach.FindPoint(meshPointCoordinate)))
	else:
		thresholdDach1Array.InsertNextValue(1)

	meanDach1Value = meanDach1Value/len(pointsPool)
	maxDach1Array.InsertNextValue(maxDach1Value)
	meanDach1Array.InsertNextValue(meanDach1Value)
# Add Dach1 Array to point data
maxDach1Array.SetName('Max Dach1 (' + str(INTERPOLATING_RADIUS) + ')') 
mesh.GetPointData().AddArray(maxDach1Array)
meanDach1Array.SetName('Mean Dach1 (' + str(INTERPOLATING_RADIUS)+ ')') 
mesh.GetPointData().AddArray(meanDach1Array)
thresholdDach1Array.SetName('Threshold Dach1 (' + str(INTERPOLATING_RADIUS)+ ')') 
mesh.GetPointData().AddArray(thresholdDach1Array)

# Convert TA_WSS and Pressure to micron units and new arrays to point data
TA_WSS = VN.vtk_to_numpy(mesh.GetPointData().GetArray('vTAWSS'))
TA_WSS_numpy = WSS_CONV_FACTOR*TA_WSS
TA_WSS_vtk = VN.numpy_to_vtk(TA_WSS_numpy)
TA_WSS_vtk.SetName('vTAWSS (dynes/cm^2)') 
mesh.GetPointData().AddArray(TA_WSS_vtk) 

pressure_avg = VN.vtk_to_numpy(mesh.GetPointData().GetArray('pressure_avg'))
pressure_numpy = pressure_avg/P_CONV_FACTOR
pressure_vtk = VN.numpy_to_vtk(pressure_numpy)
pressure_vtk.SetName('pressure_avg (mmHg)') 
mesh.GetPointData().AddArray(pressure_vtk) 


# Write a new .vtp file that has all converted values and mapped values.
w = vtk.vtkXMLPolyDataWriter()
w.SetInputData(mesh)
w.SetFileName('all_results_mapped_interpolated_threshold.vtp')
w.Write()
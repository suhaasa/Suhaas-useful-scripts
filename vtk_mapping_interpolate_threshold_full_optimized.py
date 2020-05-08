# --- VTK-PYTHON SCRIPT FOR READING MESH VTP AND DACH1 VTI TO MAP DACH1 VALUES ONTO MESH SURFACE
# --- WRITTEN BY SUHAAS ANBAZHAKAN
# --- BASED ON A SCRIPT BY AEKAANSH VERMA

import sys
import vtk
import numpy as np
from vtk.util import numpy_support as VN
def convertData(model, oldDataName, newdataName, factor):
	data = VN.vtk_to_numpy(model.GetPointData().GetArray(oldDataName))
	numpy = factor*data
	vtk = VN.numpy_to_vtk(numpy)
	vtk.SetName(newdataName) 
	mesh.GetPointData().AddArray(vtk) 

vtk_file_prefix = 'all_results'
WSS_CONV_FACTOR = 10000
P_CONV_FACTOR = .133333333333333
INTERPOLATING_RADIUS = 10
RANGE = 10
CX40THRESHOLD = 50
# First, read in the .vtp file for mesh and Dach1
datareader=vtk.vtkXMLPolyDataReader()
datareader.SetFileName(vtk_file_prefix + '.vtp')
datareader.Update()

ref = vtk.vtkXMLImageDataReader()
ref.SetFileName('E175_H2_ImageData_Combined.vti')
ref.Update()

# Read your data into another polydata variable for reading
mesh=vtk.vtkPolyData()
mesh=datareader.GetOutput()

imageData=vtk.vtkPolyData()
imageData=ref.GetOutput()

meshNormals = vtk.vtkPolyDataNormals()
meshNormals.SetInputData(mesh)
meshNormals.AutoOrientNormalsOn()
meshNormals.ComputePointNormalsOn()
meshNormals.Update()

numPts=mesh.GetNumberOfPoints()

factors = ['Dach1']

for molecularFactor in xrange(0, len(factors)):
	# Create new Dach1 Array to fill
	maxMFArray = vtk.vtkDoubleArray()
	meanMFArray = vtk.vtkDoubleArray()
	thresholdMFArray = vtk.vtkDoubleArray()

	# Loop through mesh points, find closest point in Dach1.vti, and add to Dach1 Array

	for meshPointID in xrange(0, numPts):
		meshPointCoordinate = mesh.GetPoint(meshPointID)
		x,y,z = mesh.GetPoint(meshPointID)
		pointsPool = np.array([[x,y,z]]);
		pointsPool = np.append(pointsPool,[[x+INTERPOLATING_RADIUS,y,z],[x-INTERPOLATING_RADIUS,y,z]],axis=0)
		pointsPool = np.append(pointsPool,[[x,y+INTERPOLATING_RADIUS,z],[x,y-INTERPOLATING_RADIUS,z]],axis=0)
		pointsPool = np.append(pointsPool,[[x,y,z+INTERPOLATING_RADIUS],[x,y,z-INTERPOLATING_RADIUS]],axis=0)
		maxMFValue=0
		meanMFValue=0
		x0 = pointsPool[0,0]
		y0 = pointsPool[0,1]
		z0 = pointsPool[0,2]
		meshPointCoordinate0 = x0,y0,z0
		thresholdCx40Value = imageData.GetPointData().GetArray('Scalars_').GetValue(imageData.FindPoint(meshPointCoordinate0))
		for i in xrange(0,len(pointsPool)):
			x = pointsPool[i,0]
			y = pointsPool[i,1]
			z = pointsPool[i,2]
			meshPointCoordinate = x,y,z
			meshPointID = imageData.FindPoint(meshPointCoordinate)
			MFValue = imageData.GetPointData().GetArray(factors[molecularFactor]).GetValue(meshPointID)
			meanMFValue +=MFValue

			Cx40Value = imageData.GetPointData().GetArray('Scalars_').GetValue(meshPointID)

			if maxMFValue<MFValue:
				maxMFValue = MFValue

		# if thresholdCx40Value>CX40THRESHOLD:
		# 	thresholdMFArray.InsertNextValue(imageData.GetPointData().GetArray(factors[molecularFactor]).GetValue(imageData.FindPoint(meshPointCoordinate0)))
		# else:
		# 	thresholdMFArray.InsertNextValue(0)

		meshPointNormal = meshNormals.GetOutput().GetPointData().GetNormals().GetTuple(mesh.FindPoint(meshPointCoordinate0))
		maxCx40PointID = imageData.FindPoint(meshPointCoordinate0)
		maxCx40Value = imageData.GetPointData().GetArray(factors[molecularFactor]).GetValue(imageData.FindPoint(meshPointCoordinate0))
		for i in xrange(0,RANGE):
			x = meshPointCoordinate0[0] - i/10 * meshPointNormal[0]
			y = meshPointCoordinate0[1] - i/10 * meshPointNormal[1]
			z = meshPointCoordinate0[2] - i/10 * meshPointNormal[2]
			meshPointCoordinateNormal = x,y,z
			Cx40ValueNormal = imageData.GetPointData().GetArray(factors[molecularFactor]).GetValue(imageData.FindPoint(meshPointCoordinateNormal))
			if(Cx40ValueNormal>maxCx40Value):
				maxCx40PointID = imageData.FindPoint(meshPointCoordinateNormal)
				maxCx40Value = Cx40ValueNormal
			x = meshPointCoordinate0[0] + i/10 * meshPointNormal[0]
			y = meshPointCoordinate0[1] + i/10 * meshPointNormal[1]
			z = meshPointCoordinate0[2] + i/10 * meshPointNormal[2]
			meshPointCoordinateNormal = x,y,z
			Cx40ValueNormal = imageData.GetPointData().GetArray(factors[molecularFactor]).GetValue(imageData.FindPoint(meshPointCoordinateNormal))
			if(Cx40ValueNormal>maxCx40Value):
				maxCx40PointID = imageData.FindPoint(meshPointCoordinateNormal)
				maxCx40Value = Cx40ValueNormal
		if(maxCx40Value>20):
			MF_val = imageData.GetPointData().GetArray(factors[molecularFactor]).GetValue(maxCx40PointID)
			thresholdMFArray.InsertNextValue(MF_val)
		else:
			thresholdMFArray.InsertNextValue(-1)
		meanMFValue = meanMFValue/len(pointsPool)
		maxMFArray.InsertNextValue(maxMFValue)
		meanMFArray.InsertNextValue(meanMFValue)

	# Add Dach1 Array to point data
	maxMFArray.SetName('Max ' + factors[molecularFactor] + ' (' + str(INTERPOLATING_RADIUS) + ')') 
	mesh.GetPointData().AddArray(maxMFArray)
	meanMFArray.SetName('Mean ' + factors[molecularFactor] + ' (' + str(INTERPOLATING_RADIUS)+ ')') 
	mesh.GetPointData().AddArray(meanMFArray)
	thresholdMFArray.SetName('Threshold ' + factors[molecularFactor] + ' (' + str(CX40THRESHOLD)+ ')') 
	mesh.GetPointData().AddArray(thresholdMFArray)

# Convert TA_WSS and Pressure to micron units and new arrays to point data
convertData(mesh,'vTAWSS','vTAWSS (dynes/cm^2)',WSS_CONV_FACTOR)
convertData(mesh,'pressure_avg','pressure_avg (mmHg)',1/P_CONV_FACTOR)


# Write a new .vtp file that has all converted values and mapped values.
w = vtk.vtkXMLPolyDataWriter()
w.SetInputData(mesh)
w.SetFileName('all_results_mapped_interpolated_full_threshold.vtp')
w.Write()
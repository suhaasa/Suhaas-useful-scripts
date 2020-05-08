# Bivariate Spatial Association
import sys
import vtk
import numpy as np

def buildConnectivityMatrix(model):

	numPts = model.GetNumberOfPoints()
	print(numPts)
	c = np.zeros((numPts,numPts))
	for i in xrange(0,numPts):

		cell_list = vtk.vtkIdList()
		cell_list.NewInstance()
		
		model.GetPointCells(i,cell_list)
		
		for j in xrange(0,cell_list.GetNumberOfIDs()):
			pt_list = vtk.vtkIdList()
			model.GetPointID(j,pt_list)
			for k in xrange(0,pt_list.GetNumberOfPoints()):
				if (c[i][pt_list.GetId(j)] == 0):
					c[i][pt_list.GetId(j)] = 1
	return c

def getConnectedVerticesAll(model, seedPt):
	cell_list = vtk.vtkIdList()
	connectedPts_list = vtk.vtkIdList()
	model.GetPointCells(seedPt,cell_list)
	for j in xrange(0,cell_list.GetNumberOfIds()):
		pt_list = vtk.vtkIdList()
		pt_list = model.GetCell(cell_list.GetId(j)).GetPointIds()
		for k in xrange(0,pt_list.GetNumberOfIds()):
			connectedPts_list.InsertUniqueId(pt_list.GetId(k))
	return connectedPts_list

def getConnectedVerticesNotIncludingSeed(model, seedPt):
	cell_list = vtk.vtkIdList()
	connectedPts_list = vtk.vtkIdList()
	model.GetPointCells(seedPt,cell_list)
	for j in xrange(0,cell_list.GetNumberOfIds()):
		pt_list = vtk.vtkIdList()
		pt_list = model.GetCell(cell_list.GetId(j)).GetPointIds()
		for k in xrange(0,pt_list.GetNumberOfIds()):
			if (pt_list.GetId(k) != seedPt):
				connectedPts_list.InsertUniqueId(pt_list.GetId(k))
	return connectedPts_list

def calcDistance2Points(model, pt1,pt2):
	x1,y1,z1 = model.GetPoint(pt1)
	x2,y2,z2 = model.GetPoint(pt2)
	distance = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**(.5)
	return distance

def maxDistanceBetweenPoints(model, seedPt, connectedPts_list):
	max = 0
	for i in xrange(0,connectedPts_list.GetNumberOfIds()):
		distance = calcDistance2Points(model, seedPt,connectedPts_list.GetId(i))
		if(distance > max):
			max = distance
	return max

def getConnectedVerticesWithinRadius(model, seedPt, radius):
	cell_list = vtk.vtkIdList()
	connectedPts_list = vtk.vtkIdList()
	model.GetPointCells(seedPt,cell_list)
	radiusReached = 0
	while (radiusReached == 0):
		for j in xrange(0,cell_list.GetNumberOfIds()):
			pt_list = vtk.vtkIdList()
			pt_list = model.GetCell(cell_list.GetId(j)).GetPointIds()
			for k in xrange(0,pt_list.GetNumberOfIds()):
				connectedPts_list.InsertUniqueId(pt_list.GetId(k))
		if (maxDistanceBetweenPoints(model, seedPt,connectedPts_list) > radius):
			radiusReached = 1
		else:
			connectedCell_list = vtk.vtkIdList()
			for i in xrange(0,pt_list.GetNumberOfIds()):
				model.GetPointCells(pt_list.GetId(i),connectedCell_list)
				for j in xrange(0,connectedCell_list.GetNumberOfIds()):
					cell_list.InsertUniqueId(connectedCell_list.GetId(j))


	return connectedPts_list

def calcMagnitude(vector):
	mag = 0
	for i in xrange(0,len(vector)):
		mag += vector[i]**2
	return mag ** (.5)

def calcZ(value,mean,std):

	return (value-mean) / std

def calcMean(model,factor):
	sum = 0
	count = 0
	numPts = model.GetNumberOfPoints()
	for i in xrange(0,numPts):
		vector = model.GetPointData().GetArray(factor).GetTuple(i)
		mag = calcMagnitude(vector)
		if (mag > 0 and mag < 200):
			sum += mag
			count = count + 1
	return sum/count

def calcStdev(model, factor):
	squaredsum = 0
	mean = calcMean(model, factor)
	numPts = model.GetNumberOfPoints()
	count = 0
	for i in xrange(0,numPts):
		vector = model.GetPointData().GetArray(factor).GetTuple(i)
		mag = calcMagnitude(vector)
		if (mag > 0 and mag < 200):
			squaredsum += (mag-mean) ** 2
			count = count + 1

	return (squaredsum/(count-1)) ** (.5)

def getValue(model, ptID, factor):
	return calcMagnitude(model.GetPointData().GetArray(factor).GetTuple(ptID))

vtk_file_prefix = 'all_results_mapped_interpolated_threshold_correlated'
Mfactors = ['Threshold Dach1 (20)']
Hfactors = ['vTAWSS (dynes/cm^2)']
RADIUS = 5
# First, read in the .vtp file containing your quantities of interest
datareader=vtk.vtkXMLPolyDataReader()
datareader.SetFileName(vtk_file_prefix + '.vtp')
datareader.Update()
# Read your data into another polydata variable for manipulation
model=vtk.vtkPolyData()
model=datareader.GetOutput()
print('Start')
#c = buildConnectivityMatrix(model)
numPts = model.GetNumberOfPoints()
for MF in xrange(0,len(Mfactors)):
	meanMF = calcMean(model,Mfactors[MF])
	stdMF = calcStdev(model,Mfactors[MF])
	print('mean ' + Mfactors[MF] + ': ' + str(meanMF))
	print('std ' + Mfactors[MF] + ': ' + str(stdMF))
	for HF in xrange(0,len(Hfactors)):
		meanHF = calcMean(model,Hfactors[HF])
		stdHF = calcStdev(model,Hfactors[HF])
		print('mean ' + Hfactors[HF] + ': ' + str(meanHF))
		print('std ' + Hfactors[HF] + ': ' + str(stdHF))
		R_MF_Array = vtk.vtkDoubleArray()
		R_HF_Array = vtk.vtkDoubleArray()
		R_MF_HF_Array = vtk.vtkDoubleArray()
		R_MF0_HF0_Array = vtk.vtkDoubleArray()
		for seedPtID in xrange(0,numPts):
			if(getValue(model,seedPtID,Mfactors[MF])>0):
				MF0_Z = calcZ(getValue(model,seedPtID,Mfactors[MF]),meanMF,stdMF)
				HF0_Z = calcZ(getValue(model,seedPtID,Hfactors[HF]),meanHF,stdHF)
				connectedPts_list = getConnectedVerticesWithinRadius(model,seedPtID,RADIUS)
				print(connectedPts_list)
				weighted_MF_Z = 0
				MF_neighbors = 0
				weighted_HF_Z = 0
				HF_neighbors = 0
				for connectedPt_index in xrange(0,connectedPts_list.GetNumberOfIds()):
					MF_val = getValue(model,connectedPts_list.GetId(connectedPt_index),Mfactors[MF])
					if (MF_val>0):
						weighted_MF_Z += calcZ(MF_val,meanMF,stdMF)
						MF_neighbors = MF_neighbors + 1
					HF_val = getValue(model,connectedPts_list.GetId(connectedPt_index),Hfactors[HF])
					if (HF_val>0):
						weighted_HF_Z += calcZ(HF_val,meanHF,stdHF)
						HF_neighbors = HF_neighbors + 1
				if (MF_neighbors > 0 and HF_neighbors > 0):
					R_MF_Array.InsertNextValue(MF0_Z*(weighted_HF_Z/HF_neighbors))
					R_HF_Array.InsertNextValue(HF0_Z*(weighted_MF_Z/MF_neighbors))
					R_MF_HF_Array.InsertNextValue((weighted_MF_Z/MF_neighbors)*(weighted_HF_Z/HF_neighbors))
					R_MF0_HF0_Array.InsertNextValue(MF0_Z*HF0_Z)
				else:
					R_MF_Array.InsertNextValue(0)
					R_HF_Array.InsertNextValue(0)
					R_MF_HF_Array.InsertNextValue(0)
					R_MF0_HF0_Array.InsertNextValue(0)
			else:
				R_MF_Array.InsertNextValue(0)
				R_HF_Array.InsertNextValue(0)
				R_MF_HF_Array.InsertNextValue(0)
				R_MF0_HF0_Array.InsertNextValue(0)
		R_MF_Array.SetName('R ' + Mfactors[MF] + '(' + Hfactors[HF] + ')') 
		model.GetPointData().AddArray(R_MF_Array)
		R_HF_Array.SetName('R ' + Hfactors[HF] + '(' + Mfactors[MF] + ')')
		model.GetPointData().AddArray(R_HF_Array)
		R_MF_HF_Array.SetName('R ' + Hfactors[HF] + '-' + Mfactors[MF] )
		model.GetPointData().AddArray(R_MF_HF_Array)
		R_MF0_HF0_Array.SetName('R0 ' + Hfactors[HF] + '-' + Mfactors[MF] )
		model.GetPointData().AddArray(R_MF0_HF0_Array)

w = vtk.vtkXMLPolyDataWriter()
w.SetInputData(model)
w.SetFileName('all_results_mapped_interpolated_threshold_cross_correlation.vtp')
w.Write()





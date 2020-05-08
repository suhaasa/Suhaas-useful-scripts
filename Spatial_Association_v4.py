# Spatial Association

import sys
import vtk
import numpy as np 
import random as rand
from scipy import stats
import scipy.stats as st
import scipy.stats.mstats as zsc

def determineDiscriminingStastics(HF_vector,MF_vector):
	sum = 0
	for i in xrange(0, len(HF_vector)):
		sum = sum + abs(HF_vector[i]-MF_vector[i])
	return sum / len(HF_vector)

def calcSpatialSmoothingScalar(vector, SL_vector, mean):
	Var = 0
	SL_Var = 0
	for i in xrange(0, len(vector)):
		if(vector[i] != 0 and SL_vector[i] != 0):
			Var += (vector[i]-mean)**2
			SL_Var += (SL_vector[i]-mean)**2
	return (SL_Var/Var)

def calcPearsonsCoeff(X, Y, X_mean, X_stdev, Y_mean, Y_stdev):
	N = len(X)
	XY = 0
	for i in xrange(0,N):
		XY += X[i] * Y[i]
	X_2 = 0
	for i in xrange(0,N):
		X_2 += X[i]**2
	Y_2 = 0
	for i in xrange(0,N):
		Y_2 += Y[i]**2
	Xsum = 0
	for i in xrange(0,N):
		Xsum += X[i]
	Ysum = 0
	for i in xrange(0,N):
		Ysum += Y[i]
	num = N*XY - Xsum * Ysum
	dem = ((N*X_2 - Xsum**2)*(N*Y_2 - Ysum**2))**(.5)
	r = num/dem
		# print(X[i])
		# print(X_mean)
		# print(X_stdev)
	return r

def determineSpatialAssociationStatistic(MF_vector,HF_vector,SL_MF_vector,SL_HF_vector, \
	MF_mean,MF_stdev,HF_mean,HF_stdev):
	HF_SSS = calcSpatialSmoothingScalar(HF_vector,SL_HF_vector,HF_mean)
	MF_SSS = calcSpatialSmoothingScalar(MF_vector,SL_MF_vector,MF_mean)
	R_MF_HF = calcPearsonsCoeff(SL_MF_vector,SL_HF_vector,MF_mean,MF_stdev,HF_mean,HF_stdev)
	# print('HF_SSS=' + str(HF_SSS))
	# print('MF_SSS=' + str(MF_SSS))
	# print('R=' + str(R_MF_HF))
	return HF_SSS**(.5) * MF_SSS**(.5) * R_MF_HF

def getSamplePts(numSamplePts, validPts):
	return rand.sample(validPts, numSamplePts)

def getValidPtIDs(model, MF_quantity, MF_threshold, WSSG_threshold, option):
	numPts=model.GetNumberOfPoints()
	validIDs = []
	for ID in xrange(0,numPts):
		MF_val = model.GetPointData().GetArray(MF_quantity).GetValue(ID)
		Gradients = getValue(model,ID, 'Gradients')
		if(MF_val>MF_threshold):
			if(option == 1 and Gradients > WSSG_threshold):
				validIDs.append(ID)
			elif(option == 3 and Gradients < WSSG_threshold):
				validIDs.append(ID)
			elif(option == 2):
				validIDs.append(ID)
	return validIDs

def getSurrogateData(HF_data):
	HF_shuffled_data = np.zeros((len(HF_data),len(HF_data[0])))
	for i in xrange(0,len(HF_data)):
		HF_vector = np.transpose(HF_data[i])
		HF_shuffled = np.random.permutation(HF_vector)
		for j in xrange(0,len(HF_data[0])):
			HF_shuffled_data[i][j] = HF_shuffled[j]
	return HF_shuffled_data

def getSurrogateData2(HF_data, SL_HF_data):
	HF_shuffled_data = np.zeros((2,len(HF_data),len(HF_data[0])))
	for i in xrange(0,len(HF_data)):
		HF_vector = np.transpose(HF_data[i])
		SL_HF_vector = np.transpose(SL_HF_data[i])
		indicies = rand.sample(np.arange(len(HF_data[0])),len(HF_data[0]))
		for j in xrange(0,len(HF_data[0])):
			HF_shuffled_data[0][i][j] = HF_vector[indicies[j]]
			HF_shuffled_data[1][i][j] = SL_HF_vector[indicies[j]]

	return HF_shuffled_data

def randomizeData(data, indicies):
	shuffled_data = np.zeros((len(data),len(data[0])))
	for i in xrange(0,len(data)):
		for j in xrange(0,len(data[0])):
			shuffled_data[i][j] = data[i][indicies[j]]
	return shuffled_data

def randomizePointData(model, Mfactors, Hfactors):
	numPts = model.GetNumberOfPoints()
	randomIndicies = rand.sample(np.arange(numPts),numPts)
	random_MF_PointDataArray_list = []
	random_HF_PointDataArray_list = []
	for MF in xrange(0,len(Mfactors)):
		random_MF_PointDataArray_list.append(vtk.vtkDoubleArray())
	for HF in xrange(0,len(Hfactors)):
		random_HF_PointDataArray_list.append(vtk.vtkDoubleArray())

	for i in xrange(0,len(randomIndicies)):
		for MF in xrange(0,len(Mfactors)):
			random_MF_PointDataArray_list[MF].InsertNextValue(getValue(model,randomIndicies[i],Mfactors[MF]))
			for HF in xrange(0,len(Hfactors)):
				random_HF_PointDataArray_list[HF].InsertNextValue(getValue(model,randomIndicies[i],Hfactors[HF]))
	for MF in xrange(0,len(Mfactors)):
		if(model.GetPointData().HasArray('Randomized ' + Mfactors[MF]) == 1):
			model.GetPointData().RemoveArray('Randomized ' + Mfactors[MF])
		random_MF_PointDataArray_list[MF].SetName('Randomized ' + Mfactors[MF])
		model.GetPointData().AddArray(random_MF_PointDataArray_list[MF])

	for HF in xrange(0,len(Hfactors)):
		if(model.GetPointData().HasArray('Randomized ' + Hfactors[HF]) == 1):
			model.GetPointData().RemoveArray('Randomized ' + Hfactors[HF])
		random_HF_PointDataArray_list[HF].SetName('Randomized ' + Hfactors[HF])
		model.GetPointData().AddArray(random_HF_PointDataArray_list[HF])



def getData(model,IDs,quantities):
	data = np.zeros((len(quantities),len(IDs)))
	for ptID in xrange(0,len(IDs)):
		for iq in xrange(0,len(quantities)):
			data[iq][ptID] = getValue(model,IDs[ptID],quantities[iq])
	return data

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

def addSpatialLag(model,factors,PtIDs):
	ConnectionCounts = np.zeros(len(factors))
	numPts = model.GetNumberOfPoints()
	for f in xrange(0,len(factors)):
		SL = vtk.vtkDoubleArray()
		for PtID in xrange(0,len(PtIDs)):
			connectedPts_list = getConnectedVerticesWithinRadius(model,PtIDs[PtID],5)
			sum = 0
			count = 0
			for i in xrange(0,connectedPts_list.GetNumberOfIds()):
				val = getValue(model,connectedPts_list.GetId(i),factors[f])
				if (val > 0):
					ConnectionCounts[f] += 1
					sum += val
					count += 1
			if(count != 0):
				SL.InsertNextValue(sum/count)
			else:
				SL.InsertNextValue(0)
		ConnectionCounts[f] /= numPts
		SL.SetName('SL ' + factors[f])
		if(model.GetPointData().HasArray('SL ' + factors[f]) == 1):
			model.GetPointData().RemoveArray('SL ' + factors[f])
		model.GetPointData().AddArray(SL)
	return ConnectionCounts

def calcZ(value,mean,std):
	return (value-mean) / std

def calcMean(vector):
	sum = 0
	count = 0
	for i in xrange(0,len(vector)):
		sum += vector[i]
		count = count + 1
	return sum/count

def calcStdev(vector):
	squaredsum = 0
	mean = calcMean(vector)
	count = 0
	for i in xrange(0,numPts):
		squaredsum += (mag-mean) ** 2
		count = count + 1
	return (squaredsum/(count-1)) ** (.5)

def getValue(model, ptID, factor):
	return calcMagnitude(model.GetPointData().GetArray(factor).GetTuple(ptID))

def calcMagnitude(vector):
	mag = 0
	for i in xrange(0,len(vector)):
		mag += vector[i]**2
	return mag ** (.5)

def SLconversion(factors):
	SLfactors = []
	for i in  xrange(0,len(factors)):
		SLfactor = 'SL ' + factors[i]
		SLfactors.append(SLfactor)
	return SLfactors

def Randomizedconversion(factors):
	SLfactors = []
	for i in  xrange(0,len(factors)):
		SLfactor = 'Randomized ' + factors[i]
		SLfactors.append(SLfactor)
	return SLfactors

def parseOption(num):
	if(num == 1):
		return 'high'
	elif(num == 2):
		return 'all'
	return 'low'

vtk_file_prefix = 'all_results_mapped_interpolated_threshold_correlated'
Mfactors = ['Threshold Dach1 (20)']
Hfactors = ['vTAWSS (dynes/cm^2)','Gradients','pressure_avg (mmHg)']
SL_Mfactors = SLconversion(Mfactors)
SL_Hfactors = SLconversion(Hfactors)
Rand_Mfactors = Randomizedconversion(Mfactors)
Rand_Hfactors = Randomizedconversion(Hfactors)
SL_Rand_Mfactors = SLconversion(Rand_Mfactors)
SL_Rand_Hfactors = SLconversion(Rand_Hfactors)
MF_THRESHOLD = 20
WSSG_THRESHOLD = 1
# Above - 1, all - 2, below - 3
AboveAllBelow = 3
# First, read in the .vtp file containing your quantities of interest
datareader=vtk.vtkXMLPolyDataReader()
datareader.SetFileName(vtk_file_prefix + '.vtp')
datareader.Update()

# Read your data into another polydata variable for manipulation
model=vtk.vtkPolyData()
model=datareader.GetOutput()

numPts=model.GetNumberOfPoints()

numSamples = 1
numSurrogates = 1

All_stats_sample = np.zeros((numSamples,len(Mfactors),len(Hfactors)))
validPts = getValidPtIDs(model,Mfactors[0],MF_THRESHOLD,WSSG_THRESHOLD,AboveAllBelow)
numSamplePts = len(validPts)/1
All_stats_collection = np.zeros((numSamples, numSurrogates, len(Mfactors),len(Hfactors)))

MF_connectivity = addSpatialLag(model,Mfactors,np.arange(numPts))
HF_connectivity = addSpatialLag(model,Hfactors,np.arange(numPts))
print('Start')
for i_Smp in xrange(0,numSamples):
	SamplePtIDs = getSamplePts(numSamplePts,validPts)
	# Get the data of all factors for the sample selected
	HF_data_original = getData(model,SamplePtIDs,Hfactors)
	MF_data_original = getData(model,SamplePtIDs,Mfactors)
	SL_HF_data_original = getData(model,SamplePtIDs,SL_Hfactors)
	SL_MF_data_original = getData(model,SamplePtIDs,SL_Mfactors)
	# Caclulate means and stdevs for factors
	MF_means = np.mean(MF_data_original, axis=1, dtype=np.float64)
	MF_stdev = np.std(MF_data_original, axis=1, dtype=np.float64)
	HF_means = np.mean(HF_data_original, axis=1, dtype=np.float64)
	HF_stdev = np.std(HF_data_original, axis=1, dtype=np.float64)

	Stats_original = np.zeros((len(Mfactors),len(Hfactors)))
	All_stats_surrogate = np.zeros((numSurrogates+1,len(Mfactors),len(Hfactors)))
	HF_SSS = np.zeros((numSurrogates+1,len(Hfactors)))
	MF_SSS = np.zeros((numSurrogates+1,len(Mfactors)))
	# Determine the statistic for the original data
	for MF in xrange(0,len(Mfactors)):
		for HF in xrange(0,len(Hfactors)):
			Stats_original[MF][HF] = determineSpatialAssociationStatistic( \
				MF_data_original[MF],HF_data_original[HF],SL_MF_data_original[MF],SL_HF_data_original[HF], \
				MF_means[MF],MF_stdev[MF],HF_means[HF],HF_stdev[HF])
			All_stats_surrogate[0][MF][HF] = Stats_original[MF][HF]
			HF_SSS[0][MF] = calcSpatialSmoothingScalar(HF_data_original[HF],SL_HF_data_original[HF],HF_means[HF])
			MF_SSS[0][MF] = calcSpatialSmoothingScalar(MF_data_original[MF],SL_MF_data_original[MF],MF_means[MF])
	print(Stats_original)
	# Determine the statistic for the surrogate data
	for i_Sur in xrange(1,numSurrogates+1):
		randomIndicies = rand.sample(np.arange(len(SL_MF_data_original[0])),len(SL_MF_data_original[0]))
		#print('Random Indicies: ')
		#print(randomIndicies)
		# SL_MF_data_surrogate = randomizeData(SL_MF_data_original,randomIndicies)
		# SL_HF_data_surrogate = randomizeData(SL_HF_data_original,randomIndicies)
		randomizePointData(model, Mfactors, Hfactors)
		addSpatialLag(model, Rand_Mfactors, np.arange(numPts))
		addSpatialLag(model, Rand_Hfactors, np.arange(numPts))
		HF_data_surrogate = getData(model,SamplePtIDs,Rand_Hfactors)
		MF_data_surrogate = getData(model,SamplePtIDs,Rand_Mfactors)
		SL_HF_data_surrogate = getData(model,SamplePtIDs,SL_Rand_Hfactors)
		SL_MF_data_surrogate = getData(model,SamplePtIDs,SL_Rand_Mfactors)
		# all_HF_data_surrogate = getSurrogateData2(HF_data_original,SL_HF_data_original)
		# HF_data_surrogate = all_HF_data_surrogate[0]
		# SL_HF_data_surrogate = all_HF_data_surrogate[1]
		for MF in xrange(0,len(Mfactors)):
			for HF in xrange(0,len(Hfactors)):
				All_stats_surrogate[i_Sur][MF][HF] = determineSpatialAssociationStatistic( \
					MF_data_surrogate[MF],HF_data_surrogate[HF],SL_MF_data_surrogate[MF],SL_HF_data_surrogate[HF], \
					MF_means[MF],MF_stdev[MF],HF_means[HF],HF_stdev[HF])
				HF_SSS[i_Sur][HF] = calcSpatialSmoothingScalar(HF_data_surrogate[HF],SL_HF_data_surrogate[HF],HF_means[HF])
				MF_SSS[i_Sur][MF] = calcSpatialSmoothingScalar(MF_data_surrogate[MF],SL_MF_data_surrogate[MF],MF_means[MF])
		print(i_Sur)

	mean = np.mean(All_stats_surrogate, axis=0)
	stdev = np.std(All_stats_surrogate, axis=0)
	print(mean)
	print(stdev)
	for MF in xrange(0,len(Mfactors)):
		for HF in xrange(0,len(Hfactors)):
			zscore = (All_stats_surrogate[0][MF][HF] - mean[MF][HF]) / stdev[MF][HF]
			All_stats_sample[i_Smp][MF][HF] = zscore

	for i_Sur in xrange(0,numSurrogates):
		for MF in xrange(0,len(Mfactors)):
			for HF in xrange(0,len(Hfactors)):
				All_stats_collection[i_Smp][i_Sur][MF][HF] = All_stats_surrogate[i_Sur][MF][HF]
	print(i_Smp)

SSS_max = [np.amax(MF_SSS),np.amax(HF_SSS)]
SSS_mean = [np.mean(MF_SSS),np.mean(HF_SSS)]
print(All_stats_sample)
mean_Z = np.mean(All_stats_sample, axis=0)
print(mean_Z)

# Make data file
option = parseOption(AboveAllBelow)
outfile_name = 'Spatial_Association_data_' + str(option) + '_WSSG_' + str(WSSG_THRESHOLD*10)
outfile = open(outfile_name + '.csv','w')
out_string = ''

for i_Sur in xrange(0,numSurrogates):
	for MF in xrange(0,len(Mfactors)):
		for HF in xrange(0,len(Hfactors)):
			for i_Smp in xrange(0,numSamples):
				out_string += str(All_stats_collection[i_Smp][i_Sur][MF][HF]) + ', '
	out_string += '\n'

outfile.write(out_string)
outfile.close()

# Make summary file
summary_outfile_name = 'Summary_' + outfile_name
summary_outfile = open(outfile_name + '.txt','w')
summary_out_string  = 'Threshold (WSSG): ' + str(WSSG_THRESHOLD) + '\n'
summary_out_string += 'Threshold (MF): ' + str(MF_THRESHOLD) + '\n'
summary_out_string += 'Fraction of total area: ' + str(float(len(validPts))/numPts) + '\n'
summary_out_string += 'Average SSS: ' + str(SSS_mean) + '\n'
summary_out_string += 'Max SSS: ' + str(SSS_max) + '\n'
summary_out_string += 'MF Means: ' + str(MF_means) + '\n'
summary_out_string += 'MF Standard Deviations: ' + str(MF_stdev) + '\n'
summary_out_string += 'HF Means: ' + str(HF_means) + '\n'
summary_out_string += 'HF Standard Deviations: ' + str(HF_stdev) + '\n'
summary_out_string += 'Connectivity: ' + str(MF_connectivity) + str(HF_connectivity) + '\n'
summary_out_string += 'Zscore: ' + str(mean_Z) + '\n'

summary_outfile.write(summary_out_string)
summary_outfile.close()

w = vtk.vtkXMLPolyDataWriter()
w.SetInputData(model)
w.SetFileName('all_results_mapped_interpolated_threshold__spatial_assocation.vtp')
w.Write()
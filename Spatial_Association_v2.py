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
		Var += (vector[i]-mean)**2
		SL_Var += (SL_vector[i]-mean)**2
	return (SL_Var/Var)

def calcPearsonsCoeff(X, Y, mean, stdev):
	r = 0
	for i in xrange(0,len(X)):
		r += calcZ(X[i],mean,stdev) * calcZ(Y[i],mean,std)
	return r

def determineSpatialAssociationStatistic(HF_vector,MF_vector,SL_HF_vector,SL_MF_vector,mean,stdev):
	HF_SSS = calcSpatialSmoothingScalar(HF_vector,SL_HF_vector,mean)
	MF_SSS = calcSpatialSmoothingScalar(MF_vector,SL_MF_vector,mean)
	R_MF_HF = calcPearsonsCoeff(SL_MF_vector,SL_HF_vector)
	return HF_SSS**(.5) * MF_SSS**(.5) * R_MF_HF

def getSamplePts(numSamplePts, validPts):
	return rand.sample(validPts, numSamplePts)

def getValidPtIDs(model, MF_quantity, threshold):
	numPts=model.GetNumberOfPoints()
	validIDs = []
	for ID in xrange(0,numPts):
		MF_val = model.GetPointData().GetArray(MF_quantity).GetValue(ID)
		Gradients = getValue(model,'Gradients',ID)
		if(MF_val>threshold and Gradients > 2):
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

def getData(model,IDs,quantities):
	data = np.zeros((len(quantities),len(IDs)))
	for ptID in xrange(0,len(IDs)):
		for iq in xrange(0,len(quantities)):
			data[iq][ptID] = getValue(model,quantities[iq],IDs[ptID])
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

def addSpatialLag(model,factors):
	numPts = model.GetNumberOfPoints()
	for f in xrange(0,len(factors))
		SL = vtk.vtkDoubleArray()

		for PtID in xrange(0,numPts):
			connectedPts_list = getConnectedVerticesNotIncludingSeed(model,PtID)
			sum = 0
			count = 0
			for i in xrange(0,connectedPts_list.GetNumberOfIds()):
				val = getValue(model,connectedPts_list.GetId(i),factors[f])
				if (val > 0):
					sum += val
					count += 1
			SL.InsertNextValue(sum/count)
		SL.SetName('SL ' + factors[f])
		model.GetPointData().AddArray(SL)

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

def calcMagnitude(vector):
	mag = 0
	for i in xrange(0,len(vector)):
		mag += vector[i]**2
	return mag ** (.5)

vtk_file_prefix = 'all_results_mapped_interpolated_threshold_correlated'
Mfactors = ['Threshold Dach1 (20)']
Hfactors = ['vTAWSS (dynes/cm^2)']

# First, read in the .vtp file containing your quantities of interest
datareader=vtk.vtkXMLPolyDataReader()
datareader.SetFileName(vtk_file_prefix + '.vtp')
datareader.Update()

# Read your data into another polydata variable for manipulation
model=vtk.vtkPolyData()
model=datareader.GetOutput()

numPts=model.GetNumberOfPoints()

numSamples = 10
numSurrogates = 10000

All_stats_sample = np.zeros((numSamples,len(Mfactors),len(Hfactors)))
validPts = getValidPtIDs(model,Mfactors[0],0)
numSamplePts = len(validPts)/2
All_stats_collection = np.zeros((numSamples, numSurrogates, len(Mfactors),len(Hfactors)))
for i_Smp in xrange(0,numSamples):
	SamplePtIDs = getSamplePts(numSamplePts,validPts)
	HF_data_original = getData(model,SamplePtIDs,Hfactors)
	MF_data_original = getData(model,SamplePtIDs,Mfactors)
	Stats_original = np.zeros((len(Mfactors),len(Hfactors)))
	All_stats_surrogate = np.zeros((numSurrogates+1,len(Mfactors),len(Hfactors)))
	for MF in xrange(0,len(Mfactors)):
		for HF in xrange(0,len(Hfactors)):
			Stats_original[MF][HF] = determineDiscriminingStastics(HF_data_original[HF],MF_data_original[MF])
			All_stats_surrogate[0][MF][HF] = Stats_original[MF][HF]

	for i_Sur in xrange(1,numSurrogates):
		HF_data_surrogate = getSurrogateData(HF_data_original)
		for MF in xrange(0,len(Mfactors)):
			for HF in xrange(0,len(Hfactors)):
				All_stats_surrogate[i_Sur][MF][HF] = determineDiscriminingStastics(HF_data_surrogate[HF],MF_data_original[MF])

	mean = np.mean(All_stats_surrogate, axis=0)
	stdev = np.std(All_stats_surrogate, axis=0)
	for MF in xrange(0,len(Mfactors)):
		for HF in xrange(0,len(Hfactors)):
			zscore = (All_stats_surrogate[0][MF][HF] - mean[MF][HF]) / stdev[MF][HF]
			All_stats_sample[i_Smp][MF][HF] = zscore

	for i_Sur in xrange(0,numSurrogates):
		for MF in xrange(0,len(Mfactors)):
			for HF in xrange(0,len(Hfactors)):
				All_stats_collection[i_Smp][i_Sur][MF][HF] = All_stats_surrogate[i_Sur][MF][HF]
	print(i_Smp)

print(All_stats_sample)
mean = np.mean(All_stats_sample, axis=0)
print(mean)
outfile_name = 'Spatial_Association_data.csv'
outfile = open(outfile_name,'w')
out_string = ''

for i_Sur in xrange(0,numSurrogates):
	for MF in xrange(0,len(Mfactors)):
		for HF in xrange(0,len(Hfactors)):
			for i_Smp in xrange(0,numSamples):
				out_string += str(All_stats_collection[i_Smp][i_Sur][MF][HF]) + ', '
	out_string += '\n'

outfile.write(out_string)
outfile.close()


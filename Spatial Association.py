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

def getSamplePts(numSamplePts, validPts):
	return rand.sample(validPts, numSamplePts)

def getValidPtIDs(model, MF_quantity, threshold):
	numPts=model.GetNumberOfPoints()
	validIDs = []
	for ID in xrange(0,numPts):
		MF_val = model.GetPointData().GetArray(MF_quantity).GetValue(ID)
		Gradients = model.GetPointData().GetArray('Gradients').GetValue(ID)
		if(MF_val>threshold):
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
			vector = model.GetPointData().GetArray(quantities[iq]).GetTuple(IDs[ptID])
			mag = 0
			for i in xrange(0,len(vector)):
				mag = mag + vector[i]**2
			data[iq][ptID] = mag**(.5)
	return data


vtk_file_prefix = 'all_results_mapped_interpolated_threshold_correlated'
Mfactors = ['Threshold Dach1 (20)']
Hfactors = ['vTAWSS (dynes/cm^2)','Gradients','pressure_avg (mmHg)']

# First, read in the .vtp file containing your quantities of interest
datareader=vtk.vtkXMLPolyDataReader()
datareader.SetFileName(vtk_file_prefix + '.vtp')
datareader.Update()

# Read your data into another polydata variable for manipulation
model=vtk.vtkPolyData()
model=datareader.GetOutput()

numPts=model.GetNumberOfPoints()

numSamples = 10
numSurrogates = 500

All_stats_sample = np.zeros((numSamples,len(Mfactors),len(Hfactors)))
validPts = getValidPtIDs(model,Mfactors[0],0)
numSamplePts = len(validPts)/20
All_stats_surrogate_collection =  np.zeros((numSamples,numSurrogates,len(Mfactors),len(Hfactors)))
for i_Smp in xrange(0,numSamples):
	SamplePtIDs = getSamplePts(numSamplePts,validPts)
	HF_data_original = getData(model,SamplePtIDs,Hfactors)
	MF_data_original = getData(model,SamplePtIDs,Mfactors)
	Stats_original = np.zeros((len(Mfactors),len(Hfactors)))
	All_stats_surrogate = np.zeros((numSurrogates,len(Mfactors),len(Hfactors)))
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
	for MF in xrange(0,len(Mfactors)):
		for HF in xrange(0,len(Hfactors)):
			for i_Sur in xrange(0,numSurrogates):
				All_stats_surrogate_collection[i_Smp][i_Sur][MF][HF] = All_stats_surrogate[i_Sur][MF][HF]
	print(i_Smp)

print(All_stats_sample)
mean = np.mean(All_stats_sample, axis=0)
print(mean)

output_filename = 'Spatial Association Statistics'
outfile = open(output_filename, 'w')

out_string = ''
for i in xrange(0, numSurrogates):
	for MF in xrange(0, len(Mfactors)):
		for HF in xrange(0, len(Hfactors)):
			for j in xrange(0, numSamples):
				out_string += str(All_stats_surrogate_collection[j][i][MF][HF]) + ', '
		out_string += ', '
	out_string += '\n'
outfile.write(out_string)
outfile.close()






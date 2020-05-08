import sys
import vtk
import numpy as np
from vtk.util import numpy_support as VN
# Model surface
vtk_file_prefix = 'all_results_mapped_interpolated_threshold_gradients'

# Molecular factors
Mfactors = ['Dach1']
MF_THRESHOLD_FACTOR = [.1,.15,.2]
# Hemodynamic factors
Hfactors = ['vTAWSS (dynes/cm^2)','Gradients','pressure_avg (mmHg)']
HF_THRESHOLD_FACTOR = [.04,.002,.8]
# First, read in the .vtp file containing your quantities of interest
datareader=vtk.vtkXMLPolyDataReader()
datareader.SetFileName(vtk_file_prefix + '.vtp')
datareader.Update()
# Read your data into another polydata variable for manipulation
model=vtk.vtkPolyData()
model=datareader.GetOutput()

numPts=model.GetNumberOfPoints()

for TF in xrange(0,len(MF_THRESHOLD_FACTOR)):
	# initalize MF storage arrays
	MF_STDEV = np.zeros(len(Mfactors))
	MF_VAR = np.zeros(len(Mfactors))
	MF_SDSQ = np.zeros(len(Mfactors))
	MF_AV = np.zeros(len(Mfactors))
	MF_tot = np.zeros(len(Mfactors))
	MF_count = np.zeros(len(Mfactors))
	MF_tot_2 = np.zeros(len(Mfactors))
	MF_max = np.zeros(len(Mfactors))

	for MF in xrange(0,len(Mfactors)):
		output_filename = 'Point_correlation_' + Mfactors[MF] + '.dat'
		

		HF_tot = np.zeros(len(Hfactors))
		HF_tot_2 = np.zeros(len(Hfactors))
		HF_count = np.zeros(len(Hfactors))
		HF_max = np.zeros(len(Hfactors))
		MF_HF = np.zeros(len(Hfactors))
		
		
		# Calculates totals and counts
		for pointID in xrange(0,numPts):
			MF_val = abs(model.GetPointData().GetArray('Threshold ' + Mfactors[MF] + ' (20)').GetValue(pointID))
			if(MF_val>0):
				MF_tot[MF] = MF_tot[MF] + MF_val
				MF_tot_2[MF] = MF_tot_2[MF] + MF_val ** 2
				MF_count[MF] = MF_count[MF] + 1
				if(MF_val>MF_max[MF]):
					MF_max[MF] = MF_val
				for HF in xrange(0,len(Hfactors)):
					HF_vector = model.GetPointData().GetArray(Hfactors[HF]).GetTuple(pointID)
					HF_val = 0
					for i in xrange(0,len(HF_vector)):
						HF_val = HF_val + HF_vector[i]**2
					HF_val = HF_val**(.5)
					HF_tot[HF] = HF_tot[HF] + HF_val
					HF_tot_2[HF] = HF_tot_2[HF] + HF_val ** 2
					HF_count[HF] = HF_count[HF] + 1
					MF_HF[HF] = MF_HF[HF] + (MF_val * HF_val)
					if(HF_val>HF_max[HF]):
						HF_max[HF] = HF_val
		# Sample Pearson correlation coefficient
		R_coeff = np.zeros(len(Hfactors))

		for HF in xrange(0,len(Hfactors)):
			top = HF_count[HF]*MF_HF[HF] - MF_tot[MF] * HF_tot[HF]
			bot = (HF_count[HF]*HF_tot_2[HF] - HF_tot[HF]**2)*(MF_count[MF]*MF_tot_2[MF] - MF_tot[MF]**2)
			R_coeff[HF] = top/(bot**(.5))
			print(R_coeff[HF])

		# Calculate averages
		if(MF_count[MF]>0):
			MF_AV[MF] = MF_tot[MF]/MF_count[MF]

		HF_AV = np.zeros(len(Hfactors))

		for HF in xrange(0, len(Hfactors)):
			if(HF_count[HF]>0):
				HF_AV[HF] = HF_tot[HF]/HF_count[HF]
		# Calculate sum of the difference squared
		HF_SDSQ = np.zeros(len(Hfactors))

		for pointID in xrange(0,numPts):
			MF_val = abs(model.GetPointData().GetArray('Threshold ' + Mfactors[MF] + ' (20)').GetValue(pointID))
			if(MF_val>0):
				MF_SDSQ[MF] = MF_SDSQ[MF] + (MF_val - MF_AV[MF]) ** 2
				for HF in xrange(0,len(Hfactors)):
					HF_vector = model.GetPointData().GetArray(Hfactors[HF]).GetTuple(pointID)
					HF_val = 0
					for i in xrange(0,len(HF_vector)):
						HF_val = HF_val + HF_vector[i]**2
					HF_val = HF_val**(.5)
					HF_SDSQ[HF] = HF_SDSQ[HF] + (HF_val - HF_AV[HF]) ** 2

		if(MF_count[MF]>0):
			MF_STDEV[MF] = (MF_SDSQ[MF]/MF_count[MF]) ** (.5)
		# Calculate standard deviation
		HF_STDEV = np.zeros(len(Hfactors))

		for HF in xrange(0, len(Hfactors)):
			if(HF_count[HF]>0):
				HF_STDEV[HF] = (HF_SDSQ[HF]/HF_count[HF]) ** (.5)
		# Set Threshold values
		MF_threshold = MF_AV
		HF_threshold = HF_AV
		MF_scale = MF_STDEV
		HF_scale = HF_STDEV
		MF_threshold[MF] = MF_THRESHOLD_FACTOR[TF] * MF_max[MF]
		MF_scale[MF] = MF_max[MF]
		#HF_threshold = .04 * HF_max
		#HF_scale[HF] = MF_max[MF]
		print(HF_AV/HF_max)
		print(MF_AV/MF_max)
		# Calculate correlation coeff at each point
		MF_HF_CORR = np.zeros((len(Hfactors), numPts))

		for pointID in xrange(0,numPts):
			MF_val = abs(model.GetPointData().GetArray('Threshold ' + Mfactors[MF] + ' (20)').GetValue(pointID))
			if(MF_val>0):
				for HF in xrange(0,len(Hfactors)):
					HF_vector = model.GetPointData().GetArray(Hfactors[HF]).GetTuple(pointID)
					HF_val = 0
					for i in xrange(0,len(HF_vector)):
						HF_val = HF_val + HF_vector[i]**2
					HF_val = HF_val**(.5)
					MF_HF_CORR[HF,pointID] = ((MF_val-MF_threshold[MF])*(HF_val-HF_threshold[HF])) \
					/(MF_scale[MF]*HF_scale[HF])
		# Add the correlation point array to the data
		for HF in xrange(0,len(Hfactors)):
			HfactorArray = vtk.vtkDoubleArray()
			for pointID in xrange(0,numPts):
				if(MF_HF_CORR[HF,pointID]>0 or MF_HF_CORR[HF,pointID]<0):
					HfactorArray.InsertNextValue(MF_HF_CORR[HF,pointID])
				else:
					HfactorArray.InsertNextValue(0)
			HfactorArray.SetName(Mfactors[MF] + ' - ' + Hfactors[HF] + ' - ' + str(MF_THRESHOLD_FACTOR[TF]*100)) 
			model.GetPointData().AddArray(HfactorArray)

w = vtk.vtkXMLPolyDataWriter()
w.SetInputData(model)
w.SetFileName('all_results_mapped_interpolated_threshold_correlated.vtp')
w.Write()
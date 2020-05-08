import sys
import vtk
import numpy as np

vtk_file_prefix = 'all_results_mapped_interpolated_threshold_gradients'
Mfactors = ['Dach1','Jag1']
Hfactors = ['vTAWSS (dynes/cm^2)','Gradients','pressure_avg (mmHg)']
# First, read in the .vtp file containing your quantities of interest
datareader=vtk.vtkXMLPolyDataReader()
datareader.SetFileName(vtk_file_prefix + '.vtp')
datareader.Update()
# Read your data into another polydata variable for manipulation
model=vtk.vtkPolyData()
model=datareader.GetOutput()

numPts=model.GetNumberOfPoints()

for MF in xrange(0,len(Mfactors)):
	output_filename = 'Point_correlation_' + Mfactors[MF] + '.dat'
	MF_EXP = np.arange(256)
	MF_HF = np.zeros((len(Hfactors),256))
	MF_HF_count = np.zeros((len(Hfactors),256))
	MF_HF_max = np.zeros((len(Hfactors),256))
	

	for pointID in xrange(0,numPts):
		MFValue = model.GetPointData().GetArray('Threshold ' + Mfactors[MF] + ' (20)').GetValue(pointID)
		Gradient = model.GetPointData().GetArray('Gradients').GetValue(pointID)
		if(MFValue>0 and Gradient < .2):
			for HF in xrange(0,len(Hfactors)):
				HF_vector = model.GetPointData().GetArray(Hfactors[HF]).GetTuple(pointID)
				HF_val = 0
				for i in xrange(0,len(HF_vector)):
					HF_val = HF_val + HF_vector[i]**2
				HF_val = HF_val**(.5)
				if(HF_val<200):
					MF_HF[HF, MFValue] = MF_HF[HF, MFValue] + HF_val
					MF_HF_count[HF, MFValue] = MF_HF_count[HF, MFValue] + 1
					if (HF_val>MF_HF_max[HF, MFValue]):
						MF_HF_max[HF, MFValue] = HF_val

	outfile = open(output_filename, 'w')

	
	MF_HF_AV = np.zeros((len(Hfactors), 256))

	for i in xrange(0, 255):
		for j in xrange(0, len(Hfactors)):
			if(MF_HF_count[j,i]>3):
				average = MF_HF[j,i]/MF_HF_count[j,i]
				MF_HF_AV[j,i] = average
# First print a header that tells what each integrated quantity of interest is
	header = ''
	subheader = ''
	for i in xrange(0, len(Hfactors)):
		header = header + Hfactors[i] + ', , '
		subheader = subheader + 'Average, Max, '
	header = header + '\n'
	subheader = subheader + '\n'
	outfile.write(header)
	outfile.write(subheader)

	
	for i in xrange(0, 255):
		out_string = str(MF_EXP[i])
		for j in xrange(0, len(Hfactors)):
			if(MF_HF_AV[j,i]>0):
				out_string = out_string + ', ' + str(MF_HF_AV[j,i]) 
				out_string = out_string + ', ' + str(MF_HF_max[j,i])
		if(out_string != str(MF_EXP[i])):
			out_string = out_string + '\n'
			outfile.write(out_string)

	outfile.close()


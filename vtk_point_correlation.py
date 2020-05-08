import sys
import vtk
import numpy as np

vtk_file_prefix = 'all_results_mapped_interpolated_threshold_gradients'
output_filename = 'out.dat'
output_collection = []
dach1_EXP = np.arange(256)
dach1_WSS = np.zeros(256)
dach1_WSSG = np.zeros(256)
jag1_EXP = np.arange(256)
jag1_WSS = np.zeros(256)
jag1_WSSG = np.zeros(256)
jag1_WSS_count = np.zeros(256)
dach1_WSS_count = np.zeroes(256)
jag1_WSSG_count = np.zeros(256)
dach1_WSSG_count = np.zeroes(256)
# First, read in the .vtp file containing your quantities of interest
datareader=vtk.vtkXMLPolyDataReader()
datareader.SetFileName(vtk_file_prefix + '.vtp')
datareader.Update()

# Read your data into another polydata variable for manipulation
model=vtk.vtkPolyData()
model=datareader.GetOutput()

numPts=model.GetNumberOfPoints()

for pointID in xrange(0,numPts):
	dach1Value = model.GetPointData().GetArray('Threshold Dach1 (10)').GetValue(pointID)
	if(dach1Value>0):
		dach1_WSS[dach1Value] = dach1_WSS[dach1Value] + (model.GetPointData().GetArray('vWSS').GetValue(pointID))
		dach1_WSS_count[dach1Value] = dach1_WSS_count[dach1Value] + 1
		if((model.GetPointData().GetArray('Gradients').GetValue(pointID))<200):
			dach1_WSSG[dach1Value] = dach1_WSSG[dach1Value] + (model.GetPointData().GetArray('Gradients').GetValue(pointID))
			dach1_WSSG_count[dach1Value] = dach1_WSSG_count[dach1Value] + 1

	jag1Value = model.GetPointData().GetArray('Threshold Jag1 (10)').GetValue(pointID)
	if(jag1Value>0):
		jag1_WSS[jag1Value] = jag1_WSS[jag1Value] + (model.GetPointData().GetArray('vWSS').GetValue(pointID))
		jag1_WSS_count[jag1Value] = jag1_WSS_count[jag1Value] + 1
		if((model.GetPointData().GetArray('Gradients').GetValue(pointID))<200):
			jag1_WSSG[jag1Value] = jag1_WSSG[jag1Value] + (model.GetPointData().GetArray('Gradients').GetValue(pointID))
			jag1_WSSG_count[jag1Value] = jag1_WSSG_count[jag1Value] + 1

outfile = open(output_filename, 'w')

# First print a header that tells what each integrated quantity of interest is
dach1_WSS_AV = []
dach1_WSSG_AV = []
jag1_WSS_AV = []
jag1_WSSG_AV = []

for i in xrange(0, len(dach1_WSS)):
	if(dach1_WSS_count[i]>0):
		average = dach1_WSS[i]/dach1_WSS_count[i]
		dach1_WSS_AV.append(average)

for i in xrange(0, len(dach1_WSSG)):
	if(dach1_WSSG_count[i]>0):
		average = dach1_WSSG[i]/dach1_WSSG_count[i]
		dach1_WSSG_AV.append(average)

for i in xrange(0, len(jag1_WSS)):
	if(jag1_WSS_count[i]>0):
		average = jag1_WSS[i]/jag1_WSS_count[i]
		jag1_WSS_AV.append(average)

for i in xrange(0, len(jag1_WSSG)):
	if(jag1_WSSG_count[i]>0):
		average = jag1_WSSG[i]/jag1_WSSG_count[i]
		jag1_WSSG_AV.append(average)

for i in xrange(0, 255):
	out_string = ''
	out_string = out_string + 'Dach1, '
	if(i<len(dach1_WSS_AV)):
		out_string = out_string + str(dach1_EXP[i]) + ', ' + str(dach1_WSS[i]) + ', ' + str(dach1_WSSG[i])
	out_string = out_string + 'Jag1, '
	out_string = out_string + str(jag1_EXP[i]) + ', ' + str(jag1_WSS[i]) + ', ' + str(jag1_WSSG[i])
	out_string = out_string + '\n'
	outfile.write(out_string)


for i in xrange(0, len(jag1_EXP)):
	out_string = 'Jag1, '
	out_string = out_string + str(jag1_EXP[i]) + ', ' + str(jag1_WSS[i]) + ', ' + str(jag1_WSSG[i])
	out_string = out_string + '\n'
	outfile.write(out_string)

outfile.close()
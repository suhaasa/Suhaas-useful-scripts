import sys
import os
import vtk
import numpy as np
import time
import glob
import math
import argparse
import csv
import tqdm
import matplotlib.pyplot as plt

def intializeVTU(filename):
	datareader=vtk.vtkXMLUnstructuredGridReader()
	datareader.SetFileName(filename)
	datareader.Update()

	mesh=datareader.GetOutput()
	print('Loaded .vtu file.')
	return mesh

def writeVTU(mesh,filename):
	print('Writing .vtu file...')
	w = vtk.vtkXMLUnstructuredGridWriter()
	w.SetInputData(mesh)
	w.SetFileName(filename)
	w.Write()
	print('done.')

def determinePerfusionVolumes(heart,caps):
	numPts = heart.GetNumberOfPoints()

	#generate separate data array for each perfusion volume
	#calculate the volume of each perfused area of each cap
	volumes = []
	cap_pt_list = vtk.vtkIdList()
	for i in range(0,len(caps)):
		for ip in range(0,numPts):
			if(heart.GetPointData().GetArray('PerfusionVolumes').GetValue(ip)==i):
				cap_pt_list.InsertNextId(ip)
		Mass = extractRegionVolume(heart,cap_pt_list)
		volumes.append(Mass.GetVolume())
		cap_pt_list.Reset()

	return volumes

def extractRegionVolume(mesh,selection_nodes):
	#Intialize variables
	ids = vtk.vtkIdTypeArray()
	cell_nodes = vtk.vtkIdList()
	cell_vtk_Id_list = vtk.vtkIdList()
	cellIds = vtk.vtkIdTypeArray()
	ids.SetNumberOfComponents(1)
	
	#Determines the cells enclosed by selection_nodes (which are points)
	vprint('Number of nodes in this volume: ', selection_nodes.GetNumberOfIds())
	for i in range(0,selection_nodes.GetNumberOfIds()):
		ids.InsertNextValue(selection_nodes.GetId(i))
		mesh.GetPointCells(selection_nodes.GetId(i), cell_nodes)
		for j in range(0,cell_nodes.GetNumberOfIds()):
			cell_vtk_Id_list.InsertUniqueId(cell_nodes.GetId(j))

	#Converts the vtkIdList into vtkIdTypeArray
	for i in range(0,cell_vtk_Id_list.GetNumberOfIds()):
		cellIds.InsertNextValue(cell_vtk_Id_list.GetId(i))
	vprint('Number of cells in this volume: ', cell_vtk_Id_list.GetNumberOfIds())
	#Creates the selection object to extract the subset of cells from the mesh
	region=vtk.vtkExtractSelection()
	region.SetInputData(0,mesh)
	tempCells = vtk.vtkSelectionNode()
	tempCells.SetFieldType(vtk.vtkSelectionNode.CELL)
	tempCells.SetContentType(vtk.vtkSelectionNode.INDICES)
	tempCells.SetSelectionList(cellIds)
	tempSelection = vtk.vtkSelection()
	tempSelection.AddNode(tempCells)
	region.SetInputData(1,tempSelection)
	region.Update()

	#Outputs the mesh as an Mass object
	output = vtk.vtkPolyData()
	output.ShallowCopy(region.GetOutput())
	vprint(region.GetOutput().GetNumberOfCells())
	dssf = vtk.vtkDataSetSurfaceFilter()
	dssf.SetInputConnection(region.GetOutputPort())
	dssf.Update()
	Mass = vtk.vtkMassProperties()
	Mass.SetInputData(dssf.GetOutput())
	Mass.Update()
	return Mass

def createParser():
	parser = argparse.ArgumentParser(description='Maps flow from each outlet to heart tissue.')
	parser.add_argument('heart', type=str, help='the input model heart')
	parser.add_argument('flow_data_loc', type=str, help='the flow_data filename location (include file ext)')
	parser.add_argument('out', type=str, help='the output filename (include file ext)')
	parser.add_argument('-v', '-verbose', type=int, nargs='?', const=1, default=0, help='turn on verbosity')
	return parser

def labelHeart(heart,flows,norm_flows,name_of_array):
	data_arrays = []
	numArrays = heart.GetPointData().GetNumberOfArrays()
	numPts = heart.GetNumberOfPoints()
	for i in range(0,numArrays):
		data_arrays.append(heart.GetPointData().GetArrayName(i))
	flow_ids = dict()
	for array in data_arrays:
		prefix = array.split('_')[0]
		if prefix.isdigit():
			#array_index.append(int(array.split('_')[0]))
			#vprint('a=',array[len(array)-len(caps[i]):],'c=',caps[i])
			flow_ids[int(array.split('_')[0])] = '_'.join(array.split('_')[1:])
			#vprint(flows)
	# get volumes of each region
	#vprint(len(array_flows))
	#vprint(array_flows)
	vtk_data_array = vtk.vtkDoubleArray()
	#volumes = determinePerfusionVolumes(heart,array_names)
	#tot_vol = 0
	#for i_vol in volumes:
#		tot_vol += i_vol
	#vprint(len(array_flows))
	#vprint(len(volumes))
	print(flows)
	flowsToIntegrate = {'LCA_11','LCA_11_1','LCA_11_2','LCA_11_3','LCA_11_4','LCA_11_5','LCA_11_6','LCA_11_7','LCA_12','LCA_13','LCA_14','LCA_14_1','LCA_14_2','LCA_14_3','LCA_14_4','LCA_14_5','LCA_14_6','LCA_14_7','LCA_14_8','LCA_14_9','LCA_14_10','LCA_14_11','LCA_14_12','LCA_14_1','LCA_14_14','LCA_14_15','LCA_14_16','LCA_14_17','LCA_15','LCA_16','LCA_16_1','LCA_16_2','LCA_16_3','LCA_16_4','LCA_16_5','LCA_16_6','LCA_16_7','LCA_16_8','LCA_16_9','LCA_16_10','LCA_16_11','LCA_16_12','LCA_17','LCA_18','LCA_19','LCA_20','LCA_21','LCA_22','LCA_23','LCA_24','LCA_24_1','LCA_24_2','LCA_24_3','LCA_25','LCA_26','LCA_27','LCA_28','LCA_29','LCA_main','LCA_main_75','LCA_main_80','LCA_main_85','LCA_main_90','LCA_main_95','LCA_main_99'}
	for ip in range(0,numPts):
		cap_index = heart.GetPointData().GetArray('caps_all').GetValue(ip)
		if(cap_index>=0 and flow_ids[int(cap_index)] in flows and flow_ids[int(cap_index)] in norm_flows):
			if norm_flows[flow_ids[int(cap_index)]]>0:
				vtk_data_array.InsertNextValue(flows[flow_ids[int(cap_index)]]/norm_flows[flow_ids[int(cap_index)]])
			else:
				vtk_data_array.InsertNextValue(-1)	
		else:
			vtk_data_array.InsertNextValue(-1)
	vtk_data_array.SetName(name_of_array+'_FlowPerfusion')
	heart.GetPointData().AddArray(vtk_data_array)
	writeToIschemicVolumeFile(heart,name_of_array,flow_ids,flows,norm_flows)

def writeToIschemicVolumeFile(heart,name_of_array,flow_ids,flows,norm_flows):
	thresholds = [.01,.02,.03,.04,.05,.06,.07,.08,.09,.1,.12,.14,.16,.18,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95,.96,.97,.98,.99,1]
	numPts = heart.GetNumberOfPoints()
	#get flow ids and flow percent recovery
	flow_recovery_fraction = dict()

	for i in norm_flows:
		if i in flows and norm_flows[i]>0:
			flow_recovery_fraction[i] = flows[i]/norm_flows[i]
	#if percent recovery lower than threshold, then assign it to ischemic volume and add it to growing sum
	ischemic_volume = dict()
	for it in thresholds:
		ischemic_nodes = vtk.vtkIdList()
		for ip in range(0,numPts):
			if heart.GetPointData().GetArray(name_of_array+'_FlowPerfusion').GetValue(ip) < it and heart.GetPointData().GetArray(name_of_array+'_FlowPerfusion').GetValue(ip) > 0:
				ischemic_nodes.InsertNextId(ip)
		region = vtk.vtkThreshold()
		region.SetInputData(heart)
		region.SetInputArrayToProcess(0,0,0,vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,name_of_array+'_FlowPerfusion')
		region.ThresholdBetween(0,it)
		region.Update()
		dssf = vtk.vtkDataSetSurfaceFilter()
		dssf.SetInputConnection(region.GetOutputPort())
		dssf.Update()
		Mass = vtk.vtkMassProperties()
		Mass.SetInputData(dssf.GetOutput())
		Mass.Update()
		ischemic_volume[it] = Mass.GetVolume()
		ischemic_nodes.Reset()

	file = open('Summary_ischemic_collateral.csv','a+')
	out_string = name_of_array
	for i in thresholds:
		out_string = out_string + ',' + str(ischemic_volume[i])
	file.write(out_string + '\n')
	file.close()


def writeToFlowFile(Sim_name,flows):
	with open('Summary_flows_collateral.csv') as f:
		first_line = f.readline() 
	first_line = first_line.split(',')[1:]
	out_list = []
	for i in first_line:
		if i in flows:
			out_list.append(str(flows[i]))
		else:
			out_list.append('N/A')
	file = open('Summary_flows_collateral.csv','a+')
	out_string = Sim_name
	for i in out_list:
		out_string = out_string + ',' + i
	file.write(out_string + '\n')
	file.close()

def main(args):
	FLOW_DIR = args.flow_data_loc
	heart = intializeVTU(args.heart)
	global vprint
	if args.v:
	    def vprint(*args):
	      # Print each argument separately so caller doesn't need to
	      # stuff everything to be printed into a single string
	      for arg in args:
	        print(arg),
	else:
		vprint = lambda *a: None

	# vtp files of the caps to calculate
	# cap_names = []
	# for file in os.listdir(args.caps):
	# 	if file.endswith('.vtp'):
	# 		cap_names.append(file[:len(file)-4])

	#vprint('Found ' + str(len(cap_names)) + ' caps.')
	# read csv file
	for flow_file in os.listdir(FLOW_DIR):
		if flow_file.split('_')[2]=='Sim00':
			#read flow file
			norm_flows = dict()
			with open(os.path.join(FLOW_DIR,flow_file)) as f:
				for line in f:
					if(len(line.split('\t')[0:2])==2):
						(key, val) = line.split('\t')[0:2]
						norm_flows[key] = abs(float(val))

	for flow_file in os.listdir(FLOW_DIR):
		#read flow file
		flows = dict()
		with open(os.path.join(FLOW_DIR,flow_file)) as f:
			for line in f:
				if(len(line.split('\t')[0:2])==2):	
					(key, val) = line.split('\t')[0:2]
					flows[key] = abs(float(val))
		# make list of cap names and flow values
		
		# print(data)
		# for i in range(0,len(data)):
		# 	if(len(data[i][0])>0):
		# 		flows[data[i][0]] = abs(float(data[i][1]))
		# vprint(flows)
		#labelHeart
		#write downstream outlet flows 
	



		# for i_row row in data[1:]:
		# 	vprint(row[0].split('\t')[1:])
		# 	flow_names = row[0].split('\t')[0]
		# 	for i in row[0].split('\t')[1:]
		# 		flows.append(float(i))
	
		#label heart
		labelHeart(heart,flows,norm_flows,flow_file.split('_')[2])
		writeToFlowFile(flow_file.split('_')[2],flows)
		#plt.hist(list(flows.values()),bins='auto')	
		#plt.savefig(flow_file.split('_')[2]+'.png')

	
	writeVTU(heart,args.out)
	# vprint(len(cap_names))
	# vprint(len(flows))
	flows_dict = {}
	# for i in range(0,len(cap_names)):
	# 	flows_dict[cap_names[i]] = flows[0][i]

	return flows_dict

if __name__ == '__main__':
	parser = createParser()
	args = parser.parse_args()
	main(args)

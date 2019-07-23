import sys
import os
import vtk
import numpy as np
import time
import glob
import math
import argparse
import csv

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
	parser.add_argument('flow_data', type=str, help='the flow_data filename (include file ext)')
	parser.add_argument('out', type=str, help='the output filename (include file ext)')
	parser.add_argument('-v', '-verbose', type=int, nargs='?', const=1, default=0, help='turn on verbosity')
	parser.add_argument('caps', type=str, help='the input model cap locations')
	return parser

def labelHeart(heart,flow_names,flows,caps):
	data_arrays = []
	numArrays = heart.GetPointData().GetNumberOfArrays()
	numPts = heart.GetNumberOfPoints()
	for i in range(0,numArrays):
		data_arrays.append(heart.GetPointData().GetArrayName(i))
	array_index = []
	array_names = []
	array_flows = []
	for array in data_arrays:
		try:
			array_index.append(int(array.split('_')[0]))
			for i in range(0,len(caps)):
				vprint('a=',array[len(array)-len(caps[i]):])
				if(array[len(array)-len(caps[i]):]==caps[i]):
					vprint('a=',array[len(array)-len(caps[i]):],'c=',caps[i])
					array_names.append(caps[i])
					vprint(flows)
					array_flows.append(flows[0][i])	
		except ValueError:
			vprint()
	# get volumes of each region
	vprint(len(array_flows))
	vprint(array_flows)
	vtk_data_array = vtk.vtkDoubleArray()
	volumes = determinePerfusionVolumes(heart,array_names)
	tot_vol = 0
	for i_vol in volumes:
		tot_vol += i_vol
	vprint(len(array_flows))
	vprint(len(volumes))
	for ip in range(0,numPts):
		cap_index = heart.GetPointData().GetArray('PerfusionVolumes').GetValue(ip)
		if(cap_index>=0):
			vtk_data_array.InsertNextValue(array_flows[int(cap_index)]/volumes[int(cap_index)])
		else:
			vtk_data_array.InsertNextValue(-1)
	vtk_data_array.SetName('FlowPerfusion')
	heart.GetPointData().AddArray(vtk_data_array)


def main(args):
	DATA_FILENAME = args.flow_data
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
	cap_names = []
	for file in os.listdir(args.caps):
		if file.endswith('.vtp'):
			cap_names.append(file[:len(file)-4])

	vprint('Found ' + str(len(cap_names)) + ' caps.')
	# read csv file
	
	with open(DATA_FILENAME) as csvfile:
		csvreader = csv.reader(csvfile, delimiter=' ')
		data = []
		for row in csvreader:
			data.append(row[0].split('\t'))
		# make list of cap names and flow values
		caps = []
		for cap in data[0]:
			if(len(cap)>0):
				caps.append(cap)
		vprint(caps)
		flows = np.zeros((len(data[1:]),len(data[0])))
		for i in range(1,len(data)):
			flow_names = data[i][0]
			for j in range(1,len(data[i])):
				if(len(data[i][j])>0):
					flows[i-1][j] = float(data[i][j])

		vprint(flows)



		# for i_row row in data[1:]:
		# 	vprint(row[0].split('\t')[1:])
		# 	flow_names = row[0].split('\t')[0]
		# 	for i in row[0].split('\t')[1:]
		# 		flows.append(float(i))
	
	#label heart
	labelHeart(heart,flow_names,flows,caps)
	writeVTU(heart,args.out)
	vprint(len(cap_names))
	vprint(len(flows))
	flows_dict = {}
	for i in range(0,len(cap_names)):
		flows_dict[cap_names[i]] = flows[0][i]
	return flows_dict

if __name__ == '__main__':
	parser = createParser()
	args = parser.parse_args()
	main(args)
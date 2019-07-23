import sys
import os
import vtk
import numpy as np
import time
import glob
import math
import argparse
import csv
import labelFlows as lf
import vtk_VolumePerfusion as vp


def intializeVTP(filename):
  datareader=vtk.vtkXMLPolyDataReader()
  datareader.SetFileName(filename)
  datareader.Update()

  mesh=vtk.vtkDataSetMapper()
  mesh=datareader.GetOutput()
  vprint('Loaded %s file.' % filename)
  return mesh

def intializeVTU(filename):
	datareader=vtk.vtkXMLUnstructuredGridReader()
	datareader.SetFileName(filename)
	datareader.Update()

	mesh=datareader.GetOutput()
	vprint('Loaded %s file.' % filename)
	return mesh

def writeVTU(mesh,filename):
	print('Writing .vtu file...')
	w = vtk.vtkXMLUnstructuredGridWriter()
	w.SetInputData(mesh)
	w.SetFileName(filename)
	w.Write()
	print('done.')

def writeVTP(mesh,filename):
  w = vtk.vtkXMLUnstructuredDataWriter()
  w.SetInputData(mesh)
  w.SetFileName(filename + '.vtp')
  w.Write()

def calculateFlow(slice):
	vn=[0.0,0.0,0.0]
	slice_area=0.0
	slice_flow=0.0
	slice_pressure=0.0
	cell = slice.GetCell(0)
	p0 = cell.GetPoints().GetPoint(0)
	p1 = cell.GetPoints().GetPoint(1)
	p2 = cell.GetPoints().GetPoint(2)
	vtk.vtkTriangle().ComputeNormal(p0,p1,p2,vn)
	pressure_array = slice.GetPointData().GetArray("pressure")
	velocity_array = slice.GetPointData().GetArray("velocity")
	for cellId in range(slice.GetNumberOfCells()):
		cell = slice.GetCell(cellId)
		p0 = cell.GetPoints().GetPoint(0)
		p1 = cell.GetPoints().GetPoint(1)
		p2 = cell.GetPoints().GetPoint(2) 
		cell_area=vtk.vtkTriangle().TriangleArea(p0, p1, p2)
		slice_area=slice_area+cell_area
		nodeids = cell.GetPointIds()
		v1 = np.array(velocity_array.GetTuple3(nodeids.GetId(0)))
		v2 = np.array(velocity_array.GetTuple3(nodeids.GetId(1)))
		v3 = np.array(velocity_array.GetTuple3(nodeids.GetId(2)))  
		slice_flow=slice_flow+sum((v1+v2+v3)/3.0*vn)*cell_area # 1 point Gauss quad rule
		p0 = np.array(pressure_array.GetTuple(nodeids.GetId(0)))
		p1 = np.array(pressure_array.GetTuple(nodeids.GetId(1)))
		p2 = np.array(pressure_array.GetTuple(nodeids.GetId(2)))
		slice_pressure=slice_pressure+np.mean([p0,p1,p2])*cell_area
	slice_pressure=slice_pressure/slice_area
	return slice_flow

def isolateCap(model,cap,caps_location):
	cap_vtp = intializeVTP(caps_location + '/' + cap + '.vtp')
	cellIds = vtk.vtkIdTypeArray()
	cap_cell_set = set()
	for i_cell in range(0,cap_vtp.GetNumberOfCells()):
		cap_cell_set.add(cap_vtp.GetCellData().GetArray('GlobalElementID').GetValue(i_cell))
	for i_cell in range(0,model.GetNumberOfCells()):
		if(model.GetCellData().GetArray('GlobalElementID').GetValue(i_cell) in cap_cell_set):
			cellIds.InsertNextValue(i_cell)

	#Creates the selection object to extract the subset of cells from the mesh
	region=vtk.vtkExtractSelection()
	region.SetInputData(0,model)
	tempCells = vtk.vtkSelectionNode()
	tempCells.SetFieldType(vtk.vtkSelectionNode.CELL)
	tempCells.SetContentType(vtk.vtkSelectionNode.INDICES)
	tempCells.SetSelectionList(cellIds)
	tempSelection = vtk.vtkSelection()
	tempSelection.AddNode(tempCells)
	region.SetInputData(1,tempSelection)
	region.Update()

	#Outputs the mesh as an polyData object
	dssf = vtk.vtkDataSetSurfaceFilter()
	dssf.SetInputConnection(region.GetOutputPort())
	dssf.Update()
	return dssf.GetOutput()

def writeFlowFile(model,caps,caps_location,flow_data):
	outfile = open(flow_data, 'w')
	flows = {}
	# loop over each cap and write the outflow to a text file
	for i_cap in range(0,len(caps)):
		#isolate vtp model of cap from model
		cap = isolateCap(model,caps[i_cap],caps_location)
		#calculate flow of cap
		flow = abs(calculateFlow(cap))
		#write to file
		outfile.write(caps[i_cap] + '\t' + str(flow) + '\n')
		flows[caps[i_cap]] = flow
	outfile.close()
	return flows

def readVolumes(volume_filename):
	volumes_dict = {}
	with open(volume_filename) as csvfile:
		csvreader = csv.reader(csvfile, delimiter=' ')
		data = []
		for row in csvreader:
			data.append(row[0].split(','))
		# make list of cap names and volumes values	
		for i in range(1,len(data)):
			volumes_dict[data[i][0]] = float(data[i][1])
	return volumes_dict

# reads the solver file to assign resistance values to a list of sorted caps
def readResistances(cap_names,solver_filename):
	KW = 'Resistance Values:'
	with open(solver_filename) as file:
		line = file.readline()
		while line:
			if(line[:len(KW)]==KW):
				data = line.split('\t')[1:]
				# make dict of resistance values
				resistances_dict = {}
				resistances_dict['aorta_2'] = data[len(data)-1]
				data = data[:len(data)-1]
				for i in range(0,len(data)):
					resistances_dict[cap_names[i]] = float(data[i])

			line = file.readline()
	return resistances_dict

def Tune(cap_names,flows,target_perfusion,target_flow,tuning_param,solver_filename,new_solver_filename,volume_filename):
	KW = 'Resistance Values:'
	#Read Resistances
	resistances = readResistances(cap_names,solver_filename)
	#Read Volumes
	volumes = readVolumes(volume_filename)
	#calculate perfusions
	perfusions = {}
	for i in flows:
		if i in flows and i in volumes:
			perfusions[i] = flows[i]/volumes[i]
		elif i not in flows:
			print('ERROR: %s not found in flows.' % i)
		elif i not in volumes:
			print('ERROR: %s not found in volumes.' % i)
	#Scale resistances by perfusion errors (from target perfusion)
	scaled_resistances = {}
	for i in perfusions:
		vprint(resistances[i])
		vprint(target_perfusion)
		vprint(perfusions[i])
		scaled_resistances[i] = 1/(1/float(resistances[i])*(float(target_perfusion)/float(perfusions[i]))**2)
	scaled_resistances['aorta_2'] = resistances['aorta_2']
	#Write new solver.inp with updated resistances
	old_solv = open(solver_filename, 'r')
	new_solv = open(new_solver_filename,'w')
	
	line = old_solv.readline()
	while line:
		if(line[:len(KW)]==KW):
			line = line.split('\t')[0]
			for i in cap_names:
				line += '\t' + str(scaled_resistances[i])
			line += '\t' + str(scaled_resistances['aorta_2'])	
			line += '\n'
		new_solv.write(line)
		line = old_solv.readline()
	old_solv.close()
	new_solv.close()

def createParser():
	parser = argparse.ArgumentParser(description='Postprocess data.')
	parser.add_argument('heart', type=str, help='the input model heart (vtu)')
	parser.add_argument('model', type=str, help='the input model vessels (vtp)')
	parser.add_argument('flow_data', type=str, help='the flow_data filename (include file ext)')
	parser.add_argument('out', type=str, help='the output filename (include file ext)')
	parser.add_argument('solver', type=str, help='the solver filename (include file ext)')
	parser.add_argument('new_solver', type=str, help='the new solver filename (include file ext)')
	parser.add_argument('volumes', type=str, help='the volumes filename (include file ext)')
	parser.add_argument('caps_loc', type=str, help='the input model cap locations')
	parser.add_argument('target_perf', type=str, help='the target perfusion')
	parser.add_argument('target_flow', type=str, help='the target flow')
	parser.add_argument('tune_param', type=str, default=0, help='tune flow (0) or perfusion (1)')
	parser.add_argument('-v', '-verbose', type=int, nargs='?', const=1, default=0, help='turn on verbosity')
	return parser

def main(args):
	global vprint
	if args.v:
		def vprint(*args):
			# Print each argument separately so caller doesn't need to
			# stuff everything to be printed into a single string
			for arg in args:
				print(arg),	
	else:
		vprint = lambda *a: None
	cap_names = []
	caps = []	
	for file in os.listdir(args.caps_loc):
		if file.endswith('.vtp'):
			cap_names.append(file[:len(file)-4])
			cap = intializeVTP(args.caps_loc + file)
			caps.append(cap)
	cap_names.sort()

	flows = writeFlowFile(intializeVTP(args.model),cap_names,args.caps_loc,args.flow_data)
	# lf_parser = lf.createParser()
	# lf_args = lf_parser.parse_args([args.heart,args.flow_data,args.out,args.caps_loc])
	# flows = lf.main(lf_args)

	Tune(cap_names, flows,args.target_perf,args.target_flow,args.tune_param,args.solver,args.new_solver,args.volumes)

if __name__ == '__main__':
	parser = createParser()
	args = parser.parse_args()
	main(args)
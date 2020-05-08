#Setup_Sims
import sys
import os
import vtk
import numpy as np
import time
import glob
import math
import argparse
import csv
from tqdm import tqdm

def intializeVTP(filename):
  datareader=vtk.vtkXMLPolyDataReader()
  datareader.SetFileName(filename)
  datareader.Update()

  mesh=vtk.vtkDataSetMapper()
  mesh=datareader.GetOutput()
  #vprint('Loaded %s file.' % filename)
  return mesh

def intializeVTU(filename):
	datareader=vtk.vtkXMLUnstructuredGridReader()
	datareader.SetFileName(filename)
	datareader.Update()

	mesh=datareader.GetOutput()
	#vprint('Loaded %s file.' % filename)
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

def calculateFlowAndPressure(slice):
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
	return slice_flow,slice_pressure

def isolateCap(model,cap,caps_location):
	cap_vtp = intializeVTP(caps_location + '/' + cap + '.vtp')
	cellIds = vtk.vtkIdTypeArray()
	cap_cell_set = set()
	for i_cell in range(0,cap_vtp.GetNumberOfCells()):
		cap_cell_set.add(cap_vtp.GetCellData().GetArray('GlobalElementID').GetValue(i_cell))
	for i_cell in range(0,model.GetNumberOfCells()):
		if(model.GetCellData().GetArray('GlobalElementID').GetValue(i_cell) in cap_cell_set):
			cellIds.InsertNextValue(i_cell)
	#print(cellIds.GetNumberOfValues())
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
	return dssf.GetOutput(),cellIds

def calculateArea(cap):
	masser = vtk.vtkMassProperties()
	masser.SetInputData(cap)
	masser.Update()
	return masser.GetSurfaceArea()

def writeFlowFile(model,caps,caps_location,flow_data):
	outfile = open(flow_data, 'w')
	flows = {}
	pressures = {}
	# loop over each cap and write the outflow to a text file
	vprint('Analizing caps...')
	all_cap_global_ids = dict()
	for i_cap in tqdm(range(0,len(caps))):
		#isolate vtp model of cap from model
		cap,cap_global_ids = isolateCap(model,caps[i_cap],caps_location)
		all_cap_global_ids[caps[i_cap]] = cap
		#calculate area of cap
		area = abs(calculateArea(cap))
		#calculate flow of cap
		flow,pressure = calculateFlowAndPressure(cap)
		#write to file
		outfile.write(caps[i_cap] + '\t' + str(flow) + '\t' + str(pressure) + '\t' + str(area) + '\n')
		flows[caps[i_cap]] = flow
		pressures[caps[i_cap]] = pressure
	outfile.close()
	return flows,pressures,all_cap_global_ids

def writeFlowFile_capKnown(model,cap_global_ids,flow_data):
	outfile = open(flow_data, 'w')
	flows = {}
	pressures = {}
	# loop over each cap and write the outflow to a text file
	vprint('Analizing caps...')
	for i_cap in tqdm(cap_global_ids):
		#calculate flow of cap
		flow,pressure = calculateFlowAndPressure(cap_global_ids[i_cap])
		#write to file
		outfile.write(i_cap + '\t' + str(flow) + '\t' + str(pressure) + '\n')
		flows[caps[i_cap]] = flow
		pressures[caps[i_cap]] = pressure
	outfile.close()
	return flows,pressures

def createParser():
	parser = argparse.ArgumentParser(description='Postprocess downstream flow.')
	parser.add_argument('S', type=str, default='', help='the simulation directory (from cwd ref)')
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
	master_file = 'Simulations_log.log'
	#Simulations = ['./Sim03_v10_lcac85/', './Sim05_v10_lcac95/',	'./Sim06_v10_lcac99/','./Sim20_v10_lcac_6col_20um/','./Sim23_v10_lcac85_6col_20um/','./Sim25_v10_lcac95_6col_20um/','./Sim26_v10_lcac99_6col_20um/','./Sim26_v10_lcac99_6col_20um/','./Sim30_v10_lcac_12col_20um/',	'./Sim30_v10_lcac_12col_20um/',	'./Sim33_v10_lcac85_12col_20um/',	'./Sim35_v10_lcac95_12col_20um/',	'./Sim36_v10_lcac99_12col_20um/',	'./Sim40_v10_lcac_6col_28um/',	'./Sim43_v10_lcac85_6col_28um/',	'./Sim45_v10_lcac95_6col_28um/',	'./Sim46_v10_lcac99_6col_28um/',	'./Sim50_v10_lcac_3col_40um/',	'./Sim50_v10_lcac_3col_40um/',	'./Sim50_v10_lcac_3col_40um/',	'./Sim53_v10_lcac85_3col_40um/',	'./Sim55_v10_lcac95_3col_40um/',	'./Sim56_v10_lcac99_3col_40um/',	'./Sim60_v10_lcac_9col_20um/',	'./Sim63_v10_lcac85_9col_20um/',	'./Sim65_v10_lcac95_9col_20um/',	'./Sim66_v10_lcac99_9col_20um/']
	#Simulations = [	'./Sim33_v10_lcac85_12col_20um/',	'./Sim35_v10_lcac95_12col_20um/',	'./Sim36_v10_lcac99_12col_20um/',	'./Sim40_v10_lcac_6col_28um/',	'./Sim43_v10_lcac85_6col_28um/',	'./Sim45_v10_lcac95_6col_28um/',	'./Sim46_v10_lcac99_6col_28um/',	'./Sim50_v10_lcac_3col_40um/',	'./Sim50_v10_lcac_3col_40um/',	'./Sim50_v10_lcac_3col_40um/',	'./Sim53_v10_lcac85_3col_40um/',	'./Sim55_v10_lcac95_3col_40um/',	'./Sim56_v10_lcac99_3col_40um/',	'./Sim60_v10_lcac_9col_20um/',	'./Sim63_v10_lcac85_9col_20um/',	'./Sim65_v10_lcac95_9col_20um/',	'./Sim66_v10_lcac99_9col_20um/']
	Simulations = os.path.join(os.getcwd(),args.S)
	Models = []
	for file in os.listdir(Simulations):
				if file.endswith('.vtp') and file.startswith('all_'):
					Models.append(os.path.join(Simulations,file))
	Models.sort()
	model_cap_gIds = dict()
	for model in Models:
		Sim_name = os.path.basename(model).split('_')
		Sim_name = Sim_name[2:len(Sim_name)-1]
		print(Sim_name)
		caps_loc = os.path.join(os.getcwd(),'Collateral_Sims','_'.join(Sim_name),'mesh-complete','mesh-surfaces')
		flowsToIntegrate = {'LCA_11','LCA_11_1','LCA_11_2','LCA_11_3','LCA_11_4','LCA_11_5','LCA_11_6','LCA_11_7','LCA_12','LCA_13','LCA_14','LCA_14_1','LCA_14_2','LCA_14_3','LCA_14_4','LCA_14_5','LCA_14_6','LCA_14_7','LCA_14_8','LCA_14_9','LCA_14_10','LCA_14_11','LCA_14_12','LCA_14_1','LCA_14_14','LCA_14_15','LCA_14_16','LCA_14_17','LCA_15','LCA_16','LCA_16_1','LCA_16_2','LCA_16_3','LCA_16_4','LCA_16_5','LCA_16_6','LCA_16_7','LCA_16_8','LCA_16_9','LCA_16_10','LCA_16_11','LCA_16_12','LCA_17','LCA_18','LCA_19','LCA_20','LCA_21','LCA_22','LCA_23','LCA_24','LCA_24_1','LCA_24_2','LCA_24_3','LCA_25','LCA_26','LCA_27','LCA_28','LCA_29','LCA_main','LCA_main_75','LCA_main_80','LCA_main_85','LCA_main_90','LCA_main_95','LCA_main_99'}
		aorta = {'aorta_v2','aorta_v2_2'}
		cap_names = []
		caps = []	
		for file in os.listdir(caps_loc):
			if file.endswith('.vtp') and not file.startswith('wall') and not file.startswith('C'):
				cap_names.append(file[:len(file)-4])
				cap = intializeVTP(os.path.join(caps_loc,file))
				caps.append(cap)
		cap_names.sort()
		#if Sim_name in model_cap_gIds.keys():
		#	flows,pressures = writeFlowFile_capKnown(intializeVTP(model),model_cap_gIds[Sim_name],os.path.basename(model).split('.')[0]+'_flow-pressure')
		#else:
		flows,pressures,all_cap_global_ids = writeFlowFile(intializeVTP(model),cap_names,caps_loc,os.path.normcase(os.path.join(os.getcwd(),'Collateral_Sims','_'.join(Sim_name),os.path.basename(model).split('.')[0]+'_flow-pressure')))
		#model_cap_gIds[Sim_name] = all_cap_global_ids
		downstream_flow = 0
		tot_flow = 0
		LCA_flow = 0
		RCA_flow = 0
		RSA_flow = 0
		for i in flows:
			if i in flowsToIntegrate:
				downstream_flow += flows[i]
			if i not in aorta:
				tot_flow += flows[i]
			if i.startswith('LCA'):
				LCA_flow += flows[i]
			if i.startswith('RCA'):
				RCA_flow += flows[i]
			if i.startswith('RSA'):
				RSA_flow += flows[i]
		file = open('Summary_flows.csv','a+')
		file.write(os.path.basename(model) + ','+ str(downstream_flow) + ','+ str(tot_flow) + ',' + str(LCA_flow) + ',' + str(RCA_flow) + ',' + str(RSA_flow) + ',' + str(RSA_flow+RCA_flow) + ',' + str(flows['aorta_v2']) + ',' + str(flows['aorta_v2_2']) + '\n')
		file.close()

if __name__ == '__main__':
	parser = createParser()
	args = parser.parse_args()
	main(args)

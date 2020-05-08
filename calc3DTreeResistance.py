from collections import defaultdict
import vtk
from tqdm import tqdm
import argparse
import os
from vtk_functions import read_geo, write_geo, calculator, cut_plane, connectivity, get_points_cells, clean, Integration
import numpy as np
import pickle

# Read a vtp file and return the polydata
def read_polydata(filename, datatype=None):
    """
    Load the given file, and return a vtkPolyData object for it.
    Args:
        filename (str): Path to input file.
        datatype (str): Additional parameter for vtkIdList objects.
    Returns:
        polyData (vtkSTL/vtkPolyData/vtkXMLStructured/
                    vtkXMLRectilinear/vtkXMLPolydata/vtkXMLUnstructured/
                    vtkXMLImage/Tecplot): Output data.
    """

    # Check if file exists
    if not os.path.exists(filename):
        raise RuntimeError("Could not find file: %s" % filename)

    # Check filename format
    fileType = filename.split(".")[-1]
    if fileType == '':
        raise RuntimeError('The file does not have an extension')

    # Get reader
    if fileType == 'stl':
        reader = vtk.vtkSTLReader()
        reader.MergingOn()
    elif fileType == 'vtk':
        reader = vtk.vtkPolyDataReader()
    elif fileType == 'vtp':
        reader = vtk.vtkXMLPolyDataReader()
    elif fileType == 'vts':
        reader = vtk.vtkXMinkorporereLStructuredGridReader()
    elif fileType == 'vtr':
        reader = vtk.vtkXMLRectilinearGridReader()
    elif fileType == 'vtu':
        reader = vtk.vtkXMLUnstructuredGridReader()
    elif fileType == "vti":
        reader = vtk.vtkXMLImageDataReader()
    elif fileType == "np" and datatype == "vtkIdList":
        result = np.load(filename).astype(np.int)
        id_list = vtk.vtkIdList()
        id_list.SetNumberOfIds(result.shape[0])
        for i in range(result.shape[0]):
            id_list.SetId(i, result[i])
        return id_list
    else:
        raise RuntimeError('Unknown file type %s' % fileType)

    # Read
    reader.SetFileName(filename)
    reader.Update()
    polydata = reader.GetOutput()

    return polydata

def write_polydata(input_data, filename, datatype=None):
    """
    Write the given input data based on the file name extension.
    Args:
        input_data (vtkSTL/vtkPolyData/vtkXMLStructured/
                    vtkXMLRectilinear/vtkXMLPolydata/vtkXMLUnstructured/
                    vtkXMLImage/Tecplot): Input data.
        filename (str): Save path location.
        datatype (str): Additional parameter for vtkIdList objects.
    """
    # Check filename format
    fileType = filename.split(".")[-1]
    if fileType == '':
        raise RuntimeError('The file does not have an extension')

    # Get writer
    if fileType == 'stl':
        writer = vtk.vtkSTLWriter()
    elif fileType == 'vtk':
        writer = vtk.vtkPolyDataWriter()
    elif fileType == 'vts':
        writer = vtk.vtkXMLStructuredGridWriter()
    elif fileType == 'vtr':
        writer = vtk.vtkXMLRectilinearGridWriter()
    elif fileType == 'vtp':
        writer = vtk.vtkXMLPolyDataWriter()
    elif fileType == 'vtu':
        writer = vtk.vtkXMLUnstructuredGridWriter()
    elif fileType == "vti":
        writer = vtk.vtkXMLImageDataWriter()
    elif fileType == "np" and datatype == "vtkIdList":
        output_data = np.zeros(input_data.GetNumberOfIds())
        for i in range(input_data.GetNumberOfIds()):
            output_data[i] = input_data.GetId(i)
        output_data.dump(filename)
        return
    else:
        raise RuntimeError('Unknown file type %s' % fileType)

    # Set filename and input
    writer.SetFileName(filename)
    writer.SetInputData(input_data)
    writer.Update()

    # Write
    writer.Write()

def calculateFlowPressureAverage(slice,start,stop,step):
	vn=[0.0,0.0,0.0]
	slice_pressure = dict()
	slice_area = dict()
	slice_flow = dict()
	cell = slice.GetCell(0)
	p0 = cell.GetPoints().GetPoint(0)
	p1 = cell.GetPoints().GetPoint(1)
	p2 = cell.GetPoints().GetPoint(2)
	vtk.vtkTriangle().ComputeNormal(p0,p1,p2,vn)
	num_timesteps = (stop-start+step)/step
	timesteps = np.linspace(start,stop,num_timesteps)
	for t in timesteps:
		slice_area[t]=0.0
		slice_flow[t]=0.0
		#print(str(int(t)))
		slice_pressure[t]=0.0
		pressure_array = slice.GetPointData().GetArray("pressure_"+str(int(t)))
		velocity_array = slice.GetPointData().GetArray("velocity_"+str(int(t)))
		for cellId in range(slice.GetNumberOfCells()):
			cell = slice.GetCell(cellId)
			p0 = cell.GetPoints().GetPoint(0)
			p1 = cell.GetPoints().GetPoint(1)
			p2 = cell.GetPoints().GetPoint(2) 
			cell_area=vtk.vtkTriangle().TriangleArea(p0, p1, p2)
			slice_area[t]=slice_area[t]+cell_area
			nodeids = cell.GetPointIds()
			v1 = np.array(velocity_array.GetTuple3(nodeids.GetId(0)))
			v2 = np.array(velocity_array.GetTuple3(nodeids.GetId(1)))
			v3 = np.array(velocity_array.GetTuple3(nodeids.GetId(2)))  
			slice_flow[t]=slice_flow[t]+sum((v1+v2+v3)/3.0*vn)*cell_area # 1 point Gauss quad rule
			p0 = np.array(pressure_array.GetTuple(nodeids.GetId(0)))
			p1 = np.array(pressure_array.GetTuple(nodeids.GetId(1)))
			p2 = np.array(pressure_array.GetTuple(nodeids.GetId(2)))
			slice_pressure[t]=slice_pressure[t]+np.mean([p0,p1,p2])*cell_area
		slice_pressure[t]=slice_pressure[t]/slice_area[t]
		# print(slice_pressure[t])
		# print(slice_flow[t])
	return slice_flow,slice_pressure,slice_area

def calculateFlowPressure(slice):
	vn=[0.0,0.0,0.0]
	slice_pressure = dict()
	slice_area = dict()
	slice_flow = dict()
	cell = slice.GetCell(0)
	p0 = cell.GetPoints().GetPoint(0)
	p1 = cell.GetPoints().GetPoint(1)
	p2 = cell.GetPoints().GetPoint(2)
	vtk.vtkTriangle().ComputeNormal(p0,p1,p2,vn)
	num_timesteps = 1
	timesteps = [0]
	for t in timesteps:
		slice_area[t]=0.0
		slice_flow[t]=0.0
		#print(str(int(t)))
		slice_pressure[t]=0.0
		pressure_array = slice.GetPointData().GetArray("pressure")
		velocity_array = slice.GetPointData().GetArray("velocity")
		for cellId in range(slice.GetNumberOfCells()):
			cell = slice.GetCell(cellId)
			p0 = cell.GetPoints().GetPoint(0)
			p1 = cell.GetPoints().GetPoint(1)
			p2 = cell.GetPoints().GetPoint(2) 
			cell_area=vtk.vtkTriangle().TriangleArea(p0, p1, p2)
			slice_area[t]=slice_area[t]+cell_area
			nodeids = cell.GetPointIds()
			v1 = np.array(velocity_array.GetTuple3(nodeids.GetId(0)))
			v2 = np.array(velocity_array.GetTuple3(nodeids.GetId(1)))
			v3 = np.array(velocity_array.GetTuple3(nodeids.GetId(2)))  
			slice_flow[t]=slice_flow[t]+sum((v1+v2+v3)/3.0*vn)*cell_area # 1 point Gauss quad rule
			p0 = np.array(pressure_array.GetTuple(nodeids.GetId(0)))
			p1 = np.array(pressure_array.GetTuple(nodeids.GetId(1)))
			p2 = np.array(pressure_array.GetTuple(nodeids.GetId(2)))
			slice_pressure[t]=slice_pressure[t]+np.mean([p0,p1,p2])*cell_area
		slice_pressure[t]=slice_pressure[t]/slice_area[t]
		# print(slice_pressure[t])
		# print(slice_flow[t])
	return slice_flow,slice_pressure,slice_area

def getConnectedVerticesNotIncludingSeed(model, seedPt):
    cell_list = vtk.vtkIdList()
    connectedPts_list = vtk.vtkIdList()
    model.GetPointCells(seedPt,cell_list)
    for j in range(0,cell_list.GetNumberOfIds()):
        pt_list = vtk.vtkIdList()
        pt_list = model.GetCell(cell_list.GetId(j)).GetPointIds()
        for k in range(0,pt_list.GetNumberOfIds()):
            if (pt_list.GetId(k) != seedPt):
                connectedPts_list.InsertUniqueId(pt_list.GetId(k))
    return connectedPts_list

def extractRegion(model,value,array_name):
	regionIds = vtk.vtkIdList()
	array = model.GetPointData().GetArray(array_name)
	numPts = model.GetNumberOfPoints()
	for i in range(0,numPts):
		if value == array.GetValue(i):
			regionIds.InsertUniqueId(i)
	return regionIds

def extractCellRegion(model,value,array_name):
	regionIds = vtk.vtkIdList()
	array = model.GetCellData().GetArray(array_name)
	numPts = model.GetNumberOfCells()
	for i in range(0,numPts):
		if value == array.GetValue(i):
			regionIds.InsertUniqueId(i)
	return regionIds

def getBranchConnectivity(model):
	connectivity = defaultdict(list)
	[BranchId_min, BranchId_max] = model.GetPointData().GetArray('BifurcationId').GetRange()
	Bifurcation_PtIds = vtk.vtkIdList()
	#loop through each bifurcationId
	for b in tqdm(range(int(BranchId_min),int(BranchId_max+1))):
		Bifurcation_PtIds = extractRegion(model,b,'BifurcationId')
		num_different_branchIds = set()
		#loop through each point in bifurcation region
		for pt in tqdm(range(0,Bifurcation_PtIds.GetNumberOfIds())):
			#find all neighbor points for each bifurcation point and see if it's connected to a different branchId
			neighbor_pts = getConnectedVerticesNotIncludingSeed(model, Bifurcation_PtIds.GetId(pt))
			for neighbor in range(0,neighbor_pts.GetNumberOfIds()):
				neighbor_branchId = model.GetPointData().GetArray('BranchId').GetValue(neighbor_pts.GetId(neighbor))
				if neighbor_branchId >= 0 and neighbor_branchId not in num_different_branchIds:
					num_different_branchIds.add(neighbor_branchId)
		#Assumes parent branch has the lowest BranchId
		if len(list(num_different_branchIds))>0:
			parent_branch = list(num_different_branchIds)[0]
			child_branches = []
			for branchId in num_different_branchIds:
				if branchId <= parent_branch:
					parent_branch = branchId
			num_different_branchIds.remove(parent_branch)
			child_branches = list(num_different_branchIds)
			if len(child_branches)>0:	
				connectivity[parent_branch] = child_branches
			else:
				print('WARNING: Parent branch '+str(int(parent_branch))+' has no child branches.')
		else:
			print('WARNING:Ignoring bifurcation '+str(int(b)))
	return connectivity

def getResistance(branch,centerline,connectivity,branch_resistances):
	global resistance
	global parallel_res
	#checks to see if parent branch is in connectivity dictionary
	if branch in connectivity.keys():
		#If yes, then loop through each child branch and recursively call the method for downstream branches summed in parallel
		parallel_res = 0
		for child_branch in connectivity[branch]:
			parallel_res += 1.0/getResistance(child_branch,centerline,connectivity,branch_resistances)
		resistance = 1.0/parallel_res + branch_resistances[branch]
	else:
		#Otherwise, return the resistance of the branch itself
		mapBranch_resistances({branch:branch_resistances[branch]},centerline,'Downstream_Resistances')
		return branch_resistances[branch]
	print({branch:resistance})
	mapBranch_resistances({branch:resistance},centerline,'Downstream_Resistances')
	return resistance

def getBranchResistances(model):
	branch_resistances = dict()
	#loop through each branchId
	[BranchId_min, BranchId_max] = model.GetPointData().GetArray('BranchId').GetRange()
	for branchId in range(int(BranchId_min),int(BranchId_max+1)):
		Branch_PtIds = extractRegion(model,branchId,'BranchId')
		#loop through each branch point and store pressure at each node
		pressures = dict()
		for pt in range(0,Branch_PtIds.GetNumberOfIds()):
			pressures[Branch_PtIds.GetId(pt)] = model.GetPointData().GetArray('pressure').GetValue(Branch_PtIds.GetId(pt))
		#loop through pressures to find min and max pressure
		pmax_node = list(pressures.keys())[0]
		pmin_node = list(pressures.keys())[0]
		pmax = pressures[list(pressures.keys())[0]]
		pmin = pressures[list(pressures.keys())[0]]
		for p in pressures:
			if p > pmax_node:
				pmax = pressures[p]
				pmax_node = p
			if p < pmin_node:
				pmin = pressures[p]
				pmin_node = p
		#Calculate pressure difference and flow at the minimum pressure node to find resistance of branch
		branch_pdiff = pmax - pmin
		print(pmin_node)
		branch_flow = model.GetPointData().GetArray('flow').GetValue(pmin_node)
		print('pdiff',branch_pdiff)
		print('p1',pmax,'p2',pmin)
		print('flow',branch_flow)
		branch_resistances[branchId] = branch_pdiff/branch_flow
	return branch_resistances

def getmin_node(vtkIdList):
	min_node = vtkIdList.GetId(0)
	for i in range(0,vtkIdList.GetNumberOfIds()):
		if min_node > vtkIdList.GetId(i):
			min_node = vtkIdList.GetId(i)
	return min_node

def getmax_node(vtkIdList):
	max_node = vtkIdList.GetId(0)
	for i in range(0,vtkIdList.GetNumberOfIds()):
		if max_node < vtkIdList.GetId(i):
			max_node = vtkIdList.GetId(i)
	return max_node


def getBranchResistancesAverage(centerline,volume,connectivity,start,stop,step,SEG_COVERAGE):
	branch_resistances = dict()
	#loop through each branchId
	[BranchId_min, BranchId_max] = centerline.GetPointData().GetArray('BranchId').GetRange()
	for branchId in tqdm(connectivity):
		Branch_PtIds = extractRegion(centerline,branchId,'BranchId')
		parent_p1_node = getmin_node(Branch_PtIds)
		parent_p2_node = getmax_node(Branch_PtIds)
		slice_1 = slice_vessel(volume, centerline.GetPoint(parent_p1_node), centerline.GetPointData().GetArray('CenterlineSectionNormal').GetTuple(parent_p1_node))
		slice1_flow,slice1_pressure,slice1_area = calculateFlowPressureAverage(slice_1,start,stop,step)
		slice_2 = slice_vessel(volume, centerline.GetPoint(parent_p2_node), centerline.GetPointData().GetArray('CenterlineSectionNormal').GetTuple(parent_p2_node))
		slice2_flow,slice2_pressure,slice2_area = calculateFlowPressureAverage(slice_2,start,stop,step)
		print(connectivity[branchId])
		for childId in connectivity[branchId]:
			Branch_PtIds = extractRegion(centerline,childId,'BranchId')
			print(childId)
			child_p1_node = getmin_node(Branch_PtIds)
			child_p2_node = int((getmax_node(Branch_PtIds)-getmin_node(Branch_PtIds))*SEG_COVERAGE)+getmin_node(Branch_PtIds)
			child_slice_1 = slice_vessel(volume, centerline.GetPoint(child_p1_node), centerline.GetPointData().GetArray('CenterlineSectionNormal').GetTuple(child_p1_node))
			child_slice1_flow,child_slice1_pressure,child_slice1_area = calculateFlowPressureAverage(child_slice_1,start,stop,step)
			child_slice_2 = slice_vessel(volume, centerline.GetPoint(child_p2_node), centerline.GetPointData().GetArray('CenterlineSectionNormal').GetTuple(child_p2_node))
			child_slice2_flow,child_slice2_pressure,child_slice2_area = calculateFlowPressureAverage(child_slice_2,start,stop,step)
			write_polydata(child_slice_1,'./Slices/slice_'+str(int(childId))+'_1.vtp')
			write_polydata(child_slice_2,'./Slices/slice_'+str(int(childId))+'_2.vtp')
			#Calculate pressure difference and flow at the minimum pressure node to find resistance of branch
			branch_pdiff = time_average(slice2_pressure) - time_average(child_slice2_pressure)
			segment_pdiff = time_average(child_slice1_pressure) - time_average(child_slice2_pressure)
			if branch_pdiff < 0:
				print(parent_p2_node,child_p2_node,branch_pdiff,'less than 0')
			print(parent_p2_node,child_p2_node)
			branch_flow = time_average(child_slice2_flow)
			print('pdiff',branch_pdiff)
			print('p1',time_average(slice2_pressure),'p2',time_average(child_slice2_pressure))
			print('flow',branch_flow)
			mapBranch_resistances({childId:segment_pdiff/branch_flow},centerline,'Segment_Resistances')
			branch_resistances[childId] = branch_pdiff/branch_flow
	Branch_PtIds = extractRegion(centerline,0,'BranchId')
	parent_p1_node = getmin_node(Branch_PtIds)+1
	parent_p2_node = int((getmax_node(Branch_PtIds)-getmin_node(Branch_PtIds))*SEG_COVERAGE)+getmin_node(Branch_PtIds)
	slice_1 = slice_vessel(volume, centerline.GetPoint(parent_p2_node), centerline.GetPointData().GetArray('CenterlineSectionNormal').GetTuple(parent_p1_node))
	slice1_flow,slice1_pressure,slice1_area = calculateFlowPressureAverage(slice_1,start,stop,step)
	slice_2 = slice_vessel(volume, centerline.GetPoint(parent_p2_node), centerline.GetPointData().GetArray('CenterlineSectionNormal').GetTuple(parent_p2_node))
	slice2_flow,slice2_pressure,slice2_area = calculateFlowPressureAverage(slice_2,start,stop,step)
	branch_pdiff = time_average(slice1_pressure) - time_average(slice2_pressure)
	branch_flow = time_average(slice2_flow)
	branch_resistances[0] = branch_pdiff/branch_flow
	return branch_resistances

def getBranchResistances2Average(centerline,volume,connectivity,start,stop,step,SEG_COVERAGE):
	branch_resistances = dict()
	#loop through each branchId
	[BranchId_min, BranchId_max] = centerline.GetPointData().GetArray('BranchId').GetRange()
	for branchId in tqdm(connectivity):
		Branch_PtIds = extractRegion(centerline,branchId,'BranchId')
		parent_p1_node = getmin_node(Branch_PtIds)
		parent_p2_node = int(Branch_PtIds.GetNumberOfIds()*SEG_COVERAGE)+getmin_node(Branch_PtIds)
		slice_1 = slice_vessel_sphere(volume, centerline.GetPoint(parent_p1_node), centerline.GetPointData().GetArray('CenterlineSectionNormal').GetTuple(parent_p1_node),centerline.GetPointData().GetArray('MaximumInscribedSphereRadius').GetValue(parent_p1_node))
		slice1_flow,slice1_pressure,slice1_area = calculateFlowPressureAverage(slice_1,start,stop,step)
		slice_2 = slice_vessel_sphere(volume, centerline.GetPoint(parent_p2_node), centerline.GetPointData().GetArray('CenterlineSectionNormal').GetTuple(parent_p2_node),centerline.GetPointData().GetArray('MaximumInscribedSphereRadius').GetValue(parent_p2_node))
		slice2_flow,slice2_pressure,slice2_area = calculateFlowPressureAverage(slice_2,start,stop,step)
		#loop through each child of the parent branch
		for childId in connectivity[branchId]:
			Branch_PtIds = extractRegion(centerline,childId,'BranchId')
			print(childId)
			child_p1_node = getmin_node(Branch_PtIds)
			child_p2_node = int(Branch_PtIds.GetNumberOfIds()*SEG_COVERAGE)+getmin_node(Branch_PtIds)
			print(child_p1_node,child_p2_node,Branch_PtIds.GetNumberOfIds())
			child_slice_1 = slice_vessel_sphere(volume, centerline.GetPoint(child_p1_node), centerline.GetPointData().GetArray('CenterlineSectionNormal').GetTuple(child_p1_node),centerline.GetPointData().GetArray('MaximumInscribedSphereRadius').GetValue(child_p1_node))
			child_slice1_flow,child_slice1_pressure,child_slice1_area = calculateFlowPressureAverage(child_slice_1,start,stop,step)
			child_slice_2 = slice_vessel_sphere(volume, centerline.GetPoint(child_p2_node), centerline.GetPointData().GetArray('CenterlineSectionNormal').GetTuple(child_p2_node),centerline.GetPointData().GetArray('MaximumInscribedSphereRadius').GetValue(child_p2_node))
			child_slice2_flow,child_slice2_pressure,child_slice2_area = calculateFlowPressureAverage(child_slice_2,start,stop,step)
			write_polydata(child_slice_1,'./Slices/slice_'+str(int(childId))+'_1.vtp')
			write_polydata(child_slice_2,'./Slices/slice_'+str(int(childId))+'_2.vtp')
			#Calculate pressure difference and flow between p2 of parent and p2 of the child to find resistance of child branch
			branch_pdiff = time_average(slice2_pressure) - time_average(child_slice2_pressure)
			segment_pdiff = time_average(child_slice1_pressure) - time_average(child_slice2_pressure)
			branch_flow = time_average(child_slice2_flow)
			#print info
			if branch_pdiff < 0:
				print(parent_p2_node,child_p2_node,branch_pdiff,'less than 0')
			print(parent_p2_node,child_p2_node)
			print('pdiff',branch_pdiff)
			print('p1',time_average(slice2_pressure),'p2',time_average(child_slice2_pressure))
			print('flow',branch_flow)
			#Add segment resistance to the model child segment
			mapBranch_resistances({childId:segment_pdiff/branch_flow},centerline,'Segment_Resistances')
			#Store resistance
			branch_resistances[childId] = branch_pdiff/branch_flow
	#For inlet calculate the resistance of the initial segment by taking p1-p2/Qout
	Branch_PtIds = extractRegion(centerline,0,'BranchId')
	parent_p1_node = getmin_node(Branch_PtIds)+1
	parent_p2_node = int(Branch_PtIds.GetNumberOfIds()*SEG_COVERAGE)+getmin_node(Branch_PtIds)
	slice_1 = slice_vessel_sphere(volume, centerline.GetPoint(parent_p1_node), centerline.GetPointData().GetArray('CenterlineSectionNormal').GetTuple(parent_p1_node),centerline.GetPointData().GetArray('MaximumInscribedSphereRadius').GetValue(parent_p1_node))
	slice1_flow,slice1_pressure,slice1_area = calculateFlowPressureAverage(slice_1,start,stop,step)
	slice_2 = slice_vessel_sphere(volume, centerline.GetPoint(parent_p2_node), centerline.GetPointData().GetArray('CenterlineSectionNormal').GetTuple(parent_p2_node),centerline.GetPointData().GetArray('MaximumInscribedSphereRadius').GetValue(parent_p2_node))
	slice2_flow,slice2_pressure,slice2_area = calculateFlowPressureAverage(slice_2,start,stop,step)
	branch_pdiff = time_average(slice1_pressure) - time_average(slice2_pressure)
	branch_flow = time_average(slice2_flow)
	branch_resistances[0] = branch_pdiff/branch_flow
	#Return stored values
	return branch_resistances

def getBranchResistances2(centerline,volume,connectivity,SEG_COVERAGE):
	branch_resistances = dict()
	#loop through each branchId
	[BranchId_min, BranchId_max] = centerline.GetPointData().GetArray('BranchId').GetRange()
	for branchId in tqdm(connectivity):
		Branch_PtIds = extractRegion(centerline,branchId,'BranchId')
		if(branchId!=0):
			parent_p1_node = getmin_node(Branch_PtIds)
		else:
			parent_p1_node = getmin_node(Branch_PtIds)+30
		parent_p2_node = int(Branch_PtIds.GetNumberOfIds()*SEG_COVERAGE)+getmin_node(Branch_PtIds)
		slice_1 = slice_vessel_sphere(volume, centerline.GetPoint(parent_p1_node), centerline.GetPointData().GetArray('CenterlineSectionNormal').GetTuple(parent_p1_node),centerline.GetPointData().GetArray('MaximumInscribedSphereRadius').GetValue(parent_p1_node))
		slice1_flow,slice1_pressure,slice1_area = calculateFlowPressure(slice_1)
		slice_2 = slice_vessel_sphere(volume, centerline.GetPoint(parent_p2_node), centerline.GetPointData().GetArray('CenterlineSectionNormal').GetTuple(parent_p2_node),centerline.GetPointData().GetArray('MaximumInscribedSphereRadius').GetValue(parent_p2_node))
		slice2_flow,slice2_pressure,slice2_area = calculateFlowPressure(slice_2)
		for childId in connectivity[branchId]:
			Branch_PtIds = extractRegion(centerline,childId,'BranchId')
			print(childId)
			child_p1_node = getmin_node(Branch_PtIds)
			child_p2_node = int(Branch_PtIds.GetNumberOfIds()*SEG_COVERAGE)+getmin_node(Branch_PtIds)
			print(child_p1_node,child_p2_node,Branch_PtIds.GetNumberOfIds())
			child_slice_1 = slice_vessel_sphere(volume, centerline.GetPoint(child_p1_node), centerline.GetPointData().GetArray('CenterlineSectionNormal').GetTuple(child_p1_node),centerline.GetPointData().GetArray('MaximumInscribedSphereRadius').GetValue(child_p1_node))
			child_slice1_flow,child_slice1_pressure,child_slice1_area = calculateFlowPressure(child_slice_1)
			child_slice_2 = slice_vessel_sphere(volume, centerline.GetPoint(child_p2_node), centerline.GetPointData().GetArray('CenterlineSectionNormal').GetTuple(child_p2_node),centerline.GetPointData().GetArray('MaximumInscribedSphereRadius').GetValue(child_p2_node))
			child_slice2_flow,child_slice2_pressure,child_slice2_area = calculateFlowPressure(child_slice_2)
			write_polydata(child_slice_1,'./Slices/slice_1_'+str(int(childId))+'.vtp')
			write_polydata(child_slice_2,'./Slices/slice_2_'+str(int(childId))+'.vtp')
			#Calculate pressure difference and flow at the minimum pressure node to find resistance of branch
			branch_pdiff = time_average(slice2_pressure) - time_average(child_slice2_pressure)
			segment_pdiff = time_average(child_slice1_pressure) - time_average(child_slice2_pressure)
			if branch_pdiff < 0:
				print(parent_p2_node,child_p2_node,branch_pdiff,'less than 0')
			print(parent_p2_node,child_p2_node)
			branch_flow = time_average(child_slice2_flow)
			print('pdiff',branch_pdiff)
			print('p1',time_average(slice2_pressure),'p2',time_average(child_slice2_pressure))
			print('flow',branch_flow)
			mapBranch_resistances({childId:segment_pdiff/branch_flow},centerline,'Segment_Resistances')
			branch_resistances[childId] = branch_pdiff/branch_flow
	Branch_PtIds = extractRegion(centerline,0,'BranchId')
	parent_p1_node = getmin_node(Branch_PtIds)+30
	parent_p2_node = int(Branch_PtIds.GetNumberOfIds()*SEG_COVERAGE)+getmin_node(Branch_PtIds)
	slice_1 = slice_vessel_sphere(volume, centerline.GetPoint(parent_p1_node), centerline.GetPointData().GetArray('CenterlineSectionNormal').GetTuple(parent_p1_node),centerline.GetPointData().GetArray('MaximumInscribedSphereRadius').GetValue(parent_p1_node))
	slice1_flow,slice1_pressure,slice1_area = calculateFlowPressure(slice_1)
	slice_2 = slice_vessel_sphere(volume, centerline.GetPoint(parent_p2_node), centerline.GetPointData().GetArray('CenterlineSectionNormal').GetTuple(parent_p2_node),centerline.GetPointData().GetArray('MaximumInscribedSphereRadius').GetValue(parent_p2_node))
	slice2_flow,slice2_pressure,slice2_area = calculateFlowPressure(slice_2)
	branch_pdiff = time_average(slice1_pressure) - time_average(slice2_pressure)
	branch_flow = time_average(slice2_flow)
	branch_resistances[0] = branch_pdiff/branch_flow
	return branch_resistances

def time_average(value):
	avg = 0
	for t in value:
		avg += value[t]
	avg = avg/len(list(value.keys()))
	return avg

def slice_vessel(inp_3d, origin, normal):
	"""
	Slice 3d geometry at certain plane
	Args:
	    inp_1d: vtk InputConnection for 1d centerline
	    inp_3d: vtk InputConnection for 3d volume model
	    origin: plane origin
	    normal: plane normal

	Returns:
	    Integration object
	"""
	# cut 3d geometry
	cut_3d = cut_plane(inp_3d, origin, normal)
	# extract region closest to centerline
	con = connectivity(cut_3d, origin)
	dssf = vtk.vtkDataSetSurfaceFilter()
	dssf.SetInputConnection(con.GetOutputPort())
	dssf.Update()
	return dssf.GetOutput()

def clip_sphere(inp_3d,origin,radius):
	"""
	Clips geometry in sphere
	Args:
	inp: InputConnection
	origin: clipping sphere origin
	normal: clipping sphere radius
	Returns:
	cut: clipping object
	"""
	# define clipping sphere
	sphere = vtk.vtkSphere()
	sphere.SetCenter(origin)
	sphere.SetRadius(radius*1.1)

	# define cutter
	clip = vtk.vtkClipPolyData()
	clip.SetInputData(inp_3d.GetOutput())
	clip.SetClipFunction(sphere)
	clip.InsideOutOn()
	clip.Update()
	return clip

def slice_vessel_sphere(inp_3d, origin, normal, radius):
	"""
	Slice 3d geometry at certain plane
	Args:
	    inp_1d: vtk InputConnection for 1d centerline
	    inp_3d: vtk InputConnection for 3d volume model
	    origin: plane origin
	    normal: plane normal

	Returns:
	    Integration object
	"""
	# cut 3d geometry
	cut_3d = cut_plane(inp_3d, origin, normal)
	clip_3d = clip_sphere(cut_3d,origin, radius)
	write_polydata(cut_3d.GetOutput(),'plane_cut.vtp')
	write_polydata(clip_3d.GetOutput(),'sphere_cut.vtp')
	# extract region closest to centerline
	con = connectivity(clip_3d, origin)
	dssf = vtk.vtkDataSetSurfaceFilter()
	dssf.SetInputConnection(con.GetOutputPort())
	dssf.Update()
	return dssf.GetOutput()

def mapVolumeDataToCenterline(volume,centerline,use_stored_values):
	numPts_centerline = centerline.GetNumberOfPoints()
	pressure = dict()
	flow = dict()
	area = dict()
	start = 10200
	stop = 12200
	step = 1000
	num_timesteps = (stop-start+step)/step
	if(not use_stored_values):
		for cPt in tqdm(range(0,numPts_centerline)):
			slice = slice_vessel(volume, centerline.GetPoint(cPt), centerline.GetPointData().GetArray('CenterlineSectionNormal').GetTuple(cPt))
			pressure[cPt],flow[cPt],area[cPt] = calculateFlowPressure(slice,start,stop,step)
		f = open("pressure.pkl","wb")
		pickle.dump(pressure,f)
		f.close()
		f = open("flow.pkl","wb")
		pickle.dump(flow,f)
		f.close()
		f = open("area.pkl","wb")
		pickle.dump(area,f)
		f.close()
	else:
		pressure = pickle.load(open('pressure.pkl','rb'))
		flow = pickle.load(open('flow.pkl','rb'))
		area = pickle.load(open('area.pkl','rb'))
	pressure_avg = dict()
	flow_avg = dict()
	area_avg = dict()
	for p in pressure:
		pressure_avg[p] = 0
		for t in pressure[p]:
			pressure_avg[p] += pressure[p][t]
		pressure_avg[p] = pressure_avg[p]/num_timesteps
	for f in flow:
		flow_avg[f] = 0
		for t in flow[f]:
			flow_avg[f] += flow[f][t]
		flow_avg[f] = flow_avg[f]/num_timesteps
	for a in area:
		area_avg[a] = 0
		for t in area[a]:
			area_avg[a] += area[a][t]
		area_avg[a] = area_avg[a]/num_timesteps
	pressure_array = vtk.vtkDoubleArray()
	pressure_array.SetName('pressure')
	pressure_array.SetNumberOfValues(numPts_centerline)
	pressure_array.Fill(-1)
	for p in pressure_avg:
		pressure_array.SetValue(p,pressure_avg[p])
	flow_array = vtk.vtkDoubleArray()
	flow_array.SetName('flow')
	flow_array.SetNumberOfValues(numPts_centerline)
	flow_array.Fill(-1)
	for f in flow_avg:
		flow_array.SetValue(f,flow_avg[f])
	area_array = vtk.vtkDoubleArray()
	area_array.SetName('area')
	area_array.SetNumberOfValues(numPts_centerline)
	area_array.Fill(-1)
	for a in area_avg:
		area_array.SetValue(a,area_avg[a])
	centerline.GetPointData().AddArray(pressure_array)
	centerline.GetPointData().AddArray(flow_array)
	centerline.GetPointData().AddArray(area_array)
	return centerline

def mapBranch_resistances(branch_resistances, centerline, array_name):
	#Check to see if array name exists on centerline and if not, add it
	if not centerline.GetPointData().HasArray(array_name):
		array = vtk.vtkDoubleArray()
		array.SetName(array_name)
		array.SetNumberOfValues(centerline.GetNumberOfPoints())
		array.Fill(-1)
		centerline.GetPointData().AddArray(array)
		print(array_name,'added')
	#loop through each branch and assign the resistance value to those node Ids
	for b in branch_resistances:
		Branch_PtIds = extractRegion(centerline,b,'BranchId')
		for i in range(0,Branch_PtIds.GetNumberOfIds()):
			centerline.GetPointData().GetArray(array_name).SetValue(Branch_PtIds.GetId(i),branch_resistances[b])
def createParser():
    parser = argparse.ArgumentParser(description='Calculates resistance of a given 3D model based on centerline.')
    parser.add_argument('centerline', type=str, help='the centerline to map diameters from')
    parser.add_argument('volume', type=str, help='the volume to map from')
    parser.add_argument('-use_stored_values', type=int, nargs='?', const=1, default=0, help='turn on use_stored_values')
    parser.add_argument('-o', '-option', type=int, nargs='?', const=1, default=0, help='assign Bifur/branch ID from blanking')
    parser.add_argument('-seg', '-seg_coverage', type=float, nargs='?', default=0.8, help='fraction of the segment to consider for pressure difference')
    parser.add_argument('-t', '-timesteps', type=int, nargs='?', const=1, default=0, help='manually edit timesteps')
    parser.add_argument('-v', '-verbose', type=int, nargs='?', const=1, default=0, help='turn on verbosity')
    return parser

def extractRegionVolume(mesh,selectionCells):
	#Intialize variables
	cellIds = vtk.vtkIdTypeArray()

	#Converts the vtkIdList into vtkIdTypeArray
	for i in range(0,selectionCells.GetNumberOfIds()):
		cellIds.InsertNextValue(selectionCells.GetId(i))
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
	print(region.GetOutput().GetNumberOfCells())
	dssf = vtk.vtkDataSetSurfaceFilter()
	dssf.SetInputConnection(region.GetOutputPort())
	dssf.Update()
	return dssf.GetOutput()

def addBranchBifurcIDs(centerline):
	#remove duplicate points
	cleaner = vtk.vtkCleanPolyData()
	cleaner.SetInputData(centerline)
	cleaner.PointMergingOn()
	cleaner.Update()
	centerline = cleaner.GetOutput()
	branch_regions_IDs = extractCellRegion(centerline,0,'Blanking')
	branch_regions = extractRegionVolume(centerline,branch_regions_IDs)
	bifurc_regions_IDs = extractCellRegion(centerline,1,'Blanking')
	bifurc_regions = extractRegionVolume(centerline,bifurc_regions_IDs)
	threshold = vtk.vtkThreshold()
	threshold.SetInputData(centerline)
	threshold.ThresholdByUpper(1)
	threshold.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, "Blanking")
	threshold.Update()
	#Initialize Branch and Bifurcation Arrays
	array = vtk.vtkDoubleArray()
	array.SetName('BifurcationID')
	array.SetNumberOfValues(centerline.GetNumberOfPoints())
	array.Fill(-1)
	centerline.GetPointData().AddArray(array)
	array = vtk.vtkDoubleArray()
	array.SetName('BranchID')
	array.SetNumberOfValues(centerline.GetNumberOfPoints())
	array.Fill(-1)
	centerline.GetPointData().AddArray(array)
	#loop through each bifurcation and branch connected region and assign_id to those points
	connectivity = vtk.vtkConnectivityFilter();
	connectivity.SetInputData(bifurc_regions);
	connectivity.SetExtractionModeToSpecifiedRegions();
	connectivity.InitializeSpecifiedRegionList()
	connectivity.Update()
	write_polydata(bifurc_regions,'bifurc.vtp')
	print('Assigning Bifurcation Ids')
	for i_region in tqdm(range(0,connectivity.GetNumberOfExtractedRegions())):
		connectivity.AddSpecifiedRegion(i_region)
		connectivity.Update()
		region = connectivity.GetOutput()
		write_polydata(region,'region.vtp')
		for pt in range(0,region.GetNumberOfPoints()):
			centerline.GetPointData().GetArray('BifurcationID').SetValue(centerline.FindPoint(region.GetPoint(pt)),i_region)
		connectivity.DeleteSpecifiedRegion(i_region)
	connectivity = vtk.vtkConnectivityFilter();
	connectivity.SetInputData(branch_regions);
	connectivity.SetExtractionModeToSpecifiedRegions();
	connectivity.InitializeSpecifiedRegionList()
	connectivity.Update()
	write_polydata(branch_regions,'branch.vtp')
	print('Assigning Branch Ids')
	for i_region in tqdm(range(0,connectivity.GetNumberOfExtractedRegions())):
		connectivity.AddSpecifiedRegion(i_region)
		connectivity.Update()
		region = connectivity.GetOutput()
		for pt in range(0,region.GetNumberOfPoints()):
			centerline.GetPointData().GetArray('BranchID').SetValue(centerline.FindPoint(region.GetPoint(pt)),i_region)
	write_polydata(centerline,'test_centerline.vtp')

def calculateTangent(pt1,pt2):
    diff = [0,0,0]
    vtk.vtkMath.Subtract(pt1,pt2,diff)
    return diff

def addCenterlineSectionNormal(centerline):
	#Initialize data arrays
	array = vtk.vtkDoubleArray()
	array.SetName('CenterlineSectionNormal')
	array.SetNumberOfComponents(3)
	array.SetNumberOfTuples(centerline.GetNumberOfPoints())
	array.Fill(-1)
	centerline.GetPointData().AddArray(array)

	for i_cell in range(0,centerline.GetNumberOfCells()):
		pt_list = vtk.vtkIdList()
		centerline.GetCellPoints(i_cell,pt_list)
		for pt in range(0,pt_list.GetNumberOfIds()-1):
			p1 = centerline.GetPoint(pt_list.GetId(pt))
			p2 = centerline.GetPoint(pt_list.GetId(pt+1))
			centerline.GetPointData().GetArray('CenterlineSectionNormal').SetTuple(pt_list.GetId(pt),calculateTangent(p2,p1))
		if(pt_list.GetNumberOfIds()>1):
			p0 = centerline.GetPoint(pt_list.GetId(pt_list.GetNumberOfIds()-2))
			p1 = centerline.GetPoint(pt_list.GetId(pt_list.GetNumberOfIds()-1))
		centerline.GetPointData().GetArray('CenterlineSectionNormal').SetTuple(pt_list.GetId(pt_list.GetNumberOfIds()-1),calculateTangent(p1,p0))
	return centerline

#Adds two arrays (BifurcationID and BranchID) to the point data of the centerline from GroupIDs and BlankingID (0:Branch, 1:Bifurc)
def assignBranchBifurcIDs(centerline):
	#remove duplicate points
	cleaner = vtk.vtkCleanPolyData()
	cleaner.SetInputData(centerline)
	cleaner.PointMergingOn()
	cleaner.Update()
	centerline = cleaner.GetOutput()
	#Initialize data arrays
	array = vtk.vtkDoubleArray()
	array.SetName('BranchId')
	array.SetNumberOfValues(centerline.GetNumberOfPoints())
	array.Fill(-1)
	centerline.GetPointData().AddArray(array)
	
	array = vtk.vtkDoubleArray()
	array.SetName('BifurcationId')
	array.SetNumberOfValues(centerline.GetNumberOfPoints())
	array.Fill(-1)
	centerline.GetPointData().AddArray(array)

	#loop through each GroupId
	branchId_counter = 0
	bifurcationId_counter = 0
	[GroupId_min, GroupId_max] = centerline.GetCellData().GetArray('GroupIds').GetRange()
	for gId in range(int(GroupId_min),int(GroupId_max+1)):
		#extract the region of the centerline with that GroupId and loop through each cell of that region
		region = extractCellRegion(centerline,gId,'GroupIds')
		if region.GetNumberOfIds()>0:
			for i_cell in range(0,region.GetNumberOfIds()):
				#if blanking==1 for that cell -> assign points of the cell to bifurcationIds data array with a value of the bifurcationId_counter
				if (centerline.GetCellData().GetArray('Blanking').GetValue(region.GetId(i_cell)) == 1):
					pt_list = vtk.vtkIdList()
					centerline.GetCellPoints(region.GetId(i_cell),pt_list)
					for ipt in range(0,pt_list.GetNumberOfIds()):
						centerline.GetPointData().GetArray('BifurcationId').SetValue(pt_list.GetId(ipt),bifurcationId_counter)
				#else assign points of the cell to branchId array with a value of the branchId_counter
				else:
					pt_list = vtk.vtkIdList()
					centerline.GetCellPoints(region.GetId(i_cell),pt_list)
					for ipt in range(0,pt_list.GetNumberOfIds()):
						centerline.GetPointData().GetArray('BranchId').SetValue(pt_list.GetId(ipt),branchId_counter)
			if (centerline.GetCellData().GetArray('Blanking').GetValue(region.GetId(0)) == 1):
				bifurcationId_counter += 1
			else:
				branchId_counter += 1
	write_polydata(centerline,'test_centerline.vtp')
	return centerline

def main(args):
	#intialize
	reader = vtk.vtkXMLUnstructuredGridReader()
	reader.SetFileName(args.volume)
	reader.Update()
	volume = reader
	centerline = read_polydata(args.centerline)
	SEG_COVERAGE = args.seg
	global vprint
	if args.v:
		def vprint(*args):
			# Print each argument separately so caller doesn't need to
			# stuff everything to be printed into a single string
			for arg in args:
				print(arg),	
			print('\n')
	else:
		vprint = lambda *a: None
	if(args.o==1):
		centerline = addCenterlineSectionNormal(centerline)
		centerline = assignBranchBifurcIDs(centerline)
	start = 10200
	stop =12200
	step = 400
	#Map flow and pressure data from volume mesh to centerline
	#centerline = mapVolumeDataToCenterline(volume,centerline,args.use_stored_values)
	#write_polydata(centerline,'mapped_centerline.vtp')
	#for every branch Id this dict contains a list of the child branch ids
	print('Calculating connectivity...')
	connectivity = defaultdict(list)
	connectivity = {0:[1,2],1:[4,5],2:[6,7],6:[8,9,10]}
	connectivity = getBranchConnectivity(centerline)
	f = open("connectivity.pkl","wb")
	pickle.dump(connectivity,f)
	f.close()
	print(connectivity)

	#for every branch Id this dict contains its resistance
	print('Calculating resistance...')
	branch_resistances = dict()
	branch_resistances = {0:2,1:4,2:6,3:8,4:10,5:12,6:14,7:16,8:18,9:20,10:22}
	if(args.o==0):
		print('Not using sphere cut, time average pressure and flow')
		branch_resistances = getBranchResistancesAverage(centerline,volume,connectivity,start,stop,step,SEG_COVERAGE)
	elif(args.o==1 and args.t==0):
		print('Using sphere cut, using just 1 pressure and flow')
		branch_resistances = getBranchResistances2(centerline,volume,connectivity,SEG_COVERAGE)
	else:
		print('Using sphere cut, time average pressure and flow')
		branch_resistances = getBranchResistances2Average(centerline,volume,connectivity,start,stop,step,SEG_COVERAGE)
	f = open("resistance.pkl","wb")
	pickle.dump(branch_resistances,f)
	f.close()
	#starting branch
	seed_branch = 0
	#the current parent branch
	parent_branch = 0

	global resistance
	resistance = 0	
	resisance = getResistance(parent_branch,centerline,connectivity,branch_resistances)
	print(resistance)
	write_polydata(centerline,'mapped_centerline.vtp')

if __name__ == '__main__':
    parser = createParser()
    args = parser.parse_args()
    main(args)

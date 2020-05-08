from collections import defaultdict
import vtk
import vtk_functions
import tqdm

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

def calculateFlowPressure(slice):
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
	numPts = model.GetPointData().GetNumberOfPoints()
	for i in range(0,numPts):
		if value == array.GetValue(i):
			regionIds.InsertUniqueId(i)
	return regionIds

def getBranchConnectivity(model):
	connectivity = defaultdict(list)
	[BranchId_min, BranchId_max] = model.GetPointData().GetArray('BifurcationId').GetRange()
	Bifurcation_PtIds = vtk.vtkIdList()
	#loop through each bifurcationId
	for b in range(BranchId_min,BranchId_max+1):
		Bifurcation_PtIds = extractRegion(model,b,'BifurcationId')
		num_different_branchIds = set()
		#loop through each point in bifurcation region
		for pt in range(0,Bifurcation_PtIds.GetNumberOfIds()):
			#find all neighbor points for each bifurcation point and see if it's connected to a different branchId
			neighbor_pts = getConnectedVerticesNotIncludingSeed(model, pt)
			for neighbor in range(0,neighbor_pts.GetNumberOfIds()):
				neighbor_branchId = model.GetPointData().GetArray('BranchId').GetValue(neighbor_pts.GetId(neighbor))
				if neighbor_branchId >= 0 and neighbor_branchId not in num_different_branchIds:
					num_different_branchIds.add(neighbor_branchId)
		#Assumes parent branch has the lowest BranchId
		parent_branch = 0
		child_branches = []
		for branchId in num_different_branchIds:
			if branchId < parent_branch:
				parent_branch = branchId
			else:
				child_branches.append(branchId)

		connectivity[parent_branch] = child_branches
	return connectivity

def getResistance(branch,connectivity,branch_resistances):
	global resistance
	global parallel_res
	#checks to see if parent branch is in connectivity dictionary
	if branch in connectivity.keys():
		print(branch)
		#If yes, then loop through each child branch and recursively call the method for downstream branches summed in parallel
		print(connectivity[branch])
		parallel_res = 0
		for child_branch in connectivity[branch]:
			parallel_res += 1.0/getResistance(child_branch,connectivity,branch_resistances)
			print('p_res',parallel_res)

		resistance = 1.0/parallel_res + branch_resistances[branch]
	else:
		#Otherwise, return the resistance of the branch itself
		return branch_resistances[branch]
	return resistance

def getBranchResistances(model):
	branch_resistances = dict()
	#loop through each branchId
	[BranchId_min, BranchId_max] = model.GetPointData().GetArray('BifurcationId').GetRange()
	for branchId in range(BranchId_min,BranchId_max+1):
		Branch_PtIds = extractRegion(model,branchId,'BranchId')
		#loop through each branch point and store pressure at each node
		pressures = dict()
		for pt in range(0,Branch_PtIds.GetNumberOfIds()):
			pressures[pt] = model.GetPointData().GetArray('pressure').GetValue(pt)
		#loop through pressures to find min and max pressure
		pmax = 0
		pmax_node = 0
		pmin = 0
		pmin_node = 0
		for p in pressures:
			if pressures[p] > pmax:
				pmax = pressures[p]
				pmax_node = p
			if pressures[p] < pmin:
				pmin = pressures[p]
				pmin_node = p
		#Calculate pressure difference and flow at the minimum pressure node to find resistance of branch
		branch_pdiff = pmax - pmin
		branch_flow = model.GetPointData().GetArray('flow').GetValue(pmin_node)
		branch_resistances[branchId] = branch_pdiff/branch_flow
	return branch_resistances

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

    return con

def mapVolumeDataToCenterline(volume,centerline):
	numPts_centerline = centerline.GetNumberOfPoints()
	pressure = dict()
	flow = dict()
	area = dict()
	for cPt in tqdm(range(0,numPts_centerline)):
		slice = slice_vessel(volume, cPt, centerline.GetPointData().GetArray('CenterlineSectionNormal').GetTuple(cPt))
		pressure[i],flow[i],area[i] = calculateFlowPressure(slice)
	pressure_array = vtk.vtkDoubleArray()
	pressure_array.SetName('pressure')
	for p in presure:
		pressure_array.SetValue(p,pressure[p])
	flow_array = vtk.vtkDoubleArray()
	flow_array.SetName('flow')
	for f in flow:
		flow_array.SetValue(f,flow[f])
	area_array = vtk.vtkDoubleArray()
	area_array.SetName('area')
	for a in area:
		pressure.SetValue(a,area[a])
	centerline.GetPointData().AddArray(pressure_array)
	centerline.GetPointData().AddArray(flow_array)
	centerline.GetPointData().AddArray(area_array)
	return centerline

def createParser():
    parser = argparse.ArgumentParser(description='Calculates resistance of a given 3D model based on centerline.')
    parser.add_argument('centerline', type=str, help='the centerline to map diameters from')
    parser.add_argument('volume', type=str, help='the volume to map from')
    parser.add_argument('-v', '-verbose', type=int, nargs='?', const=1, default=0, help='turn on verbosity')
    return parser

def main(args):
	#intialize
	reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(args.volume)
    reader.Update()
	volume = reader
	centerline = read_polydata(args.centerline)
	#Map flow and pressure data from volume mesh to centerline
	centerline = mapVolumeDataToCenterline(volume,centerline)
	write_polydata(centerline,'mapped_centerline.vtp')
	#for every branch Id this dict contains a list of the child branch ids
	connectivity = defaultdict(list)
	connectivity = {0:[1,2],1:[4,5],2:[6,7],6:[8,9,10]}
	connectivity = getBranchConnectivity(centerline)
	#for every branch Id this dict contains its resistance
	branch_resistances = dict()
	branch_resistances = {0:2,1:4,2:6,3:8,4:10,5:12,6:14,7:16,8:18,9:20,10:22}
	branch_resistances = getBranchResistances(centerline)
	#starting branch
	seed_branch = 0
	#the current parent branch
	parent_branch = 0

	global resistance
	resistance = 0	
	resisance = getResistance(parent_branch, connectivity,branch_resistances)
	print(resistance)

if __name__ == '__main__':
    parser = createParser()
    args = parser.parse_args()
    main(args)

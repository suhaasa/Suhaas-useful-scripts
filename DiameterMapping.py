import sys
import os
import vtk
import numpy as np
import time
import glob
import math
import argparse
from collections import defaultdict
from tqdm import tqdm
import pickle
import math

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

def calcDistance2Points(model, pt1,pt2):
    if(type(pt1) is int or type(pt1) is long):
        x1,y1,z1 = model.GetPoint(pt1)
    elif(type(pt1) is list):
        x1,y1,z1 = pt1[0],pt1[1],pt1[2]
    else:
        vprint(type(pt1))
    if(type(pt2) is int or type(pt2) is long):
        x2,y2,z2 = model.GetPoint(pt2)
    else:
        x2,y2,z2 = pt2[0],pt2[1],pt2[2]
    distance = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**(.5)
    return distance

#Graph class
class Graph:
    def __init__(self):
        self.nodes = set()
        self.edges = defaultdict(list)
        self.node_coord = dict()
        self.node_properties = dict()
        self.distances = {}

    def add_node(self, value):
        self.nodes.add(value)

    def add_node_coord(self, node, coord):
        self.node_coord[node] = coord

    def add_node_property(self,node,value):
        self.node_properties[node] = value

    def add_edge(self, from_node, to_node, distance):
        self.edges[from_node].append(to_node)
        self.edges[to_node].append(from_node)
        self.distances[(from_node, to_node)] = distance
        self.distances[(to_node, from_node)] = distance

    #adds a virtual node that collapses a group of nodes to one virtual node to make a source and sets distance to 0
    def add_virtual_node(self, v_node, zero_edge_nodes):
        self.nodes.add(v_node)
        self.node_coord[v_node] = [0,0,0]
        for i in zero_edge_nodes:
            self.edges[v_node].append(i)
            self.edges[i].append(v_node)
            self.distances[(v_node, i)] = 0
            self.distances[(i, v_node)] = 0
    
    #distance cutoff is a value that ensures RCA centerlines don't wrap onto the LV myocardium through the right ventricle
    def add_virtual_node_distances(self, v_node, edge_nodes, distances):
        self.nodes.add(v_node)
        self.node_coord[v_node] = [0,0,0]
        DISTANCE_CUTOFF = 1
        for i in edge_nodes:
            if(distances[i]<DISTANCE_CUTOFF):
                self.edges[v_node].append(i)
                self.edges[i].append(v_node)
                self.distances[(v_node, i)] = distances[i]
                self.distances[(i, v_node)] = distances[i]

    def get_node_coord(self,node):
        return self.node_coord[node]

    def get_num_of_nodes(self):
        return len(self.nodes)

#generates a graph of the mesh 
def generateGraph(mesh):
    print('Generating graph...')
    graph = Graph()
    print(mesh.GetNumberOfPoints())
    for i in tqdm(range(0,mesh.GetNumberOfPoints())):
        graph.add_node(i)
        graph.add_node_coord(i,mesh.GetPoint(i))
        connnectedPt_list = getConnectedVerticesNotIncludingSeed(mesh,i)
        for j in range(0,connnectedPt_list.GetNumberOfIds()):
            # new point to decide whether to add to patch, edge, or nothing (if already in edge)
            cpt = connnectedPt_list.GetId(j)
            graph.add_edge(i,cpt,calcDistance2Points(mesh,i,cpt))
    return graph

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

def dijsktra(graph, initial):
    visited = {}
    visited[initial] = 0
    path = {}
    path_nodes = set()

    nodes = set(graph.nodes)
    counter = 0
    pbar = tqdm(total=len(nodes))
    while nodes: 
        pbar.update(1)
        min_node = None
        for node in nodes:
            if node in visited:
                if min_node is None:
                    min_node = node
                elif visited[node] <= visited[min_node]:
                    min_node = node
        if min_node is None:
            break

        nodes.remove(min_node)
        current_weight = visited[min_node]

        for edge in graph.edges[min_node]:
            weight = current_weight + graph.distances[(min_node, edge)]
            if edge not in visited or weight <= visited[edge]:
                visited[edge] = weight
                path[edge] = min_node
                path_nodes.add(edge)
        counter += 1
    pbar.close()
    return visited, path

def dijsktra_expand_properties(graph, initial, properties):
    visited = {}
    visited[initial] = 0
    path = {}
    path_nodes = set()

    nodes = set(graph.nodes)
    counter = 0
    pbar = tqdm(total=len(nodes))
    while nodes: 
        pbar.update(1)
        min_node = None
        for node in nodes:
            if node in visited:
                if min_node is None:
                    min_node = node
                elif visited[node] <= visited[min_node]:
                    min_node = node
        if min_node is None:
            break

        nodes.remove(min_node)
        current_weight = visited[min_node]

        for edge in graph.edges[min_node]:
            if min_node in properties:
                prop_dict = properties[min_node]
                diff = [0,0,0]
                vtk.vtkMath.Subtract(graph.get_node_coord(edge), graph.get_node_coord(min_node), diff)
                multiplier = vtk.vtkMath.Dot(prop_dict['Tangent'],diff)
                #print('multiplier set to '+str(multiplier))
            multiplier = 1
            #print('multiplier set to 1')
            weight = current_weight + graph.distances[(min_node, edge)]*abs(multiplier)
            if edge not in visited or weight <= visited[edge]:
                visited[edge] = weight
                path[edge] = min_node
                path_nodes.add(edge)
                if min_node in properties:
                    properties[edge] = properties[min_node]
                    graph.add_node_property(edge,properties[min_node])
                    graph.add_node_property(min_node,properties[min_node])
        counter += 1
    pbar.close()
    return visited, path

def multipleSourceDistance(mesh,graph,v_node,child_nodes,distances,properties):
    graph.add_virtual_node(v_node,child_nodes)
    f = open("properties.pkl","wb")
    pickle.dump(properties,f)
    f.close()
    visited,path = dijsktra_expand_properties(graph,v_node,properties)
    data_array = vtk.vtkDoubleArray()
    data_array.SetName('distance_map')
    for i in range(0,mesh.GetNumberOfPoints()):
        if i in visited:
            data_array.InsertNextValue(visited[i])
        else:
            data_array.InsertNextValue(-1)
    mesh.GetPointData().AddArray(data_array)
    return mesh

def BranchSurface(surface, centerlines, RadiusArrayName, nameThis, nameOther):
    #output id array
    thisSurf = vtk.vtkDoubleArray();
    thisSurf.SetName(nameThis);
    thisSurf.SetNumberOfValues(surface.GetNumberOfPoints());
    thisSurf.Fill(-1);
    

    
    #output distance array
    surfDist = vtk.vtkDoubleArray();
    surfDist.SetNumberOfValues(surface.GetNumberOfPoints());
    surfDist.Fill(-1);

    #build locator for surface points
    dataset = vtk.vtkPolyData();
    dataset.SetPoints(surface.GetPoints());

    locator = vtk.vtkPointLocator();
    locator.Initialize();
    locator.SetDataSet(dataset);
    locator.BuildLocator();

    #get arrays
    centRadius = centerlines.GetPointData().GetArray(RadiusArrayName);
    normals = surface.GetPointData().GetArray("Normals");
    otherCent = centerlines.GetPointData().GetArray(nameOther);
    thisCent = centerlines.GetPointData().GetArray(nameThis);

    cellIds = vtk.vtkIdList();
    surfPointIds = vtk.vtkIdList();
    p_cent = [0,0,0]
    p_surf = [0,0,0]
    normal = [0,0,0]
    diff = [0,0,0]

    for i in tqdm(xrange(0,surface.GetNumberOfPoints())):
        #thisId = thisCent.GetTuple1(i);
        #otherId = otherCent.GetTuple1(i);

        #skip bifurcation points
        #if (otherId != -1):
        #    continue;

        #pre-select surface points within sphere (very fast)
        surface.GetPoint(i, p_cent);
        #locator.FindPointsWithinRadius(10.0 * radius, p_cent, surfPointIds);
        closestId = locator.FindClosestPoint(p_cent);
        vtk.vtkMath.Subtract(centerlines.GetPoint(closestId), p_cent, diff);
        radius = centRadius.GetValue(closestId);
        dist = vtk.vtkMath.Norm(diff);
        thisSurf.SetValue(surface.FindPoint(p_cent), radius);
        surfDist.SetValue(surface.FindPoint(p_cent), dist);
        #select surface points according to distances (slow)
        # for j in xrange(0,surfPointIds.GetNumberOfIds()):
        #     #get surface point and normal
        #     surface.GetPoint(surfPointIds.GetId(j), p_surf);
        #     #normals.GetTuple(surfPointIds.GetId(j), normal);

        #     #distance vector between surface and centerline points
        #     vtk.vtkMath.Subtract(p_surf, p_cent, diff);
        #     dist = vtk.vtkMath.Norm(diff);

        #     #signed distance in surface normal direction
        #     #dist_n = vtk.vtkMath.Dot(diff, normal);

        #     #check if centerline is inside branch (allow small tolerance for caps)
        #     #if -1.0e-2 <= dist_n:
        #     #if surface point already has an id from closer centerline, skip
        #     if (-1 < thisSurf.GetValue(surfPointIds.GetId(j))):
        #         if (dist > surfDist.GetValue(surfPointIds.GetId(j))):
        #            continue;

        #         #set BranchId and distance
        #         thisSurf.SetValue(surfPointIds.GetId(j), radius);
        #         surfDist.SetValue(surfPointIds.GetId(j), dist);
    surface.GetPointData().AddArray(thisSurf);
    surface.GetPointData().AddArray(surfDist);
    return surface

def addPropertiesFromGraph(mesh,graph):
    for i in graph.nodes:
        if i in graph.node_properties:
            node_property = graph.node_properties[i]
            for property_dict in node_property:
                for array_name in property_dict:
                    if not mesh.GetPointData().HasArray(array_name):
                        data = vtk.vtkDoubleArray()
                        data.SetName(array_name)
                        data.SetNumberOfValues(mesh.GetNumberOfPoints())
                        data.Fill(-1)
                        mesh.GetPointData().AddArray(data)
                        print(array_name + ' data array added to mesh.')
                    else:
                        mesh.GetPointData().GetArray(array_name).SetValue(i,property_dict[array_name])
    return mesh

def addPropertiesFromDict(mesh,dict):
    for p in range(0,mesh.GetNumberOfPoints()):
        node = mesh.GetPointData().GetArray('GlobalNodeID').GetValue(p)
        if node >= 0 and node < mesh.GetNumberOfPoints() and node in dict:
            property_dict = dict[node]
            for array_name in property_dict:
                if not mesh.GetPointData().HasArray(array_name):
                    data = vtk.vtkDoubleArray()
                    data.SetName(array_name)
                    if(type(property_dict[array_name]) is float or type(property_dict[array_name]) is int):
                        data.SetNumberOfComponents(1)
                    else:
                        data.SetNumberOfComponents(len(property_dict[array_name]))
                    data.SetNumberOfTuples(mesh.GetNumberOfPoints())
                    data.Fill(-1)
                    mesh.GetPointData().AddArray(data)
                    print(array_name + ' data array added to mesh.')
                    print(data.GetNumberOfTuples())
                if(type(property_dict[array_name]) is float or type(property_dict[array_name]) is int):
                    mesh.GetPointData().GetArray(array_name).SetValue(p,property_dict[array_name])
                else:
                    mesh.GetPointData().GetArray(array_name).SetTuple(p,property_dict[array_name])
    return mesh

def calculateTangent(pt1,pt2):
    diff = [0,0,0]
    vtk.vtkMath.Subtract(pt1,pt2,diff)
    return diff

def getCenterlineTangents(centerline):
    tangents = vtk.vtkDoubleArray()
    tangents.SetNumberOfComponents(3)
    tangents.SetNumberOfTuples(centerline.GetNumberOfPoints())
    tangents.SetName('Tangent')
    tangents.Fill(-1)
    connectedPts_list = vtk.vtkIdList()
    for i_cell in xrange(0,centerline.GetNumberOfCells()):
        pt_list = vtk.vtkIdList()
        pt_list = centerline.GetCell(i_cell).GetPointIds()
        ptIDs = []
        for k in range(0,pt_list.GetNumberOfIds()):
            ptIDs.append(pt_list.GetId(k))
        ptIDs.sort()
        if len(ptIDs)>1:
            tangents.SetTuple(ptIDs[0],calculateTangent(centerline.GetPoint(ptIDs[0]),centerline.GetPoint(ptIDs[1])))
            for ipt in xrange(1,len(ptIDs)):
                tangents.SetTuple(ptIDs[ipt-1],calculateTangent(centerline.GetPoint(ptIDs[ipt-1]),centerline.GetPoint(ptIDs[ipt])))
    centerline.GetPointData().AddArray(tangents)

def getWeights(pt,mesh,NclosestPts):
    weights = dict()
    sum = 0
    diff = [0,0,0]
    for i in range(0,NclosestPts.GetNumberOfIds()):
        vtk.vtkMath.Subtract(pt,mesh.GetPoint(NclosestPts.GetId(i)),diff)
        weights[NclosestPts.GetId(i)] = vtk.vtkMath.Norm(diff)
        sum += weights[NclosestPts.GetId(i)]
    for i in weights:
        weights[i] = weights[i]/sum
    return weights

def mapToNClosestCenterline(centerline,mesh,n):
    for array in range(0,centerline.GetPointData().GetNumberOfArrays()): 
        data = vtk.vtkDoubleArray()
        data.SetName(centerline.GetPointData().GetArray(array).GetName())
        data.SetNumberOfComponents(centerline.GetPointData().GetArray(array).GetNumberOfComponents())
        data.SetNumberOfTuples(mesh.GetNumberOfPoints())
        data.Fill(-1)
        mesh.GetPointData().AddArray(data)
        print(centerline.GetPointData().GetArray(array).GetName() + ' data array added to mesh.')
    properties = dict()
    for i in tqdm(range(0,mesh.GetNumberOfPoints())):
        pt = mesh.GetPoint(i)
        NclosestPts = vtk.vtkIdList()
        #Make locator object of centerline points
        dataset = vtk.vtkPolyData();
        dataset.SetPoints(centerline.GetPoints());

        locator = vtk.vtkPointLocator();
        locator.Initialize();
        locator.SetDataSet(dataset);
        locator.BuildLocator();
        locator.FindClosestNPoints(n,pt,NclosestPts);
        property_dict = dict()
        #make weights dictionary of ids and distance norms to closest points
        weights = getWeights(pt,mesh,NclosestPts)
        for array in range(0,centerline.GetPointData().GetNumberOfArrays()):
            summed_weights = np.zeros(centerline.GetPointData().GetArray(array).GetNumberOfComponents())
            for w in weights:
                summed_weights += [weights[w]*value for value in centerline.GetPointData().GetArray(array).GetTuple(w)] 
            mesh.GetPointData().GetArray(centerline.GetPointData().GetArray(array).GetName()).SetTuple(i,summed_weights)
            property_dict[str(centerline.GetPointData().GetArray(array).GetName())] = summed_weights                       
        properties[mesh.FindPoint(centerline.GetPoint(i))] = property_dict
    return mesh, properties

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

def smoothBifurcations(centerline,mesh):
    bifurcation_regionIds = extractCellRegion(centerline,1,'Blanking')
    dataset = vtk.vtkPolyData();
    dataset.SetPoints(mesh.GetPoints());
    locator = vtk.vtkPointLocator();
    locator.Initialize();
    locator.SetDataSet(dataset);
    locator.BuildLocator();
    for i_cell in tqdm(range(0,bifurcation_regionIds.GetNumberOfIds())):
        cellpoints = vtk.vtkIdList()
        centerline.GetCellPoints(bifurcation_regionIds.GetId(i_cell),cellpoints)
        for p in range(0,cellpoints.GetNumberOfIds()):
            ipt = bifurcation_regionIds.GetId(p)
            closePoints = vtk.vtkIdList()
            radius = centerline.GetPointData().GetArray('MaximumInscribedSphereRadius').GetValue(ipt)
            locator.FindPointsWithinRadius(radius*1, centerline.GetPoint(ipt), closePoints);
            for cpt in range(0,closePoints.GetNumberOfIds()):
                if mesh.GetPointData().GetArray('MaximumInscribedSphereRadius').GetValue(closePoints.GetId(cpt))<radius:
                    mesh.GetPointData().GetArray('MaximumInscribedSphereRadius').SetValue(closePoints.GetId(cpt),radius)

    return mesh

def smoothBifurcations2(centerline,mesh):
    bifurcation_regionIds = extractRegion(centerline,1,'CenterlineSectionBifurcation')
    dataset = vtk.vtkPolyData();
    dataset.SetPoints(mesh.GetPoints());
    locator = vtk.vtkPointLocator();
    locator.Initialize();
    locator.SetDataSet(dataset);
    locator.BuildLocator();
    for p in tqdm(range(0,bifurcation_regionIds.GetNumberOfIds())):
        ipt = bifurcation_regionIds.GetId(p)
        closePoints = vtk.vtkIdList()
        radius = centerline.GetPointData().GetArray('MaximumInscribedSphereRadius').GetValue(ipt)
        locator.FindPointsWithinRadius(radius*1.5, centerline.GetPoint(ipt), closePoints);
        for cpt in range(0,closePoints.GetNumberOfIds()):
            if mesh.GetPointData().GetArray('MaximumInscribedSphereRadius').GetValue(closePoints.GetId(cpt))<radius:
                mesh.GetPointData().GetArray('MaximumInscribedSphereRadius').SetValue(closePoints.GetId(cpt),radius)
                    
    return mesh

def connectivity(inp, origin):
    """
    If there are more than one unconnected geometries, extract the closest one
    Args:
        inp: InputConnection
        origin: region closest to this point will be extracted
    Returns:
        con: connectivity object
    """
    con = vtk.vtkConnectivityFilter()
    con.SetInputData(inp.GetOutput())
    con.SetExtractionModeToClosestPointRegion()
    con.SetClosestPoint(origin[0], origin[1], origin[2])
    con.Update()
    return con

def extract_surface(inp):
    """
    Extract surface from 3D geometry
    Args:
        inp: InputConnection
    Returns:
        extr: vtkExtractSurface object
    """
    extr = vtk.vtkDataSetSurfaceFilter()
    extr.SetInputData(inp.GetOutput())
    extr.Update()
    return extr

def maxDistanceBetweenPoints(model, seedPt, connectedPts_list):
    max = 0
    for i in xrange(0,connectedPts_list.GetNumberOfIds()):
        distance = calcDistance2Points(model, seedPt,connectedPts_list.GetId(i))
        if(distance > max):
            max = distance
    return max

def getConnectedVerticesWithinRadius(model, seedPt, radius):
    cell_list = vtk.vtkIdList()
    connectedPts_list = vtk.vtkIdList()
    model.GetPointCells(seedPt,cell_list)
    radiusReached = 0
    max_iter = 500
    prev_iter = -1
    while (radiusReached == 0 and connectedPts_list.GetNumberOfIds()<max_iter and connectedPts_list.GetNumberOfIds()>prev_iter):
        prev_iter = connectedPts_list.GetNumberOfIds()
        for j in xrange(0,cell_list.GetNumberOfIds()):
            pt_list = vtk.vtkIdList()
            pt_list = model.GetCell(cell_list.GetId(j)).GetPointIds()
            for k in xrange(0,pt_list.GetNumberOfIds()):
                connectedPts_list.InsertUniqueId(pt_list.GetId(k))
        if (maxDistanceBetweenPoints(model, seedPt,connectedPts_list) > radius):
            radiusReached = 1
        else:
            connectedCell_list = vtk.vtkIdList()
            for i in xrange(0,pt_list.GetNumberOfIds()):
                model.GetPointCells(pt_list.GetId(i),connectedCell_list)
                for j in xrange(0,connectedCell_list.GetNumberOfIds()):
                    cell_list.InsertUniqueId(connectedCell_list.GetId(j))




    return connectedPts_list

def getConnectedVerticesNotIncludingSeed(model, seedPt):
    cell_list = vtk.vtkIdList()
    connectedPts_list = vtk.vtkIdList()
    model.GetPointCells(seedPt,cell_list)
    for j in xrange(0,cell_list.GetNumberOfIds()):
        pt_list = vtk.vtkIdList()
        pt_list = model.GetCell(cell_list.GetId(j)).GetPointIds()
        for k in xrange(0,pt_list.GetNumberOfIds()):
            if (pt_list.GetId(k) != seedPt):
                connectedPts_list.InsertUniqueId(pt_list.GetId(k))
    return connectedPts_list

def smoothAll(mesh):
    numPts = mesh.GetNumberOfPoints()
    data = vtk.vtkDoubleArray()
    data.SetName('Averaged_Radius')
    data.SetNumberOfValues(mesh.GetNumberOfPoints())
    data.Fill(-1)
    mesh.GetPointData().AddArray(data)

    dataset = vtk.vtkPolyData();
    dataset.SetPoints(mesh.GetPoints());
    locator = vtk.vtkPointLocator();
    locator.Initialize();
    locator.SetDataSet(dataset);
    locator.BuildLocator();

    for p in tqdm(range(0,numPts)):
        radius = mesh.GetPointData().GetArray('MaximumInscribedSphereRadius').GetValue(p)
        if radius > .5:
            radius = .5
        closePoints = getConnectedVerticesNotIncludingSeed(mesh, p)
        #closePoints = getConnectedVerticesWithinRadius(mesh, p, 1.5*radius)
        #locator.FindPointsWithinRadius(radius*1.5, mesh.GetPoint(p), closePoints)
        radius_avg = 0
        for cpt in range(0,closePoints.GetNumberOfIds()):
            radius_avg += mesh.GetPointData().GetArray('MaximumInscribedSphereRadius').GetValue(closePoints.GetId(cpt))
        radius_avg = radius_avg/closePoints.GetNumberOfIds()
        mesh.GetPointData().GetArray('Averaged_Radius').SetValue(p,radius_avg)
    return mesh

def AssignWallIds(wall_dir,volume):
    walls = []
    for file in os.listdir(wall_dir):
        walls.append[read_polydata(file)]
    wallId = dict()
    for i_wall in range(0,len(walls)):
        numPts_wall = walls[i_wall].GetNumberOfPoints()
        for p in range(0,numPts_wall):
            gId = walls[i_wall].GetPointData('GlobalNodeID').GetValue(p)
            if wallId[gId]==-1:
                wallId[gId] = i_wall
            else:
                wallId[gId] = -1

    numPts = volume.GetNumberOfPoints()
    data = vtk.vtkDoubleArray()
    data.SetName('WallIds')
    data.SetNumberOfValues(volume.GetNumberOfPoints())
    data.Fill(-1)
    volume.GetPointData().AddArray(data)

    numPts_volume = volume.GetNumberOfPoints()
    for p in range(0,numPts_volume):
        gId = volume.GetPointData().GetArray('GlobalNodeID').GetValue(p)
        volume.GetPointData().GetArray('WallIds').SetValue(p,wallId[gId])
    return volume

def generateSurfaceNormals(surface):
    surfaceNormals = vtk.vtkPolyDataNormals();
    surfaceNormals.SetInputData(surface);
    surfaceNormals.SplittingOff();
    surfaceNormals.AutoOrientNormalsOn();
    surfaceNormals.ComputePointNormalsOn();
    surfaceNormals.ConsistencyOn();
    surfaceNormals.Update();

    surface.DeepCopy(surfaceNormals.GetOutput());
    return surface

def checkSurfaceInterior(surface,surface_point,test_point):
    diff = [0,0,0]
    normal = [0,0,0]
    p_surf = [0,0,0]
    #get surface point and normal
    surface.GetPoint(surface_point, p_surf);
    normals = surface.GetPointData().GetArray('Normals')
    normals.GetTuple(surface_point, normal);

    #distance vector between surface and centerline points
    vtk.vtkMath.Subtract(p_surf, test_point, diff);
    dist = vtk.vtkMath.Norm(diff);

    #signed distance in surface normal direction
    dist_n = vtk.vtkMath.Dot(diff, normal);
    return dist_n
    #check if centerline is inside branch (allow small tolerance for caps)
    # if (-1.0e-2 <= dist_n):
    #     return True
    # else:
    #     return False

def AssignWallIds(wall_dir,volume,centerline):
    dataset = vtk.vtkPolyData();
    dataset.SetPoints(centerline.GetPoints());
    locator = vtk.vtkPointLocator();
    locator.Initialize();
    locator.SetDataSet(dataset);
    locator.BuildLocator();
    properties = dict()
    walls = []
    for file in os.listdir(wall_dir):
        if file.startswith('wall_'):
            walls.append(os.path.join(wall_dir,file))
    wallId = dict()

    filler = vtk.vtkFillHolesFilter()

    for i_wall in tqdm(range(0,len(walls))):
        wall = read_polydata(walls[i_wall])
        filler.SetInputData(wall)
        filler.Update()
        write_polydata(filler.GetOutput(),'test_filled_wall.vtp')
        wall = generateSurfaceNormals(wall)
        numPts_wall = wall.GetNumberOfPoints()
        centerline_pts = vtk.vtkIdList()
        for p in tqdm(range(0,numPts_wall)):
            surface_point = wall.GetPoint(p)
            radius = centerline.GetPointData().GetArray('MaximumInscribedSphereRadius').GetValue(locator.FindClosestPoint(surface_point))
            norm_pts = []
            seed_pts = vtk.vtkIdList()
            seed_pts.InsertNextId(locator.FindClosestPoint(surface_point))
            prev_numIds = -1
            multiplier = 2
            max_norm = checkSurfaceInterior(wall,p,centerline.GetPoint(seed_pts.GetId(0)))
            max_norm_pt = seed_pts.GetId(0)
            max_iter = 1
            curr_iter =0
            while len(norm_pts)==0 and seed_pts.GetNumberOfIds()!=prev_numIds and curr_iter<max_iter:
                prev_numIds = seed_pts.GetNumberOfIds()
                #locator.FindPointsWithinRadius(multiplier*radius,centerline.GetPoint(locator.FindClosestPoint(surface_point)),centerline_pts)
                locator.FindPointsWithinRadius(multiplier*radius,surface_point,centerline_pts)
                for seedPt in range(0,seed_pts.GetNumberOfIds()):
                    connected_pts = getConnectedVerticesNotIncludingSeed(centerline, seed_pts.GetId(seedPt))
                    for cpt in range(0,connected_pts.GetNumberOfIds()):
                        seed_pts.InsertUniqueId(connected_pts.GetId(cpt))
                
                
                
                for cpt in range(0,seed_pts.GetNumberOfIds()):
                    diff = [0,0,0]
                    normal = [0,0,0]
                    #norm = checkSurfaceInterior(wall,p,centerline.GetPoint(seed_pts.GetId(cpt)))
                    tangent = centerline.GetPoint(seed_pts.GetId(cpt))
                    normals = wall.GetPointData().GetArray('Normals')
                    normals.GetTuple(p, normal);
                    vtk.vtkMath.Subtract(surface_point, centerline.GetPoint(seed_pts.GetId(cpt)), diff)
                    norm = vtk.vtkMath.Dot(tangent, normal);
                    if norm>=-.01 and norm<0:
                        norm_pts.append(seed_pts.GetId(cpt))
                    if norm>max_norm and norm<0:    
                        max_norm = norm
                        max_norm_pt = seed_pts.GetId(cpt)
                curr_iter += 1
                multiplier += 2
            if len(norm_pts)>0:
                min_norm_pt = norm_pts[0]
                min_distance = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(centerline.GetPoint(norm_pts[0]),surface_point))
                for npt in norm_pts:
                    distance = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(centerline.GetPoint(npt),surface_point))
                    if distance < min_distance:
                        min_distance = distance
                        min_norm_pt = npt
            else:
                min_norm_pt = max_norm_pt
            property_dict = dict()
            for j in range(0,centerline.GetPointData().GetNumberOfArrays()):
                property_dict[str(centerline.GetPointData().GetArray(j).GetName())] = centerline.GetPointData().GetArray(j).GetTuple(min_norm_pt)
            gId = wall.GetPointData().GetArray('GlobalNodeID').GetValue(p)
            properties[gId] = property_dict
        properties = cleanSurfaceData(wall,properties)

    return properties

def getUnitVector(vec):
    unit_vec = [0]*len(vec)
    mag = 0
    for i in vec:
        mag += i**2
    for i in vec:
        unit_vec = i/math.sqrt(mag)
    return unit_vec

def cleanSurfaceData(mesh,properties):
    locator = vtk.vtkCellLocator()
    locator.Initialize();
    locator.SetDataSet(mesh)
    locator.BuildLocator();
    properties_cleaned = dict()
    for i_cell in range(mesh.GetNumberOfCells()):
        i_cellPts = vtk.vtkIdList()
        mesh.GetCellPoints(i_cell,i_cellPts)
        for i_pt in range(0,i_cellPts.GetNumberOfIds()):
            i_gId = mesh.GetPointData().GetArray('GlobalNodeID').GetValue(i_cellPts.GetId(i_pt))
            normal = mesh.GetPointData().GetArray('Normals').GetTuple(i_cellPts.GetId(i_pt))
            i_radius = properties[i_gId]['MaximumInscribedSphereRadius']
            unit_normal = getUnitVector(normal)
            p1 = mesh.GetPoint(i_cellPts.GetId(i_pt))
            p2 = p1 + np.array(unit_normal)*i_radius*2*-1
            cellsAlongLine = vtk.vtkIdList()
            locator.FindCellsAlongLine(p1,p2,0.0000000001,cellsAlongLine)
            #print('Found '+str(cellsAlongLine.GetNumberOfIds())+' cells along the line.')
            for j_cell in range(0,cellsAlongLine.GetNumberOfIds()):
                #print(cellsAlongLine.GetId(j_cell))
                j_cellPts = vtk.vtkIdList()
                mesh.GetCellPoints(cellsAlongLine.GetId(j_cell),j_cellPts)
                for j_pt in range(0,j_cellPts.GetNumberOfIds()): 
                    j_gId = mesh.GetPointData().GetArray('GlobalNodeID').GetValue(j_cellPts.GetId(j_pt))
                    j_radius = properties[j_gId]['MaximumInscribedSphereRadius']
                    if j_radius<0.9*np.array(i_radius):
                        print('Setting '+str(i_radius)+' instead of '+str(j_radius))
                        properties_cleaned[j_gId] = properties[i_gId]
                    else:
                        properties_cleaned[j_gId] = properties[j_gId]
    #exit()
    return properties_cleaned


def createParser():
    parser = argparse.ArgumentParser(description='Maps diameter from given centerline to the surface of a given 3D model.')
    parser.add_argument('centerline', type=str, help='the centerline to map diameters from')
    parser.add_argument('surface', type=str, help='the surface to map onto')
    parser.add_argument('wall_dir', type=str, help='the folder location of the mesh surfaces')
    parser.add_argument('-f','-file', type=str, nargs='?', default = None, help='the pickle filename with data')
    parser.add_argument('-out', type=str, nargs='?', default = 'default.vtu', help='the vtu filename with data')
    parser.add_argument('-o','-option', type=int, nargs='?', const=1, default=0, help='choose different options')
    parser.add_argument('-v', '-verbose', type=int, nargs='?', const=1, default=0, help='turn on verbosity')
    return parser

def main(args):
    mesh = read_polydata(args.surface)
    centerline = read_polydata(args.centerline)
    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInputData(centerline)
    cleaner.PointMergingOn()
    cleaner.Update()
    centerline = cleaner.GetOutput()
    #surface = BranchSurface(surface, centerline, 'MaximumInscribedSphereRadius', 'BranchId', 'BifurcationId')
    # curvature = vtk.vtkCurvatures()
    # curvature.AddInputData(surface)
    # curvature.SetCurvatureTypeToMean()
    # curvature.Update()
    # surface = curvature.GetOutput()
    getCenterlineTangents(centerline)
    #print(centerline.GetPointData().GetArray('Tangent').GetTuple(10))
    if(args.f==None and args.o==0):
        numPts = mesh.GetNumberOfPoints()
        data = [0]*numPts
        centerline_points = set()
        properties = {}
        distances = {}
        centerline_seed_pts = vtk.vtkDoubleArray()
        centerline_seed_pts.SetName('centerline_seed_pts')
        centerline_seed_pts.SetNumberOfValues(mesh.GetNumberOfPoints())
        centerline_seed_pts.Fill(-1)
        for j in range(0,centerline.GetPointData().GetNumberOfArrays()):
                data = vtk.vtkDoubleArray()
                data.SetName(centerline.GetPointData().GetArray(j).GetName())
                data.SetNumberOfValues(mesh.GetNumberOfPoints())
                data.Fill(-1)
                mesh.GetPointData().AddArray(data)
        for i in tqdm(range(0,centerline.GetNumberOfPoints())):
            pt = mesh.FindPoint(centerline.GetPoint(i))
            centerline_points.add(pt)
            centerline_seed_pts.SetValue(pt,i)
            property_dict = dict()
            for j in range(0,centerline.GetPointData().GetNumberOfArrays()):
                property_dict[str(centerline.GetPointData().GetArray(j).GetName())] = centerline.GetPointData().GetArray(j).GetTuple(i)                       
            properties[mesh.FindPoint(centerline.GetPoint(i))] = property_dict
            distances[mesh.FindPoint(centerline.GetPoint(i))] = calcDistance2Points(mesh,mesh.FindPoint(centerline.GetPoint(i)),centerline.GetPoint(i))
        print(properties[mesh.FindPoint(centerline.GetPoint(10))])
        mesh.GetPointData().AddArray(centerline_seed_pts)
        write_polydata(mesh,args.surface.split('.')[0]+'_mapped.vtu')
        graph = generateGraph(mesh)
        mesh = multipleSourceDistance(mesh,graph,-1,centerline_points,distances,properties)
        f = open("file.pkl","wb")
        pickle.dump(graph.node_properties,f)
        f.close()
        mesh = addPropertiesFromGraph(mesh,graph)
    elif args.o==1:
        mesh,properties = mapToNClosestCenterline(centerline,mesh,50)
        f = open("_mapped50.pkl","wb")
        pickle.dump(properties,f)
        f.close()
    else:
        #node_properties = pickle.load(open(args.f,'rb'))
        #mesh = read_polydata(args.surface)
        mesh = generateSurfaceNormals(mesh)
        properties = AssignWallIds(args.wall_dir,mesh,centerline)
        f = open("wallIds.pkl","wb")
        pickle.dump(properties,f)
        f.close()
        mesh = addPropertiesFromDict(mesh,properties)
        #mesh = smoothBifurcations2(centerline,mesh)
        #mesh = smoothAll(mesh)
    
    write_polydata(mesh,args.out)

if __name__ == '__main__':
    parser = createParser()
    args = parser.parse_args()
    main(args)

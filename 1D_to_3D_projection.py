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

def dijsktra_destination(graph, initial, destinations):
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
        if min_node in destinations:
            break
    pbar.close()
    return visited, path, min_node


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
                #multiplier = vtk.vtkMath.Dot(prop_dict['Tangent'],diff)
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

#returns node ids and total distance of the shortest path between two points
def shortest_path(graph, origin, destination):
	visited, paths = dijkstra(graph, origin)
	full_path = deque()
	_destination = paths[destination]

	while _destination != origin:
		full_path.appendleft(_destination)
		_destination = paths[_destination]

	full_path.appendleft(origin)
	full_path.append(destination)

	return visited[destination], list(full_path)

# fast marching method that optimizes calculation of distances using a vtk mesh
def fastMarching(heart_graph,heart,seedPts):
	pt_set= set()
	numPts = heart.GetNumberOfPoints()
	for ptID in seedPts:
		pt_set.add(heart.FindPoint(ptID))

	#intialize edge list
	edgePt = set()
	temp_list = vtk.vtkIdList()
	pt_dist = {}
	for pt in  pt_set:
		connnectedPt_list = getConnectedVerticesNotIncludingSeed(heart,pt)
		for j in range(0,connnectedPt_list.GetNumberOfIds()):
			# new point to decide whether to add to patch, edge, or nothing (if already in edge)
			cpt = connnectedPt_list.GetId(j)
			pt_dist[cpt] = calcDistance2Points(heart,pt,cpt)
			temp_list.InsertNextId(cpt)
	for i in range(0,temp_list.GetNumberOfIds()):
		edgePt.add(temp_list.GetId(i))
		pt_set.add(temp_list.GetId(i))

	temp_list = vtk.vtkIdList()
	#search until all points are found
	while(len(edgePt) > 0):
		temp = set()
		for i in  edgePt:
			connnectedPt_list = getConnectedVerticesNotIncludingSeed(heart,i)
			for j in range(0,connnectedPt_list.GetNumberOfIds()):
				# new point to decide whether to add to patch, edge, or nothing (if already in edge)
				cpt = connnectedPt_list.GetId(j)
				if(cpt in pt_set and cpt in pt_dist):
					pt_set.add(i)
					pt_dist[i] = pt_dist[cpt] + calcDistance2Points(heart,i,cpt)
					heart_graph.add_edge(i,cpt,calcDistance2Points(heart,i,cpt))
				elif(connnectedPt_list.GetId(j) not in pt_set and cpt not in edgePt):
					temp.add(cpt)
		edgePt = temp
	data_array = np.zeros(numPts)
	for i in pt_dist:
		data_array[i] = pt_dist[i]

	vtk_array = vtk.vtkDoubleArray()
	for i in data_array:
		vtk_array.InsertNextValue(i)
	vtk_array.SetName('Point Distances')
	heart.GetPointData().AddArray(vtk_array)

	return pt_dist

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

def createParser():
    parser = argparse.ArgumentParser(description='Maps diameter from given centerline to the surface of a given 3D model.')
    parser.add_argument('centerline', type=str, help='the centerline to map diameters from')
    parser.add_argument('surface', type=str, help='the surface to map onto')
    parser.add_argument('volume', type=str, help='the volume to map onto')
    parser.add_argument('wall_dir', type=str, help='the folder location of the mesh surfaces')
    parser.add_argument('-f','-file', type=str, nargs='?', default = None, help='the pickle filename with data')
    parser.add_argument('-out', type=str, nargs='?', default = 'default.vtu', help='the vtu filename with data')
    parser.add_argument('-o','-option', type=int, nargs='?', default=0, help='choose different options')
    parser.add_argument('-v', '-verbose', type=int, nargs='?', const=1, default=0, help='turn on verbosity')
    return parser

def main(args):
    mesh = read_polydata(args.volume)
    surface = read_polydata(args.surface)
    centerline = read_polydata(args.centerline)
    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInputData(centerline)
    cleaner.PointMergingOn()
    cleaner.Update()
    centerline = cleaner.GetOutput()
    global vprint
    if args.v:
        def vprint(*args):
          # Print each argument separately so caller doesn't need to
          # stuff everything to be printed into a single string
          for arg in args:
            print(arg),
    else:
        vprint = lambda *a: None

    if(args.f==None and args.o==0):
        numPts = mesh.GetNumberOfPoints()
        data = [0]*numPts
        seed_pts = set()
        properties = {}
        distances = {}
        seed_pts_array = vtk.vtkDoubleArray()
        seed_pts_array.SetName('seed_pts_array')
        seed_pts_array.SetNumberOfValues(mesh.GetNumberOfPoints())
        seed_pts_array.Fill(-1)
        for j in range(0,centerline.GetPointData().GetNumberOfArrays()):
            data = vtk.vtkDoubleArray()
            data.SetName(centerline.GetPointData().GetArray(j).GetName())
            data.SetNumberOfValues(mesh.GetNumberOfPoints())
            data.Fill(-1)
            mesh.GetPointData().AddArray(data)
        #Add centerline points as seed points
        print('Adding centerline seed pts')
        for i in tqdm(range(0,centerline.GetNumberOfPoints())):
            pt = mesh.FindPoint(centerline.GetPoint(i))
            seed_pts.add(pt)
            seed_pts_array.SetValue(pt,i)
            property_dict = dict()
            for j in range(0,centerline.GetPointData().GetNumberOfArrays()):
                property_dict[str(centerline.GetPointData().GetArray(j).GetName())] = centerline.GetPointData().GetArray(j).GetTuple(i)                       
            properties[mesh.FindPoint(centerline.GetPoint(i))] = property_dict
            distances[mesh.FindPoint(centerline.GetPoint(i))] = calcDistance2Points(mesh,mesh.FindPoint(centerline.GetPoint(i)),centerline.GetPoint(i))
        #Add surface points as seed points
        # print('Adding surface seed pts')
        # for i in tqdm(range(0,surface.GetNumberOfPoints())):
        #     pt = mesh.FindPoint(surface.GetPoint(i))
        #     if(surface.GetPointData().GetArray('MaximumInscribedSphereRadius').GetValue(pt)>0):
        #         seed_pts.add(pt)
        #         seed_pts_array.SetValue(pt,i)
        #         property_dict = dict()
        #         for j in range(0,surface.GetPointData().GetNumberOfArrays()):
        #             property_dict[str(surface.GetPointData().GetArray(j).GetName())] = surface.GetPointData().GetArray(j).GetTuple(i)                       
        #         properties[mesh.FindPoint(surface.GetPoint(i))] = property_dict
        #         distances[mesh.FindPoint(surface.GetPoint(i))] = calcDistance2Points(mesh,mesh.FindPoint(surface.GetPoint(i)),surface.GetPoint(i))

        print('Found '+str(len(list(properties)))+' seed pts.')
        
        mesh.GetPointData().AddArray(seed_pts_array)
        write_polydata(mesh,args.surface.split('.')[0]+'_mapped.vtu')
        graph = generateGraph(mesh)
        mesh = multipleSourceDistance(mesh,graph,-1,seed_pts,distances,properties)
        f = open("file.pkl","wb")
        pickle.dump(graph.node_properties,f)
        f.close()
        #mesh = addPropertiesFromDict(mesh,graph.node_properties)
    elif(args.o==0):
        f = open("file.pkl","rb")
        props = pickle.load(f)
        mesh = addPropertiesFromDict(mesh,props)
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
from vtk.util import numpy_support as VN

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

def convert_np_array_to_vtk(name,np_array):
	data_array = vtk.vtkDoubleArray()
	data_array.SetName(name)
	for i in range(0,len(np_array)):
		data_array.InsertNextValue(np_array[i])
	return data_array

def createParser():
	parser = argparse.ArgumentParser(description='Segment volume of tissue perfused.')
	parser.add_argument('heart', type=str, help='the heart filename (include file ext)')
	parser.add_argument('vtu_data', type=str, help='the output vtu filename (include file ext)')
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
	heart = read_polydata(args.heart)
	numArrays = heart.GetPointData().GetNumberOfArrays()
	numPts = heart.GetNumberOfPoints()
	cap_ids = dict()
	cap_names = set()
	for i in range(0,numArrays):
		array_name = heart.GetPointData().GetArray(i).GetName().split('_')
		if len(array_name)>0 and array_name[0].isdigit():
			cap_name = '_'.join(array_name[1:])
			cap_names.add(cap_name)
			cap_ids[cap_name] = int(array_name[0])
	primary_branch_length = 2
	primary_caps = set()
	for i in cap_names:
		if len(i.split('_'))==primary_branch_length:
			primary_caps.add(i)
	group_arrays = dict()
	for i in primary_caps:
		group_arrays[i] = np.zeros(numPts)

	for i in range(0,numArrays):
		array_name = heart.GetPointData().GetArray(i).GetName().split('_')[1:]
		for j in primary_caps:
			if '_'.join(array_name).startswith(j):
				temp_array = heart.GetPointData().GetArray(i)
				for k in  range(0,numPts):
					if temp_array.GetValue(k) == 1:
						group_arrays[j][k] = 1

	sorted_cap_names = sorted(primary_caps)
	for i in range(0,len(sorted_cap_names)):
		vtk_name = 'primary_'+str(i)+'_'+sorted_cap_names[i]
		heart.GetPointData().AddArray(convert_np_array_to_vtk(vtk_name,group_arrays[sorted_cap_names[i]]))

	write_polydata(heart,args.vtu_data)
	exit()
	volumes_dict = {}
	for i in range(0,len(cap_names)):
		volumes_dict[cap_names[i]] = volumes[i]

	return volumes_dict

if __name__ == '__main__':
	parser = createParser()
	args = parser.parse_args()
	main(args)
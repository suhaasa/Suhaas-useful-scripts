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
from os import path

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
    if not path.exists(filename):
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


def intializeVTP(filename):
  datareader=vtk.vtkXMLPolyDataReader()
  datareader.SetFileName(filename)
  datareader.Update()

  mesh=vtk.vtkDataSetMapper()
  mesh=datareader.GetOutput()
  vprint('Loaded .vtp file.')
  return mesh

def intializeVTU(filename):
	datareader=vtk.vtkXMLUnstructuredGridReader()
	datareader.SetFileName(filename)
	datareader.Update()

	mesh=datareader.GetOutput()
	vprint('Loaded .vtu file.')
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

def printPoints(mesh,filename):
	outfile = open(filename,'a+')
	for i in range(0,mesh.GetNumberOfPoints()):
		x,y,z = mesh.GetPoint(i)
		out_string = str(x) + ' ' + str(y) + ' ' + str(z) + '\n'
		outfile.write(out_string)
	outfile.close()

def transform(mesh,matrix):
	t = vtk.vtkTransform()
	t.SetMatrix(matrix)
	# Apply the transformation to the grid
	tf = vtk.vtkTransformFilter()
	tf.SetInputData(mesh)
	tf.SetTransform(t)
	tf.Update()
	return tf.GetOutput()

def transformScale(mesh,scale):
	t = vtk.vtkTransform()
	t.Scale(scale)
	# Apply the transformation to the grid
	tf = vtk.vtkTransformFilter()
	tf.SetInputData(mesh)
	tf.SetTransform(t)
	tf.Update()
	return tf.GetOutput()

def scaleNormalizeMesh(mesh1,mesh2):
	x1min,x1max,y1min,y1max,z1min,z1max = mesh1.GetBounds()
	x2min,x2max,y2min,y2max,z2min,z2max = mesh2.GetBounds()
	x1 = abs(x1max-x1min)
	y1 = abs(y1max-y1min)
	z1 = abs(z1max-z1min)
	x2 = abs(x2max-x2min)
	y2 = abs(y2max-y2min)
	z2 = abs(z2max-z2min)
	scale = x1/x2,y1/y2,z1/z2
	vprint('Scaling mesh2 by: ' + str(scale))
	return transformScale(mesh2,scale)

def scaleLengthNormalizeMesh(mesh1,mesh2):
	l1 = mesh1.GetLength()
	l2 = mesh2.GetLength()
	scale = l1/l2,l1/l2,l1/l2
	scale = .1,.1,.1
	vprint('Scaling mesh2 by: ' + str(scale))
	return transformScale(mesh2,scale)

def generateTransformationMatrix(filename):
	#np_array = np.fromfile(filename)
	np_array = [[0.48948553305373493, -0.7091338901142712, 0.5074771313295836, 1.1491139275391937],[0.5487440753422649, 0.702777410906223, 0.45275142241268557, -10.780479512718845], [-0.6777048418801989, 0.05685979784828735, 0.7331323964201746, 11.215700750832157],[ 0.0, 0.0, 0.0, 1.0]]
	np_array = [[0.48948553305373493, -0.7091338901142712, 0.5074771313295836, -6.439346051233096], [0.5487440753422649, 0.702777410906223, 0.45275142241268557, -5.657732655162194], [-0.6777048418801989, 0.05685979784828735, 0.7331323964201746, 4.430907416055166], [0.0, 0.0, 0.0, 1.0]]
	np_array = [[0.7222571409790202, -0.410996740091788, 0.556261001633887, 1.1491139275391937], [0.6014632037643727, 0.7703142283189395, -0.21179708252687096, -10.780479512718845], [-0.3414478537380325, 0.4875424794654519, 0.8035643682334171, 11.215700750832157], [0.0, 0.0, 0.0, 1.0]]
	np_array = [[0.8167329246237517, -0.5235112268836243, 0.24265886582264673, -115.72449056133871, 0.40153475872208555, 0.8176626827072208, 0.41255008768151696, 138.3790934539467, -0.4143877017643343, -0.2395072705260391, 0.8780199883781922, 77.7852930881601, 0.0, 0.0, 0.0, 1.0]]
	np_array = [[0.9998957181259133, -0.0007343638144923404, 0.014422676008678557, -0.32345160531015177], [0.0008424376786968768, 0.9999716045507724, -0.007488684190785576, 0.3927814792611876],[-0.014416767051626374, 0.007500053462461068, 0.9998679442935645, 2.403121442209475], [0.0, 0.0, 0.0, 1.0]]
	matrix = vtk.vtkMatrix4x4()
	for i in range(0,len(np_array)):
		for j in range(0,len(np_array[i])):
			matrix.SetElement(i,j,np_array[i][j])
	return matrix

def createParser():
	parser = argparse.ArgumentParser(description='Transforms the file1 into file2 coords.')
	parser.add_argument('file1', type=str, help='the first file')
	parser.add_argument('file2', type=str, help='the second file')
	parser.add_argument('outfile1', type=str, help='the output file with coordinates')
	parser.add_argument('outfile2', type=str, help='the output file with coordinates')
	parser.add_argument('-t', '-transform', type=int, nargs='?', const=1, default=0, help='turn on transformation')
	parser.add_argument('-v', '-verbose', type=int, nargs='?', const=1, default=0, help='turn on verbosity')
	return parser

def main(args):
	file1 = args.file1
	file2 = args.file2
	outfile1 = args.outfile1
	outfile2 = args.outfile2
	global vprint
	if args.v:
	    def vprint(*args):
	      # Print each argument separately so caller doesn't need to
	      # stuff everything to be printed into a single string
	      for arg in args:
	        print(arg),
	else:
		vprint = lambda *a: None
	mesh1 = read_polydata(file1)
	mesh2 = read_polydata(file2)
	#ensure the two meshes have relatively the same size, so normalize the scale of mesh2 to mesh1
	mesh2 = scaleLengthNormalizeMesh(mesh1,mesh2)
	if(args.t):
		matrix = generateTransformationMatrix('mat')
		transformed_mesh2 = transform(mesh2,matrix)
		write_polydata(transformed_mesh2,'transformed_mesh2.vtu')
	printPoints(mesh1,outfile1)
	printPoints(mesh2,outfile2)
	write_polydata(mesh1,'mesh1.vtp')
	write_polydata(mesh2,'mesh2.vtp')


if __name__ == '__main__':
	parser = createParser()
	args = parser.parse_args()
	main(args)
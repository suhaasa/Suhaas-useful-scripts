import sys
import os
import vtk
import numpy as np
import time
import glob
import math
import argparse

def calcDistance2Points(model, pt1,pt2):
	x1,y1,z1 = model.GetPoint(pt1)
	x2,y2,z2 = pt2[0],pt2[1],pt2[2]
	distance = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**(.5)
	return distance

def calculateCapCenters(caps):
	cap_centers = []
	for cap in caps:
		numPts = cap.GetNumberOfPoints()
		center_x,center_y,center_z = 0
		for p in range(0,numPts):
			x,y,z = mesh.GetPoint(meshPointID)
			center_x += x
			center_y += y
			center_z += z
		center_x = center_x/numPts
		center_y = center_y/numPts
		center_z = center_z/numPts
		cap_centers.append([center_x,center_y,center_z])

def minDistanceBetweenPoints(model, seedPt, pt_list):
	min = 100000
	min_pt_index = 0
	for iPt in range(0,len(ptlist)):
		distance = calcDistance2Points(model, seedPt,pt_list[iPt])
		if(distance < min):
			min = distance
			min_pt_index = iPt
	return min,min_pt_index

def determinePerfusionVolumes(caps,cap_center_coordinates,image,threshold,vprint):
	numPts = image.GetNumberOfPoints()
	image_data = [0]*numPts
	for ip in range(0,numPts):
		update_progress(ip,numPts,vprint)
		if(image.GetPointData().GetArray(0).GetValue(ip)>threshold):
			min,min_pt_index = minDistanceBetweenPoints(model, seedPt, pt_list):
			image_data[ip] = min_pt_index

	data_array = vtk.vtkDoubleArray()
	data_array.SetName('PerfusionVolumes')
	for ptID in range(0,numPts):
		data_array.InsertNextValue(image_data[ptID])
	image.GetPointData().AddArray(data_array)

def update_progress(progress, total, vprint):  
	vprint('\r[{0:10}]{1:>2}'.format('#' * int(progress * 10 /total), progress))

def createParser():
	parser = argparse.ArgumentParser(description='Finds volume of tissue perfused by each outlet.')
	parser.add_argument('caps', type=str, help='the input model cap locations')
	parser.add_argument('image_data', type=str, help='the image filename (include file ext)')
	parser.add_argument('image_data', type=str, help='the image filename (include file ext)')
	parser.add_argument('-v', '-verbose', type=int, nargs='?', const=1, default=0, help='turn on verbosity')
	return parser

def main(args):
	IMAGE_FILENAME = args.image_data
	#Read vti file
	ref = vtk.vtkXMLPolyDataReader()
	ref.SetFileName(IMAGE_FILENAME)
	ref.Update()
	#Read your data into another polydata variable for reading
	image=vtk.vtkPolyData()
	image=ref.GetOutput()
	# vtp files of the caps to calculate
	caps = []
	for file in os.listdir(args.caps):
		if file.endswith('.vtp') and file.beginswith('cap_'):
			#Read vtp file
			ref = vtk.vtkXMLPolyDataReader()
			ref.SetFileName(file)
			ref.Update()
			#Read your data into another polydata variable for reading
			cap=vtk.vtkPolyData()
			cap=ref.GetOutput()
			caps.append(cap)

	vprint('Found ' + str(len(caps)) + ' caps.')

	cap_center_coordinates = calculateCapCenters(caps)

	determinePerfusionVolumes(caps,cap_center_coordinates,image)

if __name__ == '__main__':
	parser = createParser()
	args = parser.parse_args()
	main(args)
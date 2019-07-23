import sys
import os
import vtk
import numpy as np
import time
import glob
import math
import argparse

def writeFlowFile(filename,flow):
	outfile = open(filename,'w')
	out_string = '# Time (sec)\tFlow (micrometers^3/sec)' + '\n'
	outfile.write(out_string)
	for i in range(0,len(column1)):
		out_string = '0.000000000' + '\t' + str(flow) + '\n'
		outfile.write(out_string)
		out_string = '1.000000000' + '\t' + str(flow) + '\n'
		outfile.write(out_string)
	outfile.close()

def createParser():
	parser = argparse.ArgumentParser(description='Finds volume of tissue perfused by each outlet.')
	parser.add_argument('caps', type=str, help='the input model cap locations')
	parser.add_argument('data', type=str, help='the output filename (include file ext)')
	parser.add_argument('-v', '-verbose', type=int, nargs='?', const=1, default=0, help='turn on verbosity')
	parser.add_argument('-f', '-flip', type=int, nargs='?', const=1, default=0, help='turn on bct_flip')
	return parser

def main(args):
	DATA_FILENAME = args.data
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
	caps = []
	cap_names = []
	for file in os.listdir(args.caps):
		if file.endswith('.vtp'):
			cap_names.append(file[:len(file)-4])

	vprint('Found ' + str(len(caps)) + ' caps.')
	outfile = open(DATA_FILENAME,'w')
	for i in range(0,len(cap_names)):
		out_string = 'prescribed_velocities_vtp mesh-complete/mesh-surfaces/' + cap_names[i] + '.vtp' + '\n'
		out_string = out_string + 'bct_analytical_shape plug' + '\n'
		out_string = out_string + 'bct_period 1.000000000' + '\n'
		out_string = out_string + 'bct_point_number 201' + '\n'
		out_string = out_string + 'bct_fourier_mode_number 10' + '\n'
		if(args.f):
			out_string = out_string + 'bct_flip' + '\n'
		out_string = out_string + 'bct_create mesh-complete/mesh-surfaces/' + cap_names[i] + '.vtp ' + './flow-files/' + cap_names[i] + '.flow' + '\n'
		outfile.write(out_string)
	outfile.close()


if __name__ == '__main__':
	parser = createParser()
	args = parser.parse_args()
	main(args)

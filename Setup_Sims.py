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

def writePresolver(sim,KW):
	old = open(vos.cwd + 'Sim_v10.svpre', 'r')
	new = open(os.cwd + sim + '/Sim_v10.svpre','w')
	
	line = old.readline()
	while line:
		if(line=='set_surface_id_vtp mesh-complete/mesh-surfaces/LCAc_75.vtp 73'):
			line = 'set_surface_id_vtp mesh-complete/mesh-surfaces/LCAc' + KW + '.vtp 73'
		new_solv.write(line)
		line = old_solv.readline()
	old_solv.close()
	new_solv.close()

def createParser():
	parser = argparse.ArgumentParser(description='Maps flow from each outlet to heart tissue.')
	parser.add_argument('heart', type=str, help='the input model heart')
	return parser

def main(args):
	master_file = 'Simulations.txt'
	sim_dict
	for sim in sim_dict:
		writePresolver(sim, sim_dict[sim])
	writeSolver()

if __name__ == '__main__':
	parser = createParser()
	args = parser.parse_args()
	main(args)
# PDBDelaunay.py
# Author: Yogesh Joshi.

import os
import csv
import numpy as np
from Bio import PDB
import scipy.spatial as spatial
import scipy.spatial.distance as distance
from scipy.spatial import ConvexHull, Delaunay, Voronoi

pdbid = '1crn'
chainid = 'A'
filename = 'pdb'+pdbid+'.ent'
print (filename)


def List_to_CSV(OutFname, DataList):
	""" Writes the contents of a list to a CSV File"""
	with open(OutFname, 'w') as myfile:
		wr = csv.writer(myfile, delimiter=',')
		wr.writerows(line for line in DataList)


def PointCloudData(pdbid, chainid):
	""" Get C-alpha coordinates for the given pdbid and chainid
	along with the temperature factors and residue names.
	"""
	pc = []
	bf = []
	resnames = []

	if not os.path.exists(os.getcwd() + '/' + filename):
		pdbl = PDB.PDBList()
		pdbl.retrieve_pdb_file(pdbid, False, os.getcwd(), 'pdb', True)
	parser = PDB.PDBParser(PERMISSIVE=1)
	structure = parser.get_structure(pdbid, 'pdb'+pdbid+'.ent')
	model = structure[0]
	chain = model[chainid]
	for residue in chain:
		for atom in residue:
			if atom.get_id() == "CA":
				resnames.append(residue.get_resname())
				bf.append(atom.get_bfactor())
				pc.append(atom.get_coord())
	pointcloud = np.asarray(pc)
	return pointcloud, bf, resnames


def ProteinDelaunay(pointcloud, pdbid, chain):
	""" Generate the Delaunay Tessellation of all the points in the given point cloud.
	The point cloud is basically the Calpha coordinates of a three dimensional protein.
	The point, vertices, simplices and neighbors in the entire point cloud are obtained
	as arrays.
	"""
	# Convex Hull.
	ConHull = ConvexHull(pointcloud)
	hullArea = round(ConHull.area, 4)
	hullVolume = round(ConHull.volume, 4)




pc, bf, rsname = PointCloudData(pdbid, chainid)
ProteinDelaunay(pc, pdbid, chainid)
# print (pc, len(pc))
# print (bf, len(bf))
# print (rsname, len(rsname))




























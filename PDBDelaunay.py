# PDBDelaunay.py
# Author: Yogesh Joshi.

import os
import csv
import math
import numpy as np
from Bio import PDB
import amino_dict
import scipy.spatial as spatial
import scipy.spatial.distance as distance
from scipy.spatial import ConvexHull, Delaunay, Voronoi

pdbid = '1crn'
chainid = 'A'
filename = 'pdb' + pdbid + '.ent'


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


def DevTetra(lengths):
    AvgL = (sum(lengths) / len(lengths))
    T = 0.0
    i = 0
    while i < 5:
        j = i + 1
        while j < 6:
            T += (lengths[i] - lengths[j])**2
            j = j + 1
        i = i + 1
    den = (15.0 * (AvgL ** 2))
    return round (float(T / den), 4)


def ProteinDelaunay(pdbid, chain):
    """ Generate the Delaunay Tessellation of all the points in the given point cloud.
    The point cloud is basically the Calpha coordinates of a three dimensional protein.
    The point, vertices, simplices and neighbors in the entire point cloud are obtained
    as arrays.
    """
    Data = []
    Head = ['PDBID', 'Quad', 'SortedQuad', 'RedAlpha', 'SortRedAlpha', 'V1', 'V2', 'V3', 'V4', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'SumL', 'AvgL', 'DevL', 'DevTetra', 'Vol', 'TF1', 'TF2', 'TF3', 'TF4', 'SumTF', 'AvgTF']
    Data.append(Head)
    
    pointcloud, bf, resname = PointCloudData(pdbid, chainid)

    # Convex Hull.
    ConvxHull = ConvexHull(pointcloud)
    hullArea = round(ConvxHull.area, 4)
    hullVolume = round(ConvxHull.volume, 4)

    # Delaunay Tessellation
    delaunay_hull = Delaunay(pointcloud, furthest_site=False, incremental=False, qhull_options='Qc') # noqa E501
    delaunay_points = delaunay_hull.points
    delaunay_vertices = delaunay_hull.vertices
    delaunay_simplices = delaunay_hull.simplices
    delaunay_neighbors = delaunay_hull.neighbors

    for i in delaunay_vertices:

        # Obtain the indices of the vertices.
        one, two, three, four = i[2], i[1], i[3], i[0]

        # Obtain the coordinates based on the indices.
        cordA = pointcloud[one]
        cordB = pointcloud[two]
        cordC = pointcloud[three]
        cordD = pointcloud[four]

        # Get three letter amino acid names based on indices.
        a = resname[one]
        b = resname[two]
        c = resname[three]
        d = resname[four]

        # Get the temprature factors for the amino acids.
        a_tf = bf[one]
        b_tf = bf[two]
        c_tf = bf[three]
        d_tf = bf[four]

        # Get the string of three letter amino acids
        # forming the vertices of the tetrahedra.
        amino = [a, b, c, d]
        sortAmino = sorted(amino)
        amino = '-'.join(amino)
        sortAmino = '-'.join(sortAmino)

        # Get one letter code of the amino acids
        oneA = amino_dict.replace_all(a, amino_dict.one_letter)
        oneB = amino_dict.replace_all(b, amino_dict.one_letter)
        oneC = amino_dict.replace_all(c, amino_dict.one_letter)
        oneD = amino_dict.replace_all(d, amino_dict.one_letter)
        oneLet = [oneA, oneB, oneC, oneD]
        sortOneLet = sorted(oneLet)
        oneLet = ''.join(oneLet)
        sortOneLet = ''.join(sortOneLet)

        # Get Reduced Amino Acid Representations.
        flpA = amino_dict.replace_all(oneA, amino_dict.FLP)
        flpB = amino_dict.replace_all(oneB, amino_dict.FLP)
        flpC = amino_dict.replace_all(oneC, amino_dict.FLP)
        flpD = amino_dict.replace_all(oneD, amino_dict.FLP)
        flp = [flpA, flpB, flpC, flpD]
        sortflp = sorted(flp)
        flp = (''.join(flp)).upper()
        sortflp = (''.join(sortflp)).upper()

        # Calculate distances between the tetrahedra vertices.
        AB = np.linalg.norm(cordA - cordB)
        AC = np.linalg.norm(cordA - cordC)
        AD = np.linalg.norm(cordA - cordD)
        BC = np.linalg.norm(cordB - cordC)
        BD = np.linalg.norm(cordB - cordD)
        CD = np.linalg.norm(cordC - cordD)

        # Calculate the tetrahedra Volume.
        A_prime = cordA - cordD
        B_prime = cordB - cordD
        C_prime = cordC - cordD
        primes = [A_prime, B_prime, C_prime]
        primes = np.asarray(primes)
        det = np.linalg.det(primes)
        Vol = round((abs(det) / 6), 4)

        # Sum of Edge Lengths.
        SumL = (AB + AC + AD + BC + BD + CD)
        SumL = round(SumL, 4)

        # Average Edge Lengths.
        AvgL = round((SumL / 6), 4)

        # Deviation in Edge Lengths.
        devLp = (AB - AvgL) ** 2
        devLq = (AC - AvgL) ** 2
        devLr = (AD - AvgL) ** 2
        devLs = (BC - AvgL) ** 2
        devLt = (BD - AvgL) ** 2
        devLu = (CD - AvgL) ** 2
        devLy = [devLp, devLq, devLr, devLs, devLt, devLu]
        sumDevL = sum(devLy)
        DevL = round(math.sqrt(sumDevL / 6.0), 4)

        # Deviation in Tetrahedrality
        lenArr = [AB, AC, AD, BC, BD, CD]
        DevT = DevTetra(lenArr)

        # Sum and Average Temperature Factors.
        SumTF = (a_tf + b_tf + c_tf + d_tf)
        AvgTF = round(SumTF / 4, 4)

        # Data List
        line =  [pdbid, oneLet, sortOneLet, flp, sortflp, one, two, three, four, AB, AC, AD, BC, BD, CD, SumL, AvgL, DevL, DevT, Vol, a_tf, b_tf, c_tf, d_tf, SumTF, AvgTF]
        Data.append(line)
    return Data


ProDelaunay = ProteinDelaunay(pdbid, chainid)
List_to_CSV(pdbid + '_Delaunay.csv', ProDelaunay)


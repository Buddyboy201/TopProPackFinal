import os
import sys
import numpy as np
import scipy.spatial
import networkx as nx
import csv

#graphid -> coords
#coords -> resid
#resid -> res

#TODO: extra step: graphid -> coords -> resid | set(resid) -> res

def generate_hashmaps(file):
    atom_resid = {}
    resid_res = {}
    with open(file) as pdb_file:
        for line in pdb_file:
            if line[0:4] == "ATOM":
                num = int(line[6:11].strip(" "))
                residue = line[17:20].strip(" ")
                res_id = int(line[22:26].strip(" "))
                atom_resid[num] = res_id
                resid_res[res_id] = residue
    return atom_resid, resid_res
                
def generate_hashmaps_and_cliques(file):
    #atom_resid = {}
    coords_resid = {}
    resid_res = {}
    coords = []
    with open(file) as pdb_file:
        for line in pdb_file:
            if line[0:4] == "ATOM" and line[12:16].strip(" ") != "CA" and line[12:16].strip(" ") != "O" and line[12:16].strip(" ") != "N" and line[12:16].strip(" ") != "C":
                num = int(line[6:11].strip(" "))
                residue = line[17:20].strip(" ")
                coordx = float(line[30:38].strip(" "))
                coordy = float(line[38:46].strip(" "))
                coordz = float(line[46:54].strip(" "))
                res_id = int(line[22:26].strip(" "))
                #atom_resid[num] = res_id
                coords_resid[(coordx, coordy, coordz)] = res_id
                resid_res[res_id] = residue
                coords.append((coordx, coordy, coordz))
    coords_array = np.array(coords)
    del coords
    tri = scipy.spatial.qhull.Delaunay(coords_array)
    edges = list()
    for n in tri.simplices:
        edge = sorted([n[0], n[1]])
        edges.append((edge[0], edge[1]))
        edge = sorted([n[0], n[2]])
        edges.append((edge[0], edge[1]))
        edge = sorted([n[0], n[3]])
        edges.append((edge[0], edge[1]))
        edge = sorted([n[1], n[2]])
        edges.append((edge[0], edge[1]))
        edge = sorted([n[1], n[3]])
        edges.append((edge[0], edge[1]))
        edge = sorted([n[2], n[3]])
        edges.append((edge[0], edge[1]))
    graph = nx.Graph(edges)
    cliques = np.array(list(nx.find_cliques(graph)))
    return coords_array, coords_resid, resid_res, cliques

def update_cliques(cliques_arr, coords_array, coords_resid, resid_res):
    temp_arr = []
    for i in range(len(cliques_arr)):
        for j in range(len(cliques_arr[i])):
            cliques_arr[i][j] = coords_resid[tuple(coords_array[cliques_arr[i][j]])]
        cliques_arr[i] = list(set(cliques_arr[i]))
        cliques_arr[i].sort()
        if cliques_arr[i] not in temp_arr:
            temp_arr.append(cliques_arr[i])
    cliques_arr = np.array(temp_arr)
    for i in range(len(cliques_arr)):
        for j in range(len(cliques_arr[i])):
            #cliques_arr[i][j] = (cliques_arr[i][j], resid_res[cliques_arr[i][j]])
            cliques_arr[i][j] = resid_res[cliques_arr[i][j]]
    cliques_arr = list(cliques_arr)
    cliques_arr.sort(key=len)
    cliques_arr = np.array(cliques_arr)
    return cliques_arr

def update_cliques_2(cliques_arr, coords_array, coords_resid):
    temp_arr = []
    for i in range(len(cliques_arr)):
        for j in range(len(cliques_arr[i])):
            cliques_arr[i][j] = coords_resid[tuple(coords_array[cliques_arr[i][j]])]
        cliques_arr[i] = list(set(cliques_arr[i]))
        cliques_arr[i].sort()
        if cliques_arr[i] not in temp_arr:
            temp_arr.append(cliques_arr[i])
    temp_arr.sort(key=len)
    cliques_arr = np.array(temp_arr)
    return cliques_arr
            

temp = '''def main():
    os.chdir("c:\\pdb_data\\domains")
    files = os.listdir()
    #print(files)
    for i in range(0, len(files)-12, 12):
        pdb_file = files[i+8]
        cliques_file = files[i+4]
        print(pdb_file, cliques_file)
        coords_array, coords_resid, resid_res, cliques = generate_hashmaps_and_cliques(pdb_file)
        cliques = update_cliques(cliques, coords_array, coords_resid, resid_res)
        for i in cliques:
            print(i)'''

#main()

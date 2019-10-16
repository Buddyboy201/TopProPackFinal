import os
import sys
import numpy as np
import scipy.spatial
import networkx as nx
import csv
import COM_calc as COMgen

test_file = "c:\\alpha\\4quv.pdb"

def generate_hashmaps_and_cliques(file):
    centroids, coords_resid, resid_centroids = COMgen.get_COM_all(file)
    resid_res = {}
    with open(file) as pdb_file:
        for line in pdb_file:
            if line[0:4] == "ATOM": #and line[12:16].strip(" ") != "CA" and line[12:16].strip(" ") != "O" and line[12:16].strip(" ") != "N" and line[12:16].strip(" ") != "C":
                residue = line[17:20].strip(" ")
                res_id = int(line[22:26].strip(" "))
                resid_res[res_id] = residue
    tri = scipy.spatial.qhull.Delaunay(centroids)
    edges = list()
    for n in tri.simplices:
        edge = sorted([n[0], n[1]])
        print((edge[0], edge[1]))
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
    return cliques, coords_resid, resid_res, resid_centroids

def update_cliques(cliques, coords_resid, resid_res):
    for i in range(len(cliques)):
        for j in range(len(cliques[i])):
            cliques[i][j] = (coords_resid[cliques[i][j]], resid_res[coords_resid[cliques[i][j]]])
    return cliques

cliques, coords_resid, resid_res, resid_centroids = generate_hashmaps_and_cliques(test_file)
#cliques = update_cliques(cliques, coords_resid, resid_res)

def generate_centroid_cliques(file):
    cliques, coords_resid, resid_res, resid_centroids = generate_hashmaps_and_cliques(file)
    cliques = update_cliques(cliques, coords_resid, resid_res)
    return cliques, resid_centroids
#print(cliques)

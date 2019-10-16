import os
import sys
import numpy as np
import scipy.spatial
import networkx as nx


test_file = "c:\\alpha\\4quv.pdb" #confirmed alpha-helical membrane protein


#res name
#res id
#res id -> res
#atoms -> res id
# 1) convert all atoms to res ids
# 2) create a set out of res ids
# 3) convert all res ids in set to res

def generate_hashmaps_and_cliques(file):
    atom_resid = {}
    resid_res = {}
    coords = []
    with open(file) as pdb_file:
        for line in pdb_file:
            if line[0:4] == "ATOM":
                num = int(line[6:11].strip(" "))
                residue = line[17:20].strip(" ")
                coordx = float(line[30:38].strip(" "))
                coordy = float(line[38:46].strip(" "))
                coordz = float(line[46:54].strip(" "))
                res_id = int(line[22:26].strip(" "))
                atom_resid[num] = res_id
                resid_res[res_id] = residue
                coords.append([coordx, coordy, coordz])
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
    cliques = list(nx.enumerate_all_cliques(graph))
    print(len(cliques))
    return  coords_array, atom_resid, resid_res, cliques

def update_cliques(cliques_arr, atom_resid, resid_res):
    for i in range(len(cliques_arr)):
        for j in range(len(cliques_arr[i])):
            cliques_arr[i][j] = atom_resid[cliques_arr[i][j]+1]
        cliques_arr[i] = list(set(cliques_arr[i]))
        for j in range(len(cliques_arr[i])):
            cliques_arr[i][j] = resid_res[cliques_arr[i][j]]
    return cliques_arr

def main():
    os.chdir("c:\\alpha")
    #if len(sys.argv) > 2:
    #    raise Exception("Two arguments were expected; More than 2 arguments were given")
    #elif len(sys.argv) == 1:
    #    raise Exception("Two arguments were expected; Less than 2 arguments were given")
    if False: continue
    else:
        coords_array, atom_resid, resid_res, cliques = generate_hashmaps_and_cliques(test_file)
        cliques = update_cliques(cliques, atom_resid, resid_res)
    for i in cliques:
        if len(i) > 4:
            print(i)
            
main()


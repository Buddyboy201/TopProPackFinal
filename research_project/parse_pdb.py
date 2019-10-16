import os
import sys
import numpy as np
import scipy.spatial
import networkx as nx

graph = None

def update_hashmap(hashmap, num, residue, coords):
    hashmap[num] = {"res":residue, "coords":coords}
    

def generate_hashmap_and_cliques(file):
    values = {}
    coords = []
    allcliques = None
    with open(file) as pdb_file:
        for line in pdb_file:
            if line[0:4] == "ATOM":
                num = line[6:11].strip(" ")
                residue = line[17:20].strip(" ")
                coordx = float(line[30:38].strip(" "))
                coordy = float(line[38:46].strip(" "))
                coordz = float(line[46:54].strip(" "))
                update_hashmap(values, num, residue, (coordx, coordy, coordz))
                coords.append([coordx, coordy, coordz])
    coords_array = np.array(coords)
    tri = scipy.spatial.qhull.Delaunay(coords_array)
    edges = list()
    print(
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
    #print(graph.nodes)
    #cliques = nx.find_cliques(graph)
    #numb = nx.graph_number_of_cliques(graph)
    allcliques = list(nx.enumerate_all_cliques(graph))
    return values, allcliques


def update_cliques(cliques, values):
    new_cliques = []
    for clique in cliques:
        new_clique = []
        for atom in clique:
            new_clique.append(values[str(atom+1)]["res"])
        new_cliques.append(new_clique)
    return new_cliques
    


def main():
    #os.chdir("c:\\alpha")
    #if len(sys.argv) > 2:
    #    raise Exception("Two arguments were expected; More than 2 arguments were given")
    #elif len(sys.argv) == 1:
       # raise Exception("Two arguments were expected; Less than 2 arguments were given")
    if False: continue
    else:
        values, cliques = generate_hashmap_and_cliques("c:\\alpha\\4quv.pdb")
        #cliques = update_cliques(cliques, values)
        #print(values)
        #print("\n\n\n\n")
        #print(cliques[5000:5500])

main()        

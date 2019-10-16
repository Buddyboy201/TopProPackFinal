import os
import sys
import math
import numpy as np
import centroid_based_clique_generator as centroidGen

test_file = "c:\\alpha\\4quv.pdb"

def get_dist(coord1, coord2):
    return math.sqrt((coord1[0]-coord2[0])**2+(coord1[1]-coord2[1])**2+(coord1[2]-coord2[2])**2)

def getMaxMinDistances(coords):
    max_dist = 0
    min_dist = 10000
    for i in range(len(coords)-1):
        for j in range(i+1, len(coords)):
            dist = get_dist(coords[i], coords[j])
            max_dist = max([max_dist, dist])
            min_dist = min([min_dist, dist])
    return max_dist, min_dist

def convert_resid_coords(clique, resid_centroids):
    for i in range(len(clique)):
        clique[i] = resid_centroids[clique[i]]
    return clique

def getProteinMinMax(file):
    cliques, resid_centroids = centroidGen.generate_centroid_cliques(file)
    new_cliques = []
    for clique in cliques:
        temp = []
        for i in clique:
            temp.append(resid_centroids[i[0]])
        new_cliques.append(temp)
    cliques = np.array(new_cliques)
    for i in cliques:
        max_dist, min_dist = getMaxMinDistances(i)
        print("min: {} max: {}".format(min_dist, max_dist))

getProteinMinMax(test_file)



import atom
import residue
import os
import sys
import numpy as np
import scipy.spatial
import networkx as nx
import csv
import math
import time
import matplotlib.pyplot as plt
import statistics
import streamlit
import pandas as pd




def binary_search(arr, val):
    start = 0
    end = len(arr) - 1
    while True:
        mid = int((start + end) / 2)
        # print(start, end)
        if start >= end and arr[mid] != val: return -1
        if arr[mid] == val:
            return mid
        elif val > arr[mid]:
            start = mid + 1
        else:
            end = mid - 1


def get_one_var_stats(data):
    data.sort()
    x_bar = statistics.mean(data)
    med = statistics.median(data)
    std_dev = statistics.stdev(data, xbar=x_bar)
    mode = statistics.mode(data)
    data_range = data[len(data) - 1] - data[0]
    return x_bar, med, mode, data_range, std_dev


def get_dist(coord1, coord2):
    return math.sqrt((coord1[0] - coord2[0]) ** 2 + (coord1[1] - coord2[1]) ** 2 + (coord1[2] - coord2[2]) ** 2)


class CentroidProtein:
    def __init__(self, name, file_path, exclude_backbone=False):
        self.name = name
        self.exclude_backbone = exclude_backbone
        self.file_path = file_path
        self.save_file_info()
        self.residues = {}
        self.centroid_cliques = None
        self.centroid_graph = None
        self.centroid_clique_distances = None
        self.centroid_clique_frequency = None
        atom_count = 0
        res_count = -1
        prev_res = -1
        with open(self.file_path) as pdb_file:
            for line in pdb_file:
                if line[0:4] == "ATOM":
                    res_name = line[17:20].strip(" ")
                    res_id = int(line[22:26].strip(" "))
                    if prev_res != res_id:
                        prev_res = res_id
                        res_count += 1
                    # atom_id = int(line[6:11].strip(" "))
                    res_id = res_count
                    atom_id = atom_count
                    atom_count += 1
                    coordx = float(line[30:38].strip(" "))
                    coordy = float(line[38:46].strip(" "))
                    coordz = float(line[46:54].strip(" "))
                    symbol = line[76:78].strip(" ")
                    atom_name = line[12:16].strip(" ")
                    coords = (coordx, coordy, coordz)
                    atm = atom.Atom(symbol, atom_name, atom_id, coords)
                    if self.residues.get(res_id) is None:
                        self.residues[res_id] = residue.Residue(res_name, res_id, [atm])
                    else:
                        self.residues[res_id].add_atom(atm)
        for i in self.residues:
            self.residues[i].update_COM(exclude_backbone=self.exclude_backbone)

    def save_file_info(self):
        if not os.path.exists(r"C:\Users\aprak\PycharmProjects\TopProPack_v2\top_pro_pack_logs\{}".format(self.name)):
            os.makedirs(r"C:\Users\aprak\PycharmProjects\TopProPack_v2\top_pro_pack_logs\{}".format(self.name))

        basic_info = open("basic-info.txt", "w")

    def get_name(self):
        return self.name

    def get_file_path(self):
        return self.file_path

    def get_residues(self):
        return self.residues

    def generate_centroid_cliques(self, distance_cutoff=6):
        centroids = []
        centroid_res = {}
        for i in self.residues:
            centroids.append(self.residues[i].get_centroid())
            centroid_res[self.residues[i].get_centroid()] = self.residues[i]
        centroids = np.array(centroids)
        tri = scipy.spatial.qhull.Delaunay(centroids)
        edges = []
        for n in tri.simplices:
            edge = sorted([n[0], n[1]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[0], n[2]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[0], n[3]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[1], n[2]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[1], n[3]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[2], n[3]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
        graph = nx.Graph(edges)
        self.centroid_graph = graph
        self.centroid_cliques = list(nx.find_cliques(graph))
        for i in range(len(self.centroid_cliques)):
            for j in range(len(self.centroid_cliques[i])):
                self.centroid_cliques[i][j] = centroid_res[tuple(list(centroids[self.centroid_cliques[i][j]]))]
        self.centroid_cliques = np.array(self.centroid_cliques)

    def get_clique_frequency(self):
        if self.centroid_clique_frequency is not None:
            return self.centroid_clique_frequency
        if self.centroid_cliques is None:
            self.generate_centroid_cliques()
        freq_arr = [0, 0, 0, 0, 0, 0, 0]
        for i in self.centroid_cliques:
            freq_arr[len(i)] += 1
        self.centroid_clique_frequency = freq_arr
        return freq_arr

    def get_centroid_clique_distances(self):
        if self.centroid_clique_distances is not None:
            return self.centroid_clique_distances
        distances = []
        for i in range(len(self.centroid_cliques)):
            clique = self.centroid_cliques[i]
            coords = []
            for j in clique:
                coords.append(j.get_centroid())
            for x in range(len(coords) - 1):
                for y in range(x + 1, len(coords)):
                    d = get_dist(coords[x], coords[y])
                    distances.append(d)
        self.centroid_clique_distances = distances
        return self.centroid_clique_distances

    def get_heatmap_data_centroid(self):
        arr = [
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        ]
        ref = {
            "ALA": 0,
            "ARG": 1,
            "ASN": 2,
            "ASP": 3,
            "CYS": 4,
            "GLN": 5,
            "GLU": 6,
            "GLY": 7,
            "HIS": 8,
            "ILE": 9,
            "LEU": 10,
            "LYS": 11,
            "MET": 12,
            "PHE": 13,
            "PRO": 14,
            "SER": 15,
            "THR": 16,
            "TRP": 17,
            "TYR": 18,
            "VAL": 19
        }
        for clique in self.centroid_cliques:
            for i in range(len(clique)):
                for j in range(i + 1, len(clique)):
                    if ref[clique[i].get_name()] == ref[clique[j].get_name()]:
                        arr[ref[clique[i].get_name()]][ref[clique[j].get_name()]] += 1
                    else:
                        arr[ref[clique[i].get_name()]][ref[clique[j].get_name()]] += 1
                        arr[ref[clique[j].get_name()]][ref[clique[i].get_name()]] += 1
        return np.array(arr)
os.chdir(r"C:\Users\aprak\PycharmProjects\TopProPack_v2")
file_path = r"C:\alpha\7prc.pdb"
P = CentroidProtein("7prc", file_path)
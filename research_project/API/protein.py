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



#generate clique codes
#sort clique codes
#find included bench clique codes in new cliques of any type

def binary_search(arr, val):
    start = 0
    end = len(arr)-1
    while True:
        mid = int((start+end)/2)
        #print(start, end)
        if start >= end and arr[mid] != val: return -1
        if arr[mid] == val: return mid
        elif val > arr[mid]: start = mid+1
        else: end = mid-1

###############################################################################################
#TODOs:
#   make option to exclude/include backbone atoms(50% completion)
#   csv data logging
#   matplotlib/excel/google_sheets graphing methods
#   other data analysis methods[*priority*: benchmark to new analysis]
###############################################################################################

def get_dist(coord1, coord2):
    return math.sqrt((coord1[0]-coord2[0])**2+(coord1[1]-coord2[1])**2+(coord1[2]-coord2[2])**2)

class Protein:
    def read_general_info(self, file):
        gen_info_file = open(file, "r")
        data = gen_info_file.readlines()
        name = data[0][6:len(data[0])-1]
        file_path = data[1][10:len(data[1])-1]
        exclude_backbone = data[2][17:len(data[2])-1]
        num_of_residues = data[3][14:len(data[3])-1]
        #print(name, file_path, exclude_backbone, num_of_residues)
        gen_info_file.close()
        return name, file_path, exclude_backbone, num_of_residues

    def write_general_info(self, file):
        general_info_file = open(file, "w")
        general_info_file.write("name: {}\n".format(self.name))
        general_info_file.write("file_path: {}\n".format(self.file_path))
        var = "False"
        if self.exclude_backbone: var = "True"
        general_info_file.write("exclude_backbone: {}\n".format(var))
        general_info_file.write("#_of_residues: {}\n".format(len(self.residues.keys())))
        general_info_file.close()
        return

    def read_graph_info(self, file):
        graph_info_file = open(file, "r")
        edges = []
        for line in graph_info_file:
            line = line[0:len(line)-1].split(",")
            edges.append(list(map(int, line)))
        graph = nx.Graph(edges)
        graph_info_file.close()
        return graph

    def write_graph_info(self, file, edges):
        graph_info_file = open(file, "w")
        for edge in edges:
            graph_info_file.write("{},{}\n".format(edge[0], edge[1]))
        graph_info_file.close()
        return

    def read_residue_info(self, file):
        residue_info_file = open(file, "r")
        curr_index = 0
        residues = {}
        for line in residue_info_file.readlines():
            line = line[0:len(line)-1]
            line = line.split(" ")
            atm_count = 0
            res_name = None
            res_id = None
            if line[0] == "RES":                
                res_name = line[1]
                res_id = int(line[2])
                atm_count = int(line[3])
                curr_index = res_id
                residues[curr_index] = residue.Residue(res_name, res_id, [])
            else:
                atom_name = line[1]
                symbol = line[2]
                atom_id = int(line[3])
                coords = tuple(line[4])
                atm = atom.Atom(symbol, atom_name, atom_id, coords)
                residues[curr_index].add_atom(atm)
        return residues

    def read_atom_cliques(self, file):
        atom_cliques_file = open(file, "r")
        data = atom_cliques_file.readlines()
        atom_cliques = []
        for line in data:
            #line = line[0:len(line)-2]
            atom_cliques.append(list(map(int, line[0:len(line)-1].split(","))))
        atom_cliques = np.array(atom_cliques)
        #print(atom_cliques)
        atom_cliques_file.close()
        return atom_cliques

    def read_centroid_cliques(self, file):
        centroid_cliques_file = open(file, "r")
        data = centroid_cliques_file.readlines()
        centroid_cliques = []
        for line in data:
            #line = line[0:len(line)-2]
            centroid_cliques.append(list(map(int, line[0:len(line)-1].split(","))))
        centroid_cliques = np.array(centroid_cliques)
        #print(atom_cliques)
        centroid_cliques_file.close()
        return centroid_cliques

    def write_atom_cliques(self, file):
        atom_cliques_file = open(file, "w")
        for i in range(len(self.atom_cliques)):
            text = ""
            for j in range(len(self.atom_cliques[i])):
                text += str(self.atom_cliques[i][j].get_resid())
                text += ","
            text = text[0:len(text)-1]
            text += "\n"
            atom_cliques_file.write(text)
        atom_cliques_file.close()
        return

    def write_centroid_cliques(self, file):
        centroid_cliques_file = open(file, "w")
        for i in range(len(self.centroid_cliques)):
            text = ""
            for j in range(len(self.centroid_cliques[i])):
                text += str(self.centroid_cliques[i][j].get_resid())
                text += ","
            text = text[0:len(text)-1]
            text += "\n"
            centroid_cliques_file.write(text)
        centroid_cliques_file.close()
        return


            
    def __init__(self, name, file_path, exclude_backbone=False, load_files=None):
        if load_files == None:
            self.name = name
            self.exclude_backbone = exclude_backbone
            self.file_path = file_path
            self.residues = {}
            self.atom_cliques = None
            self.centroid_cliques = None
            self.atom_graph = None
            self.centroid_graph = None
            #self.distance_cutoff = 6
            self.centroid_clique_frequency = None
            self.atom_clique_frequency = None
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
                        #atom_id = int(line[6:11].strip(" "))
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
                        if self.residues.get(res_id) == None:
                            self.residues[res_id] = residue.Residue(res_name, res_id, [atm])
                        else: self.residues[res_id].add_atom(atm)
            for i in self.residues:
                self.residues[i].update_COM(exclude_backbone=self.exclude_backbone)

            os.chdir("C:\\top_pro_pack_logs\\proteins")
            if not os.path.exists(self.name): os.mkdir(self.name)
            os.chdir(self.name)
            self.write_general_info("general_info.txt")
        else:
            gen_info_file = load_files[0]
            res_info_file = load_files[1]
            atom_cliques_file = load_files[2]
            centroid_cliques_file = load_files[3]
            atom_graph_file = load_files[4]
            centroid_graph_file = load_files[5]
            name, file_path, exclude_backbone, num_of_residues = self.read_general_info(gen_info_file)
            self.name = name
            self.file_path = file_path
            self.exclude_backbone = bool(exclude_backbone.strip(" "))
            self.residues = self.read_residue_info(res_info_file)
            self.atom_cliques = self.read_atom_cliques(atom_cliques_file)
            self.centroid_cliques = self.read_centroid_cliques(centroid_cliques_file)
            self.atom_graph = self.read_graph_info(atom_graph_file)
            self.centroid_graph = self.read_graph_info(centroid_graph_file)
            self.centroid_clique_frequency = None
            self.atom_clique_frequency = None
            #for i in self.residues:
                #self.residues[i].update_COM(exclude_backbone=self.exclude_backbone)
        #self.read_general_info("general_info.txt")
        #print(self.read_residue_info("residue_info.txt"))
        #self.read_atom_cliques("atom_cliques.txt")
        #self.read_centroid_cliques("centroid_cliques.txt")


        #RES GLY 5
        #ATM name symbol id coords

        #centroid clique graph.txt/atom_clique_graph.txt
        

        residues_file = open("residue_info.txt", "w")
        for key in self.residues:
            name = self.residues[key].get_name()
            resid = self.residues[key].get_resid()
            atm_count = len(self.residues[key].get_atoms())
            residues_file.write("RES {} {} {}\n".format(name, resid, atm_count))
            for atm in self.residues[key].get_atoms():
                atm_name = atm.get_name()
                atm_symbol = atm.get_symbol()
                atm_id = atm.get_atomid()
                atm_coords = atm.get_coords()
                text = "ATM {} {} {} {}\n".format(atm_name, atm_symbol, str(atm_id), str(atm_coords))
                residues_file.write(text)
        residues_file.close()
                
        
    def convert_to_code(self, clique):
        clique = list(clique)
        clique.sort()
        s = ""
        l = list(map(str, clique))
        for i in l: s+= i
        return int(s)

    def get_sorted_clique_codes(self, cliques):
        codes = []
        for i in range(len(cliques)):
            #print(cliques[i])
            codes.append(self.convert_to_code(cliques[i]))
        codes.sort()
        return codes

    def get_resindex_resid_hash(self, bench_index_file):
        resindex_resid = {}
        with open(bench_index_file) as file:
            for line in file:
                a = line.strip(" \n").split(" ")
                resindex_resid[int(a[0])] = int(a[1])
        return resindex_resid

    def convert_resindex_resid(self, bench_cliques, resindex_resid):
        for i in range(len(bench_cliques)):
            for j in range(len(bench_cliques[i])):
                bench_cliques[i][j] = resindex_resid[bench_cliques[i][j]]
        return bench_cliques

    def get_bench_cliques(self, bench_file):
        cliques = []
        with open(bench_file) as file:
            for line in file:
                pos = line.find("|", line.find("|")+1, len(line)-1)
                clique = list(map(int, line[pos+2:].strip(" \n").split(" ")))
                cliques.append(clique)
        return cliques

    def get_included_bench_cliques(self, codes_bench, codes_atom):
        included = []
        for i in codes_bench:
            x = binary_search(codes_atom, i)
            if x == -1: continue #print("\n\n\n\n\n", x, "\n\n\n\n\n")
            else:
                included.append(codes_atom[x])
                #print("\n\n\n\n\n", x, "\n\n\n\n\n")
        return included
    
    def get_name(self):
        return self.name

    def get_file_path(self):
        return self.file_path

    def get_residues(self):
        return self.residues

    def get_cliques(self, clique_type, distance_cutoff=6): #come back to exclude_backbone as a parameter later
        if clique_type == "centroid":
            if self.centroid_cliques == None:
                self.generate_cliques("centroid", distance_cutoff=distance_cutoff)
            return self.centroid_cliques
        elif clique_type == "atom":
            if self.atom_cliques == None:
                self.generate_cliques("atom", distance_cutoff=distance_cutoff)
            return self.atom_cliques
        else: raise Exception("invalid clique type")
            
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
        self.write_graph_info("centroid_graph_info.txt", edges)
        graph = nx.Graph(edges)
        self.centroid_cliques = list(nx.find_cliques(graph))
        for i in range(len(self.centroid_cliques)):
            for j in range(len(self.centroid_cliques[i])):
                self.centroid_cliques[i][j] = centroid_res[tuple(list(centroids[self.centroid_cliques[i][j]]))]
        self.centroid_cliques = np.array(self.centroid_cliques)
        self.write_centroid_cliques("centroid_cliques.txt")#TODO: integerate write command
            

    def generate_atom_cliques(self, distance_cutoff=6):
        coords = []
        coords_resid = {}
        for i in self.residues:
            for j in self.residues[i].get_atoms():
                if self.residues[i].get_name() != "GLY" and self.exclude_backbone:
                    if j.get_name() != "CA" or j.get_name() != "C" or j.get_name() != "N" or j.get_name() != "O":
                        coords.append(j.get_coords())
                        coords_resid[j.get_coords()] = i
                else:
                    coords.append(j.get_coords())
                    coords_resid[j.get_coords()] = i
        coords_array = np.array(coords)
        del coords
        tri = scipy.spatial.qhull.Delaunay(coords_array)
        edges = []
        for n in tri.simplices:
            edge = sorted([n[0], n[1]])
            if get_dist(self.residues[coords_resid[tuple(list(coords_array[edge[0]]))]].get_centroid(), self.residues[coords_resid[tuple(list(coords_array[edge[1]]))]].get_centroid()) <= distance_cutoff:
                edges.append((edge[0], edge[1]))
            edge = sorted([n[0], n[2]])
            if get_dist(self.residues[coords_resid[tuple(list(coords_array[edge[0]]))]].get_centroid(), self.residues[coords_resid[tuple(list(coords_array[edge[1]]))]].get_centroid()) <= distance_cutoff:
                edges.append((edge[0], edge[1]))
            edge = sorted([n[0], n[3]])
            if get_dist(self.residues[coords_resid[tuple(list(coords_array[edge[0]]))]].get_centroid(), self.residues[coords_resid[tuple(list(coords_array[edge[1]]))]].get_centroid()) <= distance_cutoff:
                edges.append((edge[0], edge[1]))
            edge = sorted([n[1], n[2]])
            if get_dist(self.residues[coords_resid[tuple(list(coords_array[edge[0]]))]].get_centroid(), self.residues[coords_resid[tuple(list(coords_array[edge[1]]))]].get_centroid()) <= distance_cutoff:
                edges.append((edge[0], edge[1]))
            edge = sorted([n[1], n[3]])
            if get_dist(self.residues[coords_resid[tuple(list(coords_array[edge[0]]))]].get_centroid(), self.residues[coords_resid[tuple(list(coords_array[edge[1]]))]].get_centroid()) <= distance_cutoff:
                edges.append((edge[0], edge[1]))
            edge = sorted([n[2], n[3]])
            if get_dist(self.residues[coords_resid[tuple(list(coords_array[edge[0]]))]].get_centroid(), self.residues[coords_resid[tuple(list(coords_array[edge[1]]))]].get_centroid()) <= distance_cutoff:
                edges.append((edge[0], edge[1]))
        self.write_graph_info("atom_graph_info.txt", edges)
        graph = nx.Graph(edges)
        self.atom_cliques = list(nx.find_cliques(graph))
        temp_arr = []
        for i in range(len(self.atom_cliques)):
            for j in range(len(self.atom_cliques[i])):
                self.atom_cliques[i][j] = coords_resid[tuple(list(coords_array[self.atom_cliques[i][j]]))]
            self.atom_cliques[i] = list(set(self.atom_cliques[i]))
            self.atom_cliques[i].sort()
            if self.atom_cliques[i] not in temp_arr:
                temp_arr.append(self.atom_cliques[i])
        self.atom_cliques = np.array(temp_arr)
        for i in range(len(self.atom_cliques)):
            for j in range(len(self.atom_cliques[i])):
                self.atom_cliques[i][j] = self.residues[self.atom_cliques[i][j]]
        self.atom_cliques = np.array(self.atom_cliques)
        self.write_atom_cliques("atom_cliques.txt")

    def generate_cliques(self, clique_type, distance_cutoff=6):
        if clique_type == "centroid":
            self.generate_centroid_cliques(distance_cutoff=distance_cutoff)
        elif clique_type == "atom":
            self.generate_atom_cliques(distance_cutoff=distance_cutoff)
        else: raise Exception("Invalid clique type")

    def getMaxMinDistance(self, coords):
        max_dist = 0
        min_dist = 10000
        if len(coords) == 1: return 0, 0
        for i in range(len(coords)-1):
            for j in range(i+1, len(coords)):
                dist = get_dist(coords[i], coords[j])
                max_dist = max([max_dist, dist])
                min_dist = min([min_dist, dist])
        return max_dist, min_dist

    def distance_analysis(self, clique_type):
        if clique_type == "centroid":
            print("Min/Max CENTROID CLIQUE STATS")
            clique_max_sum = 0
            clique_min_sum = 0
            count = 0
            for i in self.centroid_cliques:
                coords = []
                for res in i:
                    coords.append(res.get_centroid())
                clique_max, clique_min = self.getMaxMinDistance(coords)
                if clique_max != 0 and clique_min != 0:
                    clique_max_sum += clique_max
                    clique_min_sum += clique_min
                    count += 1
            clique_max_avg = float(clique_max_sum)/count
            clique_min_avg = float(clique_min_sum)/count
                #print("min_dist: {} | max_dist: {}".format(clique_min, clique_max))
            print("avg_min_dist: {} | avg_max_dist: {}".format(clique_min_avg, clique_max_avg))
        elif clique_type == "atom":
            print("Min/Max ATOM CLIQUE STATS")
            clique_max_sum = 0
            clique_min_sum = 0
            count = 0
            for i in self.atom_cliques:
                coords = []
                for res in i:
                    coords.append(res.get_centroid())
                clique_max, clique_min = self.getMaxMinDistance(coords)
                if clique_max != 0 and clique_min != 0:
                    clique_max_sum += clique_max
                    clique_min_sum += clique_min
                    count += 1
            clique_max_avg = float(clique_max_sum)/count
            clique_min_avg = float(clique_min_sum)/count
                #print("min_dist: {} | max_dist: {}".format(clique_min, clique_max))
            print("avg_min_dist: {} | avg_max_dist: {}".format(clique_min_avg, clique_max_avg))
        else: raise Exception("Invalid clique type")

    def freq_analysis(self, clique_type, bench_cliques_file=None):
        freq_arr = [0,0,0,0,0,0,0]
        if clique_type == "centroid":
            for i in self.centroid_cliques:
                freq_arr[len(i)] += 1
            self.centroid_clique_frequency = freq_arr
        elif clique_type == "atom":
            for i in self.atom_cliques:
                freq_arr[len(i)] += 1
            self.atom_clique_frequency = freq_arr
        else: raise Exception("Invalid clique type")
        bench_freq_arr = [0,0,0,0,0,0,0]
        if bench_cliques_file != None:
            bench_cliques = self.get_bench_cliques(bench_cliques_file)
            for i in bench_cliques:
                bench_freq_arr[len(i)] += 1
            print("{} v. BENCH CLIQUE SIZE FREQS".format(clique_type))
            print("{}: {}".format(clique_type, freq_arr[1:]))
            print("bench: {}".format(bench_freq_arr[1:]))
        else:
            print("{} CLIQUE SIZE FREQS".format(clique_type))
            print("{}: {}".format(clique_type, freq_arr[1:]))
        #return freq_arr

    def get_included_bench_cliques(self, codes_bench, codes_new):
        included = []
        for i in codes_bench:
            x = binary_search(codes_new, i)
            if x != -1:
                included.append(codes_new[x])
        return included

    def bench_to_new_analysis(self, clique_type, bench_cliques_file, bench_index_file):
        resindex_resid = self.get_resindex_resid_hash(bench_index_file)
        cliques = []
        new_codes = None
        if clique_type == "centroid":
            for i in self.centroid_cliques:
                clique = []
                for j in i:
                    clique.append(j.get_resid())
                cliques.append(clique)
            cliques = np.array(cliques)
            new_codes = self.get_sorted_clique_codes(cliques)
        elif clique_type == "atom":
            for i in self.atom_cliques:
                clique = []
                for j in i:
                    clique.append(j.get_resid())
                cliques.append(clique)
            cliques = np.array(cliques)
            new_codes = self.get_sorted_clique_codes(cliques)
        else: raise Exception("Invalid clique type")
        bench_cliques = self.get_bench_cliques(bench_cliques_file)
        bench_cliques = self.convert_resindex_resid(bench_cliques, resindex_resid)
        bench_codes = self.get_sorted_clique_codes(bench_cliques)
        included = self.get_included_bench_cliques(bench_codes, new_codes)
        #print(included)
        print("{} CLIQUES BENCH INCLUSION STATS".format(clique_type)) 
        print("% included: {}".format((float(len(included))/len(bench_codes))*100))
        print("New: {} | Bench: {} | Included: {}".format(len(new_codes), len(bench_codes), len(included)))
        
test_bench_cliques_file = "C:\\pdb_data\\domains\\f1a4fa_.nomc.cliques"
test_bench_index_file = "C:\\pdb_data\\domains\\f1a4fa_.index"
#file = "C:\\pdb_data\\domains\\f1a4fa_.pdb"
file = "C:\\alpha\\4quv.pdb"




protein = Protein("f1a4fa_", file, True, load_files=["C:\\top_pro_pack_logs\\proteins\\f1a4fa_\\general_info.txt", "C:\\top_pro_pack_logs\\proteins\\f1a4fa_\\residue_info.txt", "C:\\top_pro_pack_logs\\proteins\\f1a4fa_\\atom_cliques.txt", "C:\\top_pro_pack_logs\\proteins\\f1a4fa_\\centroid_cliques.txt"])
#protein.get_cliques("centroid")
#protein.get_cliques("atom")
#protein.distance_analysis("atom")
#protein.distance_analysis("centroid")
#protein.freq_analysis("atom", test_bench_cliques_file)
#protein.bench_to_new_analysis("atom", test_bench_cliques_file, test_bench_index_file)

#print(protein.centroid_cliques[0][1].get_centroid())


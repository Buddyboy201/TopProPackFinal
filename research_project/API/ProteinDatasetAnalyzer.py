import atom
import residue
import protein
import os
import sys
import numpy as np
import scipy.spatial
import networkx as nx
import csv
import math
import time
import matplotlib.pyplot as plt

class ProteinDatasetAnalyzer:
    def __init__(self, proteins_names=[], protein_file_paths=[]):
        self.protein_names = proteins_names
        self.proteins = []
        for i in range(len(self.protein_names)):
            self.proteins.append(protein.Protein(self.protein_names[i], protein_file_paths[i]))


    def add_protein(self, new_protein_name, new_protein_file_path):
        new_protein = Protein(new_protein_name, new_protein_file_path)
        self.proteins.append(new_protein)

    def get_proteins(self):
        return self.proteins

    def get_protein(self, name):
        x = protein.binary_search(self.protein_names, name)
        if x == -1: return -1
        return self.proteins[x]

    def get_bench_protein(self, name, bench_protein_data):
        names = []
        for i in bench_protein_data:
            names.append(i[0])
        x = protein.binary_search(names, name)
        return x

    def generate_all_cliques(self):
        for i in range(len(self.proteins)):
            self.proteins[i].get_cliques("centroid")
            self.proteins[i].get_cliques("atom")

    def generate_cliques(self, clique_type):
        for i in range(len(self.proteins)):
            self.proteins[i].get_cliques(clique_type)
            
        

    #protein bench data format:
    # [ [name (same as in private field), clique file path, index file path], ...]

    def get_general_clique_stats(self, bench_names=None, bench_clique_files=None, bench_index_files=None):     
        self.generate_all_cliques()
        #self.generate_cliques("centroid")
        for P in self.proteins:
            print("\n{} GENERAL CLIQUE STATS:\n".format(P.get_name()))
            #P.get_cliques("centroid")
            #P.get_cliques("atom") 
            if bench_protein_data != None:
                for N, C, I in bench_protein_data:
                    x = self.get_bench_protein()
                    protein.freq_analysis("atom", test_bench_cliques_file)
                    protein.bench_to_new_analysis("atom", test_bench_cliques_file, test_bench_index_file)
            else:
                P.distance_analysis("atom")
                P.distance_analysis("centroid")
                P.freq_analysis("atom")
                P.freq_analysis("centroid")
            time.sleep(1)
            print("\n")
                

names = []
files = []
os.chdir("C:\\pdb_data\\domains")
count = 0
for i in os.listdir():
    if i[len(i)-4:] == ".pdb":
        names.append(i[0:len(i)-4])
        files.append("C:\\pdb_data\\domains\\" + i)
        count += 1
    if count == 10: break
#print(names)
#print(files)
#analyzer = ProteinDatasetAnalyzer(["7prc", "6rfc", "6rfb"], ["C:\\alpha\\7prc.pdb", "C:\\alpha\\6rfc.pdb", "C:\\alpha\\6rfb.pdb"])
analyzer = ProteinDatasetAnalyzer(names, files)
analyzer.get_general_clique_stats()
#analyzer.generate_all_cliques()

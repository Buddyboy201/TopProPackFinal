import os
import sys
import parse_pdb_and_generate_cliques_v3 as atom_clique_gen
import numpy as np
import xlsxwriter as xw
import excel_bar_test as eb_test
from time import *

def convert_to_code(clique):
    clique = list(clique)
    clique.sort()
    s = ""
    l = list(map(str, clique))
    for i in l: s += i
    return int(s)

def get_sorted_clique_codes(cliques):
    codes = []
    for i in cliques:
        codes.append(convert_to_code(i))
    codes.sort()
    return codes

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

def get_resindex_resid_hash(bench_index_file):
    resindex_resid = {}
    with open(bench_index_file) as file:
        for line in file:
            a = line.strip(" \n").split(" ")
            resindex_resid[int(a[0])] = int(a[1])
    return resindex_resid

def convert_resindex_resid(bench_cliques, resindex_resid):
    for i in range(len(bench_cliques)):
        for j in range(len(bench_cliques[i])):
            bench_cliques[i][j] = resindex_resid[bench_cliques[i][j]]
    return bench_cliques

def get_bench_cliques(bench_file):
    cliques = []
    with open(bench_file) as file:
        for line in file:
            pos = line.find("|", line.find("|")+1, len(line)-1)
            clique = list(map(int, line[pos+2:].strip(" \n").split(" ")))
            cliques.append(clique)
    return cliques

def get_included_bench_cliques(codes_bench, codes_atom):
    included = []
    for i in codes_bench:
        x = binary_search(codes_atom, i)
        if x == -1: continue #print("\n\n\n\n\n", x, "\n\n\n\n\n")
        else:
            included.append(codes_atom[x])
            #print("\n\n\n\n\n", x, "\n\n\n\n\n")
    return included

def get_included_atom_cliques(codes_atom, included_bench):
    included = []
    for i in codes_atom:
        if binary_search(included_bench, i) != -1:
            pass
def main():
    test_file = "c:\\alpha\\4quv.pdb"
    
    #print(binary_search(codes, 90115119130))
    os.chdir("c:\\pdb_data\\domains")
    workbook = xw.Workbook("c:\\excel_test\\bench_atom_include_data.xlsx")
    files = os.listdir()
    for i in range(0, len(files)-12, 12):
        pdb_file = files[i+8]
        file_name = pdb_file[0:len(pdb_file)-4]
        cliques_file = files[i+4]
        index_file = files[i+3]
        coords_array, coords_resid, resid_res, cliques_atom = atom_clique_gen.generate_hashmaps_and_cliques(pdb_file)
        cliques_atom = atom_clique_gen.update_cliques_2(cliques_atom, coords_array, coords_resid)
        codes_atom = np.array(get_sorted_clique_codes(cliques_atom))
        bench_cliques = get_bench_cliques(cliques_file)
        resindex_resid = get_resindex_resid_hash(index_file)
        bench_cliques = convert_resindex_resid(bench_cliques, resindex_resid)
        codes_bench = np.array(get_sorted_clique_codes(bench_cliques))
        included = np.array(get_included_bench_cliques(codes_bench, codes_atom))
        print("\n\n\n")
        print(len(codes_bench), len(codes_atom), len(included))
        print("\n\n\n")
        #sleep(3)
        eb_test.writeXLSXCode(workbook, file_name, [len(codes_bench)], [len(codes_atom)])
    workbook.close()

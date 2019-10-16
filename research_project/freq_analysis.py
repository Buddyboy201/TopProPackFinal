import os
import sys
import numpy as np
import scipy.spatial
import networkx as nx
import csv
import parse_pdb_and_generate_cliques_v3 as atom_clique_gen
import centroid_based_clique_generator as cen_clique_gen
import xlsxwriter as xw
import excel_bar_test as eb_test

def get_benchmark_clique_stats(cliques_file):
    freq_arr = [0,0,0,0,0,0,0,0]
    cliques = []
    with open(cliques_file) as file:
        for line in file:   
            one = line.find("|")
            two = line.find("|", one+1, len(line)-1)
            val = int(line[one+1:two-1].strip(" "))
            freq_arr[val] += 1
    return np.array(freq_arr)

def get_clique_stats(cliques):
    freq_arr = [0,0,0,0,0,0,0,0]
    for i in cliques:
        try: freq_arr[len(i)] += 1
        except: raise Exception(len(i), i)
    return np.array(freq_arr)

def format_results(name, freq_arr_bench, freq_arr_atom, freq_arr_centroid):
    print("{} RESULTS:".format(name))
    print("\tBENCH v. ATOM v. CENTROID")
    print("TOTALS: {}\t{}\t{}".format(sum(list(freq_arr_bench)), sum(list(freq_arr_atom)), sum(list(freq_arr_centroid))))
    for i in range(1,8):
        print("{}: \t{}\t{}\t{}".format(i, freq_arr_bench[i], freq_arr_atom[i], freq_arr_centroid[i]))


def test():
    file = "c:\\alpha\\6rf6.pdb"
    os.chdir("c:\\alpha")
    files = os.listdir()
    #for file in files
    cliques_centroid = cen_clique_gen.generate_centroid_cliques(file)
    count = 0
    for i in cliques_centroid:
        if len(i) == 4:
            print(i)
            count += 1
        if count > 10: break
    freq_arr_centroid = get_clique_stats(cliques_centroid)
    format_results(file, [0,0,0,0,0,0,0,0], [0,0,0,0,0,0,0,0], freq_arr_centroid)

    
def main():
    os.chdir("c:\\pdb_data\\domains")
    files = os.listdir()
    workbook = xw.Workbook("c:\\excel_test\\bench_new_freq_data.xlsx")
    total_freq_arr_bench = [0,0,0,0,0,0,0,0]
    total_freq_arr_atom = [0,0,0,0,0,0,0,0]
    total_freq_arr_centroid = [0,0,0,0,0,0,0,0]
    for i in range(0, len(files)-12, 12):
        pdb_file = files[i+8]
        file_name = pdb_file[0:len(pdb_file)-4]
        cliques_file = files[i+4]
        coords_array, coords_resid, resid_res, cliques_atom = atom_clique_gen.generate_hashmaps_and_cliques(pdb_file)
        cliques_atom = atom_clique_gen.update_cliques(cliques_atom, coords_array, coords_resid, resid_res)
        cliques_centroid = cen_clique_gen.generate_centroid_cliques(pdb_file)
        freq_arr_bench = get_benchmark_clique_stats(cliques_file)
        freq_arr_atom = get_clique_stats(cliques_atom)
        freq_arr_centroid = get_clique_stats(cliques_centroid)
        for i in range(len(freq_arr_bench)):
            total_freq_arr_bench[i] += freq_arr_bench[i]
            total_freq_arr_atom[i] += freq_arr_atom[i]
            total_freq_arr_centroid[i] += freq_arr_centroid[i]
        format_results(file_name, freq_arr_bench, freq_arr_atom, freq_arr_centroid)
        eb_test.writeXLSXProtein(workbook, file_name, list(freq_arr_bench)[1:], list(freq_arr_atom)[1:], list(freq_arr_centroid)[1:])
    workbook_total = xw.Workbook("c:\\excel_test\\bench_new_freq_data_totals.xlsx")
    eb_test.writeXLSXProtein(workbook_total, "totals", list(total_freq_arr_bench)[1:], list(total_freq_arr_atom)[1:], list(total_freq_arr_centroid)[1:])
    workbook.close()
    workbook_total.close()

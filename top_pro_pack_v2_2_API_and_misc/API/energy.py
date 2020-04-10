import os
import sys
import centroid_protein
import numpy as np
import pandas as pd
import statistics as stat
import math
import energy

class Energy:
    def __init__(self):
        self.STATIC_TOTAL_PAIRS_TABLE = pd.DataFrame(
            np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]] * 20).astype("int32"))
        self.STATIC_EPAIR_TABLE = pd.DataFrame(
            np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]] * 20).astype("float64"))
        self.up_to_date = False
        self.ref = {
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

    def initialize_static_total_pairs_table(self):
        pass

    def initialize_static_epair_table(self):
        pass

    def update_static_total_pairs_table(self, protein_pairs_matrix):
        self.STATIC_TOTAL_PAIRS_TABLE = pd.DataFrame(np.add(self.STATIC_TOTAL_PAIRS_TABLE.values.astype("int32"),
                                                            protein_pairs_matrix.astype("int32")).astype('int32'))

    def update_epair_values(self):
        epair_heat_map = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]] * 20).astype("float64")
        M_const = self.get_M(self.STATIC_TOTAL_PAIRS_TABLE)
        for aa_i in self.ref:
            for aa_j in self.ref:
                e_pair = self.get_epair(self.ref[aa_i], self.ref[aa_j], M_const, self.STATIC_TOTAL_PAIRS_TABLE)
                if aa_i == aa_j:
                    epair_heat_map[self.ref[aa_i]][self.ref[aa_j]] = e_pair
                else:
                    epair_heat_map[self.ref[aa_i]][self.ref[aa_j]] = e_pair
                    epair_heat_map[self.ref[aa_j]][self.ref[aa_i]] = e_pair
        self.STATIC_EPAIR_TABLE = pd.DataFrame(epair_heat_map)

    def get_static_total_pairs_table(self):
        return self.STATIC_TOTAL_PAIRS_TABLE

    def get_static_epair_table(self):
        return self.STATIC_EPAIR_TABLE

    def get_epair(self, aa_i, aa_j):
        return self.STATIC_EPAIR_TABLE.iloc[aa_i, aa_j]

    def get_M(self, aa_heat_map):
        total = 0
        for i in range(20):
            total += aa_heat_map.iloc[:, i].sum()
        return total

    def M_single(self, i, aa_heat_map):
        return aa_heat_map.iloc[:, i].sum()

    def M_E(self, i, j, M_const, aa_heat_map):
        return (self.M_single(i, aa_heat_map) * self.M_single(j, aa_heat_map)) / (M_const ** 2)

    def M_pair(self, i, j, M_const, aa_heat_map):
        return (aa_heat_map.iloc[j, i]) / M_const

    def get_epair(self, i, j, M_const, aa_heat_map):
        return -math.log(self.M_pair(i, j, M_const, aa_heat_map) / self.M_E(i, j, M_const, aa_heat_map), math.e)





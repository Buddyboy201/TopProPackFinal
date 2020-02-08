import os
import sys
import centroid_protein
import numpy as np
import pandas as pd
import statistics as stat

#TODO: implement init/update system to epair/pair_counts table
#TODO: implement changes to heatmap as Vladimir requested
#TODO: implement data logging system

class Energy:
    def __init__(self):
        self.STATIC_TOTAL_PAIRS_TABLE = pd.DataFrame(
            np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]] * 20).astype("int32"))
        self.STATIC_EPAIR_TABLE = pd.DataFrame(
            np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]] * 20).astype("float64"))
        self.up_to_date = False

    def initialize_static_total_pairs_table(self):
        pass

    def initialize_static_epair_table(self):
        pass

    def update_static_total_pairs_table(self, protein_pairs_matrix):
        self.STATIC_TOTAL_PAIRS_TABLE = pd.DataFrame(np.add(self.STATIC_TOTAL_PAIRS_TABLE.values.astype("int32"),
                                                            protein_pairs_matrix.astype("int32")).astype('int32'))

    def update_epair_values(self):
        pass

    def get_static_total_pairs_table(self):
        return self.STATIC_TOTAL_PAIRS_TABLE

    def get_static_epair_table(self):
        return self.STATIC_EPAIR_TABLE

    def get_epair(self, aa_i, aa_j):
        return self.STATIC_EPAIR_TABLE.iloc[aa_i, aa_j]



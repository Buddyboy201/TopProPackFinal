import centroid_protein
import atom
import residue
import energy
import os
import sys

class TPP_Engine:
    def __init__(self, exclude_backbone=False):
        os.chdir(r"C:\Users\aprak\PycharmProjects\TopProPack_v2")
        self.STATIC_TOTAL_PAIRS_TABLE = None
        self.STATIC_EPAIR_TABLE = None
        self.exclude_backbone = exclude_backbone

    def initialize_static_total_pairs_table(self):
        pass

    def initialize_static_epair_table(self):
        pass





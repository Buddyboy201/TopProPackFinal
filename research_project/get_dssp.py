from Bio.PDB import *
import os
import sys

file = "C:\\alpha\\4quv.pdb"

os.chdir("C:\\alpha")
p = PDBParser()
structure = p.get_structure("4quv", file)
model = structure[0]
dssp = DSSP(model, file, dssp="mkdssp")


    

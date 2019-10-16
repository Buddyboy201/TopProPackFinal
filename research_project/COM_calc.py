import os
import sys
import numpy as np
from mendeleev import element

element_mass = {}
def get_element_mass(atom):
    atom.capitalize()
    if element_mass.get(atom) == None:
        element_mass[atom] = element(atom).atomic_weight
    return element_mass[atom]
    


def get_COM_res(residue, resid_coords, resid_name, file):
    COM = np.array([0.0, 0.0, 0.0])
    mass_sum = 0
    for i in range(len(resid_coords[residue])):
        try: COM[0] += resid_coords[residue][i][0]*get_element_mass(resid_name[residue][i].capitalize())
        except:
            print("atom: {} | resid: {} | bad file: {}".format(resid_name[residue][i].capitalize(), residue, file))
        try: COM[1] += resid_coords[residue][i][1]*get_element_mass(resid_name[residue][i].capitalize())
        except: continue
            #print("residue: " + resid_name[residue][i], "resid: "+residue, "bad file: "+file)
        try: COM[2] += resid_coords[residue][i][2]*get_element_mass(resid_name[residue][i].capitalize())
        except: continue
            #print("residue: " + resid_name[residue][i], "resid: "+residue, "bad file: "+file)
        try: mass_sum += get_element_mass(resid_name[residue][i].capitalize())
        except: continue
            #print("residue: " + resid_name[residue][i])
    COM[0] /= float(mass_sum)
    COM[1] /= float(mass_sum)
    COM[2] /= float(mass_sum)
    return COM

def get_COM_hashes(file):
    resid_coords = {}
    resid_name = {}
    with open(file) as pdb_file:
        for line in pdb_file:
            if line[0:4] == "ATOM" and line[12:16].strip(" ") != "CA" and line[12:16].strip(" ") != "O" and line[12:16].strip(" ") != "N" and line[12:16].strip(" ") != "C":
                #num = int(line[6:11].strip(" "))
                name = line[76:78].strip(" ")
                res_id = int(line[22:26].strip(" "))
                coordx = float(line[30:38].strip(" "))
                coordy = float(line[38:46].strip(" "))
                coordz = float(line[46:54].strip(" "))
                if resid_coords.get(res_id) == None:
                    resid_coords[res_id] = [(coordx, coordy, coordz)]
                    resid_name[res_id] = [name]
                else:
                    resid_coords[res_id].append((coordx, coordy, coordz))
                    resid_name[res_id].append(name)
    for i in resid_coords:
        resid_coords[i] = np.array(resid_coords[i])
        resid_name[i] = np.array(resid_name[i])
    return resid_coords, resid_name

def get_COM_all(file):
    coords_resid = {}
    resid_centroids = {}
    resid_coords, resid_name = get_COM_hashes(file)
    centroids = []
    #print(resid_coords, resid_name)
    count = 0
    for i in resid_coords:
        COM = get_COM_res(i, resid_coords, resid_name, file)
        centroids.append(list(COM))
        coords_resid[count] = i
        resid_centroids[i] = COM
        count += 1
    return np.array(centroids), coords_resid, resid_centroids

t = '''def main():
    file = "c:\\alpha\\4quv.pdb"
    centroids = get_COM_all(file)
    print(centroids)'''
    

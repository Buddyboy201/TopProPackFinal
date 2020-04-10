import centroid_protein
import atom
import residue
import energy
import os
import sys
import json
import time
import visualizer
import requests
import sqlalchemy

#TODO: gather 100 sample proteins, do distance distribution analysis, send report
#TODO: use same 100 proteins to generate benchmark runtime speeds for processing proteins manually v. loading from indv. json files/loading from bulk json file

class TPP_Engine:
    def __init__(self, exclude_backbone=False, distance_cutoff=6):
        os.chdir(r"C:\Users\aprak\PycharmProjects\TopProPack_v2_2")
        #self.STATIC_TOTAL_PAIRS_TABLE = None
        #self.STATIC_EPAIR_TABLE = None
        self.exclude_backbone = exclude_backbone
        self.distance_cutoff = distance_cutoff
        self.proteins = []
        self.E = energy.Energy()
        self.projects = {}

    def load_all_saved_data(self):
        files = os.listdir(os.getcwd()+"\\top_pro_pack_logs")
        #times = []
        #sizes = []
        #start_time = time.time()
        #last_time = start_time
        for name in files:
            #curr_time = time.time()
            self.load_protein(name)
            #diff = curr_time-last_time
            #times.append(diff)
            #last_time = curr_time
            #size = os.stat(os.getcwd()+"\\top_pro_pack_logs\\{}\\protein_data.json".format(name)).st_size
            #sizes.append(size)
            #print(name, curr_time-start_time, diff, size)
        #print("avg. rate: {} | avg. file size: {}".format(len(self.proteins)/(curr_time - start_time), sum(sizes)/float(len(self.proteins))))

    def engine_shutdown(self):
        for i in range(len(self.proteins)):
            if self.proteins[i].updated:
                start_time = time.time()
                self.proteins[i].save_protein_data()
                print("total saving time: {}".format(time.time()-start_time))

    def add_protein(self, project_name, name, file_path, json_load=True, data_load=True, data_url="https://files.rcsb.org/download/{}"):
        out = self.init_protein(name, file_path, json_load=json_load, data_load=data_load, data_url=data_url.format(name))
        if type(out) is Exception: print(out)
        else: self.projects[project_name].append(out)

    def add_dataset(self, project_name, proteins, modifers={"json_load": True, "data_load": True, "data_url": "https://files.rcsb.org/download/{}"}):
        prev_pdb = ""
        for pdb, file_path in proteins:
            if prev_pdb != pdb: self.add_protein(project_name, pdb, file_path, json_load=modifers["json_load"], data_load=modifers["data_load"], data_url=modifers["data_url"])
            prev_pdb = pdb

    def create_new_project(self, name = "project_{}", proteins=None, modifers={"json_load": True, "data_load": True, "data_url": "https://files.rcsb.org/download/{}"}):
        if name == "project_{}":
            name = name.format(len(self.projects)+1)
        print("Creating new project {}".format(name))
        self.projects[name] = []
        if proteins is not None:
            self.add_dataset(name, proteins)
        print("Project {} created!".format(name))

    def load_protein(self, name):
        file_path = os.getcwd() +"\\top_pro_pack_logs\\{}\\protein_data.json".format(name)
        P = centroid_protein.CentroidProtein("", "", load_json=True, json_data_file_path=file_path)
        #self.proteins.append(P)
        self.E.update_static_total_pairs_table(P.get_heatmap_data_centroid())
        return P

    def init_protein(self, name, file_path, json_load=True, data_load=True, data_url="https://files.rcsb.org/download/{}"):
        if name in os.listdir(os.getcwd() + "\\top_pro_pack_logs") and json_load:
            print("Attempting to load {} from JSON".format(name))
            self.load_protein(name)
        elif len(file_path) > 0:
            print("Atempting to process {} from directly from pdb file".format(name))
            try: P = centroid_protein.CentroidProtein(name, file_path, exclude_backbone=self.exclude_backbone)
            except:
                e = sys.exc_info()[0]
                return Exception(e)
            if len(P.residues) > 0:
                P.generate_centroid_cliques(distance_cutoff=self.distance_cutoff)
                # self.proteins.append(P)
                self.E.update_static_total_pairs_table(P.get_heatmap_data_centroid())
                return P
            else:
                return Exception("{} is empty".format(P.name))
        elif data_load and data_url is not None:
            print("Attempting to download/process {} from RCSB".format(name))
            try:
                P = centroid_protein.CentroidProtein(name, "", exclude_backbone=self.exclude_backbone,
                                                     download_data=data_load, data_url=data_url)
            except sqlalchemy.orm.exc.NoResultFound:
                return Exception("{} does not exist in RCSB database".format(name))
            if len(P.residues) > 0:
                P.generate_centroid_cliques(distance_cutoff=self.distance_cutoff)
                # self.proteins.append(P)
                self.E.update_static_total_pairs_table((P.get_heatmap_data_centroid()))
                return P
            else:
                return Exception("{} is empty".format(P.name))
        else:
            print("All processing attempts failed for {}, check provided info and try again".format(name))

    def initialize_static_total_pairs_table(self):
        pass

    def initialize_static_epair_table(self):
        pass



# distance distribution
# clique size distribution
# 2 dicts of cliques, tweak old bench cliques code

#os.chdir(r"C:\top_pro_pack_v2")
#file_path = r"C:\alpha\7prc.pdb"
#P = centroid_protein.CentroidProtein("7prc", file_path)
#P.generate_centroid_cliques()
#print(P.get_name())
#P.get_centroid_clique_distances()
#P.get_clique_frequency()
#print(P.get_json_dict())
#P.save_protein_data()
#file = open(r"C:\Users\aprak\OneDrive\Desktop\multiple_json_debug_log.txt.txt", "r")
#lines = file.readlines()
#lines = lines[0:len(lines)-2]
#lines = [line.split(" ")[0] for line in lines]
#init_files = ["C:\\alpha\\{}".format(f) for f in lines]
#E = TPP_Engine(distance_cutoff=100)
#for file in init_files:
#    try:
#        E.init_protein(os.path.split(file)[1], file, auto_load=False)
#    except:
#        print("an error occurred")
#distances = []
#for i in range(len(E.proteins)):
#    dist = E.proteins[i].get_centroid_clique_distances()
#    for d in dist:
#        distances.append(d)
#try: visualizer.draw_histogram(distances, "404 Dunbrack Proteins Distance Distribution", normalized=True)
#except: print("rip")
#E.engine_shutdown()

z = '''file_path = r"C:\alpha\1b12.pdb"
P = centroid_protein.CentroidProtein("1b12.pdb", file_path)
P.generate_centroid_cliques()
P.get_centroid_clique_distances()
P.get_clique_frequency()
print(P.name, len(P.residues), len(P.centroid_cliques))
P.save_protein_data()'''
#file = open(r"C:\Users\aprak\OneDrive\Desktop\multiple_json_debug_log.txt.txt", "r")
#lines = file.readlines()
#lines = lines[0:len(lines)-2]
#lines = [line.split(" ")[0] for line in lines]
#init_files = ["C:\\alpha\\{}".format(f) for f in lines]
#E = TPP_Engine()
#start_time = time.time()
#sizes = 0
#for file in init_files:
#    E.init_protein(os.path.split(file)[1], file)
    #sizes += os.stat(r"C:\Users\aprak\PycharmProjects\TopProPack_v2_2\top_pro_pack_logs\{}\protein_data.json".format(os.path.split(file)[1])).st_size
#curr_time = time.time()
#total_time = curr_time - start_time
#print(total_time, sizes/50.0/1000)
#print("processing time: {}".format(total_time))
#E.E.update_epair_values()
#print(E.E.STATIC_EPAIR_TABLE)
#E.engine_shutdown()
def test_func():
    os.chdir(r"C:\Users\aprak\PycharmProjects\TopProPack_v2_2")
    file_path = r"C:\Users\aprak\PycharmProjects\TopProPack_v2_2\{}"
    file = open(file_path.format("list_of_pdbs.txt"), "rt")
    pdb_list = [f[:len(f)-1] for f in file.readlines()[1:]]
    base_url = "https://files.rcsb.org/download/"
    E = TPP_Engine()
    start_time = time.time()
    bad_protein_count = 0
    prev_f = ""
    for f in pdb_list:
        name = f
        print("Attempting to process " + name)
        try:
            if f != prev_f: E.init_protein(name, "", json_load=False, data_url=base_url+f)
        except:
            print("Unknown error processing " + f)
            bad_protein_count += 1
        prev_f = f
    curr_time = time.time()
    print("Processing finished in {}".format(curr_time-start_time))
    print("# bad proteins: {}".format(int(bad_protein_count)))
    E.engine_shutdown()
    print("Saving finished in {}".format(time.time()-curr_time))

def func1():
    bugs = ["4v4m.pdb", "6trg.pdb", "6tt5.pdb", "6xtj.pdb", "6xy7.pdb"]

    base_url = "https://files.rcsb.org/download/"
    file = open(r"C:\alpha_2\big_dataset_debug_log_2.txt", "rt")
    bad_files = []
    for line in file:
        if "Unknown error" in line:
            bad_files.append(line.split(" ")[3].strip("\n"))
    print(bad_files)
    file.close()

#os.chdir(r"C:\Users\aprak\PycharmProjects\TopProPack_v2_2")
#file_path = r"C:\Users\aprak\PycharmProjects\TopProPack_v2_2\{}"
#file = open(file_path.format("list_of_pdbs.txt"), "rt")
#pdb_list = [f[:len(f)-1] for f in file.readlines()[1:]]
def func2():
    bad_files = []
    prev_f = ""
    E = TPP_Engine()
    for p in bad_files:
        name = p
        print("Attempting to process " + name)
        if p != prev_f: out = E.init_protein(name, "", json_load=False, data_url=base_url + p)
        if type(out) is Exception: print(out)
    E.E.update_epair_values()
    dists = []
    for i in E.proteins:
        for j in i.get_centroid_clique_distances():
            dists.append(j)
    print(dists)
    visualizer.draw_histogram(dists, "Centroid Clique Distances", normalized=True)

os.chdir(r"C:\Users\aprak\PycharmProjects\TopProPack_v2_2")
file = open(r"C:\alpha\pruned_data_list.txt", "rt")
pdb_list = [(f.strip("\n"), r"C:\alpha\{}".format(f.strip("\n"))) for f in file.readlines()]
pdb_list = [(f.strip("\n"), "") for f in file.readlines()]
file.close()
#print(pdb_list)
bad_files = []
E = TPP_Engine()
E.create_new_project()
E.add_dataset("project_1", pdb_list)



#os.chdir(r"C:\Users\aprak\PycharmProjects\TopProPack_v2_2")
#pdb_names = os.listdir(r"C:\Users\aprak\PycharmProjects\TopProPack_v2_2\top_pro_pack_logs")

#E = TPP_Engine()
#for pdb in pdb_names:
#    print(pdb)
#    E.init_protein(pdb, "")
#print(E.proteins)

#distances = []
#for i in E.proteins:
#    distances.append(i.centroid_clique_distances)
#visualizer.draw_histogram(distances, "Plot", normalized=True)
#E.E.update_epair_values()
#print(E.E.get_static_total_pairs_table())
#visualizer.draw_heatmap(E.E.STATIC_EPAIR_TABLE)











raw_pdb_time_efficiency_testing_procedure = '''start_time = time.time()
last_time = start_time
times = []
sizes = []
for file in init_files:
    try: E.init_protein(os.path.split(file)[1], file)
    except: print("an error occurred")
    curr_time = time.time()
    diff = curr_time - last_time
    last_time = curr_time
    size = os.stat(file).st_size
    times.append(diff)
    sizes.append(size)
    print(os.path.split(file)[1], curr_time - start_time, diff, size)
print("avg. rate: {} | avg. file size: {}".format(len(E.proteins)/(curr_time-start_time), sum(sizes)/float(len(E.proteins))))
E.engine_shutdown()'''

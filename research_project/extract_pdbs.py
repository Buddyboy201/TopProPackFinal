import os

os.chdir('C:\\alpha')
pdb_files = []
for f in os.listdir():
    if f[len(f)-4:] == ".pdb":
        pdb_files.append(f[0:len(f)-4]+"\n")
pdb_files[len(pdb_files)-1] = pdb_files[len(pdb_files)-1].strip("\n")
print(pdb_files)

file = open("pdb_list.txt", "w")
file.writelines(pdb_files)
file.close()


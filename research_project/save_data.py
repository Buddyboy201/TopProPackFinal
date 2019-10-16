import os
import requests
import csv

url = "http://dunbrack.fccc.edu/Guoli/users_html/cullpdb_pc40_res4.0_R0.3_inclNOTXRAY_d190627_entries462.22016"

myfile = requests.get(url)

os.chdir("c:\\alpha")
file = open("data.csv", "wb")
file.write(myfile.content)
file.close()

def getLineData(line):
    data = line.split()
    if data[2].find(".") != -1:
        data[2] = data[2][0:len(data[2])-1]
    for i in range(len(data)):
        if data[1] == "length":
            data[i] = data[i].lower()
        else:
            if i == 3:
                if data[i] == "NA":
                    data[i] = data[i].lower()
                else:
                    data[i] = float(data[i])
            elif i == 1:
                data[i] = int(data[i])
            elif i == 4 or i == 5:
                data[i] = float(data[i])
            else:
                data[i] = data[i].lower()
    return data
#print(len(file.readlines()))

data = open("data.csv", "r")
pfile = open("pruned_data.csv", "w")
plfile = open("pruned_data_list.txt", "w")
csv_writer = csv.writer(pfile, delimiter=",", quotechar='"')
for row in data:
    row_data = getLineData(row)
    if row_data[1] == "length":
        csv_writer.writerow(row_data)
    elif row_data[2] != "nmr" and row_data[3] != "na" and row_data[3] <= 4.0:
        csv_writer.writerow(row_data)
        plfile.write(row_data[0]+".pdb\n")
data.close()
pfile.close()
plfile.close()

import numpy as np
import os,sys

def masses(eleNames):

	eleMasses = []

	eles = {
			'Au': 196.97,
			'Ir': 192.22,
			'Pt': 195.8,
			'Pd': 106.42,
			'Rh': 102.91,
			'Ag': 107.81,
			'Ni': 58.7,
			'Fe': 55.85,
			'Zn': 65.39,
			'Hg': 200.59,
			'V' : 50.94,
			}

	eleMasses = eles[eleNames]

	return eleMasses

#* centerlize the center of mass
def CoM(incar):
	mass = []
	x = 0.0
	y = 0.0
	z = 0.0

	for list in incar:
		mass.append(masses(list[0]))
		x += masses(list[0])*list[1]
		y += masses(list[0])*list[2]
		z += masses(list[0])*list[3]


	x /= sum(mass)
	y /= sum(mass)
	z /= sum(mass)

	for i in range(len(incar)):
		incar[i][1] -= x
		incar[i][2] -= y
		incar[i][3] -= z

#* cell_size = boxadd + max of the cluster
def cell_size(incar, boxadd):

	size = []
	for i in range(len(incar)):
		size.append(incar[i][1])
		size.append(incar[i][2])
		size.append(incar[i][3])
	
	box = max(size) + boxadd

	return box

#* xyz2vasp
def xyz2vasp(file_name, eleNames, boxadd):
    #* read xyz file
    with open(f"{file_name}.xyz") as xyz:
        temp_incar = []
        lines = xyz.readlines()
        for i in range(2, len(lines)):
            temp_incar.append(lines[i].split())

    #* float xyz coords
    #* sort by element
    incar = []
    elenum = []
    #* count the number of each element
    for elename in eleNames:
        ele_counter = 0 
        for i in range(len(temp_incar)):
            if temp_incar[i][0] == elename:
                ele_counter += 1
                for j in range(1,4):
                    temp_incar[i][j] = float(temp_incar[i][j])

                incar.append(temp_incar[i])

        elenum.append(ele_counter)    

    CoM(incar)

    box = cell_size(incar, boxadd)


    #* write POSCAR
    with open(f"POSCAR_{file_name}", 'w') as poscar:
        poscar.write('POSCAR\n')
        poscar.write(str(box) + '\n')

        poscar.write("1.0 0.0 0.0\n")
        poscar.write("0.0 1.0 0.0\n")
        poscar.write("0.0 0.0 1.0\n")

        poscar.write(' '.join(eleNames) + '\n')
        poscar.write(' '.join(str(i) for i in elenum) + '\n')

        poscar.write('Direct\n')

        for elename in eleNames:
            for i in range(len(incar)):
                if incar[i][0] == elename:
                    x = str((incar[i][1] + box/2) / box)
                    y = str((incar[i][2] + box/2) / box)
                    z = str((incar[i][3] + box/2) / box)
                    poscar.write(x + '\t' + y + '\t' + z + '\n')

def vasp2xyz(file_name):

    clus = []
    coords = []
    box = None

    assert (file_name == 'POSCAR') or (file_name == 'CONTCAR'), \
        'file name must be POSCAR or OCNTCAR'

    with open(file_name) as vaspfile:
        lines = vaspfile.readlines()
        box = float(lines[1])               #* box size
        eleNames = lines[5].split()         #* element names
        elenum = lines[6].split()
        elenum = [int(i) for i in elenum]   #* element number
        #* get coords
        for i in range(8, 8 + sum(elenum)):
            temp = []
            temp = lines[i].split()
            temp = [float(i) for i in temp]
            coords.append(temp)
        
        #* transform into cartesian coords
        coords = np.array(coords)
        coords = coords * box - box/2.0

    #* write into clus = [element, x, y, z]
    counter = 0
    for i in range(len(eleNames)):
        for j in range(elenum[i]):
            clus.append([eleNames[i], coords[counter][0], coords[counter][1], coords[counter][2]])
            counter += 1
    
    #* centerlize the center of mass
    CoM(clus)
    
    #* write into xyz file
    with open(f"{file_name}.xyz", 'w') as xyz:
        xyz.write(f"{sum(elenum)}\n")
        xyz.write(f"{file_name}\n")
        for i in range(len(clus)):
            xyz.write(f"{clus[i][0]}\t{clus[i][1]}\t{clus[i][2]}\t{clus[i][3]}\n")
            
    







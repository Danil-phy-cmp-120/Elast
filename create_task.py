#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import subprocess
import os
import shutil

######### Ввести директорию с рассчетом a0 #########
path0 = raw_input('Enter the path to the calculation of equilibrium lattice parameter:\n')

######### Считывание CONTCAR #########

f = open(path0 + '/CONTCAR',"r")
contcar = f.readlines()
f.close()

# Считывание векторов трансляций #
primitive_vectors = np.zeros((3, 3))
n = 0
for i in [2,3,4]:
    inp = contcar[i].split()
    for j in range(3):
        primitive_vectors[n,j] = inp[j]
    n = n + 1
primitive_vectors = primitive_vectors * float(contcar[1])

# Число атомов #
num_atoms = sum([int(x) for x in contcar[6].split()])

# Считывание базиса #
basis = []
n = 0
for line in contcar[8:8+num_atoms]:
    inp = line.split()
    for i in inp:
        basis += [float(i)]
basis = np.array(basis)
basis.shape = ((len(basis)/3, 3))
#print basis


# Рассчет длины векторов трансляций #
len_pv = np.zeros(3)
for i in range(3):
        len_pv[i] = (sum(primitive_vectors[i,:]**2))**0.5


# Определение сингонии #
if ("{:.3f}".format(len_pv[0]) == "{:.3f}".format(len_pv[1]) and "{:.3f}".format(len_pv[1]) == "{:.3f}".format(len_pv[2])) or (float("{:.3f}".format(len_pv[0])) % float("{:.3f}".format(len_pv[1])) == 0 and float("{:.3f}".format(len_pv[0])) % float("{:.3f}".format(len_pv[2])) == 0):
    crystal_system = 'cube' 
    n_max = 3
elif ("{:.3f}".format(len_pv[0]) == "{:.3f}".format(len_pv[1]) and "{:.3f}".format(len_pv[1]) != "{:.3f}".format(len_pv[2])) or ("{:.3f}".format(len_pv[1]) == "{:.3f}".format(len_pv[2]) and "{:.3f}".format(len_pv[0]) != "{:.3f}".format(len_pv[1])) or ("{:.3f}".format(len_pv[2]) == "{:.3f}".format(len_pv[0]) and "{:.3f}".format(len_pv[2]) != "{:.3f}".format(len_pv[1])):
    crystal_system = 'tetr'
    n_max = 6
else:
    print 'Calculation is possible only for cubic and tetragonal crystal system'


######### Считывание OUTCAR #########

f = open(path0 + '/OUTCAR',"r")
outcar = f.readlines()
f.close()

flag = False
magmom = np.zeros(basis.shape[0])
n = 0
for line in outcar[len(outcar)-200:len(outcar)]:
    inp = line.split()

    if len(inp) > 3 and inp[0] == 'external' and inp[1] == 'pressure':
        print 'external pressure =',  inp[3]
        if abs(float(inp[3])) > 5.0:
            print('The pressure is high, restart the ionic relaxation with a small EDIFG')

    if len(inp) > 3 and flag == True and inp[0].isdigit() == True and n < basis.shape[0]:
        magmom[n] = float(inp[len(inp)-1])
        n = n + 1
    if len(inp) > 1 and inp[0] == 'magnetization':
        flag = True


######### Изменение некоторых тегов в INCAR #########

f = open(path0 + '/INCAR',"r")
incar = f.readlines()
f.close()

for i in range(len(incar)):
    inp = incar[i].split()
    if len(inp) > 2 and inp[0] == 'NSW':
        incar[i] = 'NSW = 0\n'
    if len(inp) > 2 and inp[0] == 'MAGMOM':
        incar[i] = 'MAGMOM ='
        for j in range(magmom.size):
            incar[i] = incar[i] + ' {}'.format(magmom[j])
        incar[i] = incar[i] + '\n'

######### Создание искаженных POSCAR и запись INCAR, KPOINTS, POTCAR #########
if os.path.exists('bulk') == False:
    os.mkdir('bulk')

DELTAS = np.linspace(-0.03, 0.03, 7) # Искажения от -3 до 3 %

for n in range(n_max):

    if os.path.exists('bulk/D{}'.format(n)) == False:
        os.mkdir('bulk/D{}'.format(n))

    for delta in DELTAS:

        if os.path.exists('bulk/D{}/{}'.format(n, delta)) == False:
            os.mkdir('bulk/D{}/{}'.format(n, delta))

        if crystal_system == 'cube':
            T = np.array([[1 + delta, 0, 0],
                         [0, 1 + delta, 0],
                         [0, 0, 1 + delta],
   
                         [1 + delta, 0, 0],
                         [0, 1 - delta, 0],
                         [0, 0, 1/(1-delta**2)],
             
                         [1, delta, 0],
                         [delta, 1, 0],
                         [0, 0, 1/(1-delta**2)]])

        elif crystal_system == 'tetr':
            T = np.array([[1/(1-delta**2), 0, delta],     
                         [0, 1, 0],
                         [delta, 0, 1],
   
                         [1, delta, 0],                  
                         [delta, 1, 0],
                         [0, 0, 1/(1-delta**2)],
             
                         [1 + delta, 0, 0],              
                         [0, 1 - delta, 0],
                         [0, 0, 1/(1-delta**2)],

                         [1 + delta, 0, 0],              
                         [0, 1/(1-delta**2), 0],
                         [0, 0, 1-delta],
  
                         [1 + delta, 0, 0],              
                         [0, 1 + delta, 0],
                         [0, 0, 1 + delta],

                         [1, 0, 0],              
                         [0, 1, 0],
                         [0, 0, 1 + delta]])

        D =  T[n*3:n*3+3,:]
        primitive_vectors_new = np.dot(primitive_vectors, D)
        
        poscar = contcar
        poscar[1] = '   1.0\n'
        for i in range(3):
            poscar[2+i] = '     '
            for j in range(3):
                poscar[2+i] = poscar[2+i] + str(primitive_vectors_new[i, j]) + '    '
            poscar[2+i] = poscar[2+i] + '\n'

        f = open('bulk/D{}/{}/POSCAR'.format(n, delta), "wb") # Запись POSCAR
        for i in range(len(poscar)):
            f.write(poscar[i])
        f.close

        f = open('bulk/D{}/{}/INCAR'.format(n, delta), "wb") # Запись INCAR
        for i in range(len(incar)):
            f.write(incar[i])
        f.close

        shutil.copyfile(path0 + '/POTCAR', 'bulk/D{}/{}/POTCAR'.format(n, delta)) # Запись POTCAR

        shutil.copyfile(path0 + '/KPOINTS', 'bulk/D{}/{}/KPOINTS'.format(n, delta)) # Запись KPOINTS

    if n == 0:
        DELTAS = np.delete(DELTAS, 3)

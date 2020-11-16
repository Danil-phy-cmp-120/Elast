#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import os
import sys
from scipy import constants


def get_weight():
    MASS = {'H':1.0079, 'He':4.0026, 'Li':6.941, 'Be':9.0122, 'B':10.811, 'C':12.011, 'N':14.007, 'O':15.999, 'F':18.998, 'Ne':20.18, 'Na':22.99, 'Mg':24.305, 'Al':26.982, 'Si':28.086, 'P':30.974, 'S':32.066, 'Cl':35.453, 'Ar':39.948, 'K':39.099, 'Ca':40.078, 'Sc':44.956, 'Ti':47.88, 'V':50.942, 'Cr':51.996, 'Mn':54.938, 'Fe':55.847, 'Co':58.933, 'Ni':58.693, 'Cu':63.546, 'Zn':65.39, 'Ga':69.723, 'Ge':72.61, 'As':74.922, 'Se':78.96, 'Br':79.904, 'Kr':83.8, 'Rb':85.468, 'Sr':87.62, 'Y':88.906, 'Zr':91.224, 'Nb':92.906, 'Mo':95.94, 'Tc':98.906, 'Ru':101.07, 'Rh':102.91, 'Pd':106.42, 'Ag':107.87, 'Cd':112.41, 'In':114.82, 'Sn':118.71, 'Sb':121.76, 'Te':127.60, 'I':126.9, 'Xe':131.29, 'Cs':132.91, 'Ba':137.33, 'La':138.91, 'Ce':140.12, 'Pr':140.91, 'Nd':144.24, 'Pm':144.91, 'Sm':150.36, 'Eu':151.97, 'Gd':157.25, 'Tb':158.93, 'Dy':162.5, 'Ho':164.93, 'Er':167.26, 'Tm':168.93, 'Yb':173.04, 'Lu':174.97, 'Hf':178.49, 'Ta':180.95, 'W':183.84, 'Re':186.21, 'Os':190.23, 'Ir':192.22, 'Pt':195.08, 'Au':196.97, 'Hg':200.59, 'Tl':204.38, 'Pb':207.2, 'Bi':208.98, 'Po':209.98, 'At':209.99, 'Rn':222, 'Fr':223.02, 'Ra':226.03, 'Ac':227.03, 'Th':232.04, 'Pa':231.04, 'U':238.03, 'Np':237.05, 'Pu':244.06, 'Am':243, 'Cm':247.07, 'Bk':247, 'Cf':251, 'Es':252, 'Fm':257.1, 'Md':258, 'No':259, 'Lr':260.11, 'Rf':261.11, 'Db':262.11, 'Sg':263.12, 'Bh':262.12, 'Hs':265, 'Mt':268, 'Ds':281, 'Rg':280}

    f = open('bulk/D0/0.0/CONTCAR', "r")
    contcar = f.readlines()
    f.close()

    element = []
    inp1 = contcar[5].split()
    inp2 = contcar[6].split()

    for i in range(len(inp1)):
        element += [inp1[i]]*int(inp2[i])

    mass = 0.0
    for i in range(len(element)):
        mass = mass + MASS[element[i]]

    return mass/1000.0, sum([int(x) for x in contcar[6].split()])


def get_energy(p):
    f = open(p, "r")
    outcar = f.readlines()
    f.close()

    for line in outcar[len(outcar)-200:len(outcar)]:
        inp = line.split()
        if len(inp) > 4 and inp[4] == 'energy(sigma->0)':
            energy = float(inp[6])

    return energy

def get_volume(p):
    f = open(p, "r")
    outcar = f.readlines()
    f.close()

    for line in outcar[len(outcar)-200:len(outcar)]:
        inp = line.split()
        if len(inp) > 3 and inp[0] == 'volume' and inp[1] == 'of' and inp[2] == 'cell':
            volume = float(inp[4])

    return volume


V0 = get_volume('bulk/D0/0.0/OUTCAR') * 10**-30

path = os.listdir('bulk')
DELTAS = np.linspace(-0.03, 0.03, 7)

# Рассчет отдельных модулей упругости для кубической симметрии #
if len(path) == 3:
    energy = np.zeros((DELTAS.size, 3))
    alfa = np.zeros(3)

    for n in range(3):
       for i in range(DELTAS.size):
          if n == 0 and DELTAS[i] == 0.0:
              energy[i, :] = [get_energy('bulk/D{}/{}/OUTCAR'.format(n, DELTAS[i]))]*3
          elif DELTAS[i] != 0.0:
              energy[i, n] = get_energy('bulk/D{}/{}/OUTCAR'.format(n, DELTAS[i]))
 
    np.savetxt('energy.dat', energy, fmt='%.6f')
    energy = 1.6*10**-28 * (energy - energy[3,0])/V0

    for n in range(3):
        alfa[n] = np.polyfit(DELTAS, energy[:,n], 2)[0]

    f = open('out.dat', "wb")
    np.savetxt(f, np.column_stack((DELTAS, energy)), fmt='%.5f')
    f.write('\n')

    np.savetxt(f, alfa, fmt='%.5f')
    f.write('\n')

    B = 2*alfa[0]/9
    C = alfa[1]/2
    C11 = (3*B + 4*C)/3
    C12 = (3*B - 2*C)/3
    C44 = alfa[2]/2

    f.write('B = '+ str(B) + '       ')
    f.write('C = ' + str(C) + '\n')
    f.write('\n')
    f.write('Cij:' + '\n')
    f.write('C11 = ' + str(C11) + '\n')
    f.write('C12 = ' + str(C12) + '\n')
    f.write('C44 = ' + str(C44) + '\n')
    f.write('\n')

    BV=(C11+2*C12)/3
    BR=(C11+2*C12)/3
    GV=(C11-C12+3*C44)/5
    GR=(5*C44*(C11-C12))/(4*C44+3*(C11-C12))
    BH=(BV+BR)/2
    GH=(GV+GR)/2

    f.write('BV = ' + str(BV) + '       ')
    f.write('GV = ' + str(GV) + '\n')
    f.write('BR = ' + str(BR) + '       ')
    f.write('GR = ' + str(GR) + '\n')
    f.write('BH = ' + str(BH) + '       ')
    f.write('GH = ' + str(GH) + '\n')
    f.write('\n')
 
    mu = (3*BH-2*GH)/(6*BH+2*GH)
    E = (9*BH*GH)/(3*BH+GH)
    func = (3*(2*(2*(1+mu)/(3-6*mu))**1.5 + ((1+mu)/(3-3*mu))**1.5)**-1)**(1.0/3.0)
    Tetta = constants.hbar * ( 6*3.1415**2 * V0**0.5 * get_weight()[1] )**(1.0/3.0) * func * ((BH * 10**9 * constants.N_A) / ( constants.k**2 * get_weight()[0] ))**0.5

    f.write('mu = ' + str(mu) + '\n')
    f.write('E = ' + str(E)  + '\n')
    f.write('Tetta = ' + str(Tetta) + '\n')
    f.write('\n')

    if (C44 > 0 and C11 > abs(C12) and C11+2*C12 > 0):
        f.write('The crystal lattice is stable under distortion')
    else:
        f.write('The crystal lattice is unstable under distortion')

    f.close()

# Рассчет отдельных модулей упругости для тетрагональной симметрии #
if len(path) == 6:
    energy = np.zeros((DELTAS.size, 6))
    alfa = np.zeros(6)

    for n in range(6):
       for i in range(DELTAS.size):
          if n == 0 and DELTAS[i] == 0.0:
              energy[i, :] = [get_energy('bulk/D{}/{}/OUTCAR'.format(n, DELTAS[i]))]*6
          elif DELTAS[i] != 0.0:
              energy[i, n] = get_energy('bulk/D{}/{}/OUTCAR'.format(n, DELTAS[i]))
 
    energy = 1.6*10**-28 * (energy - energy[3,0])/V0
    #np.savetxt('energy.dat', energy, fmt='%.6f')

    for n in range(6):
        alfa[n] = np.polyfit(DELTAS, energy[:,n], 2)[0]

    f = open('out.dat', "wb")
    np.savetxt(f, np.column_stack((DELTAS, energy)), fmt='%.5f')
    f.write('\n')

    np.savetxt(f, alfa, fmt='%.5f')
    f.write('\n')

    C11 = (alfa[2] + 2*alfa[3] - 3*alfa[5] +alfa[4])/3
    C12 = (2*alfa[3] - 2*alfa[2] - 3*alfa[5] + alfa[4])/3
    C13 = (alfa[4] - 4*alfa[3] + 3*alfa[5] + alfa[2])/6
    C33 = 2*alfa[5]
    C44 = alfa[0]/2
    C66 = alfa[1]/2

    B= (C33*(C11+C12)-2*C13**2)/(C11+C12+2*C33-4*C13)

    f.write('B = '+ str(B) + '       ')
    f.write('\n')
    f.write('Cij' + '\n')
    f.write('C11 = ' + str(C11) + '\n')
    f.write('C12 = ' + str(C12) + '\n')
    f.write('C13 = ' + str(C13) + '\n')
    f.write('C33 = ' + str(C33) + '\n')
    f.write('C44 = ' + str(C44) + '\n')
    f.write('C66 = ' + str(C66) + '\n')
    f.write('\n')

    BV=(2*(C11+C12)+4*C13+C33)/9
    BR=(C33*(C11+C12)-2*C13**2)/(C11+C12+2*C33-4*C13)
    GV=(12*C44+12*C66+C11+C12+2*C33-4*C13)/30
    GR=(5*C44*C66*(C33*(C11+C12)-2*C13**2))/(2*(C44+C66)*(C33*(C11+C12)-2*C13**2)+3*BV*C44*C66)
    BH=(BV+BR)/2
    GH=(GV+GR)/2

    f.write('BV = ' + str(BV) + '       ')
    f.write('GV = ' + str(GV) + '\n')
    f.write('BR = ' + str(BR) + '       ')
    f.write('GR = ' + str(GR) + '\n')
    f.write('BH = ' + str(BH) + '       ')
    f.write('GH = ' + str(GH) + '\n')
    f.write('\n')

    mu = (3*BH-2*GH)/(6*BH+2*GH)
    E = (9*BH*GH)/(3*BH+GH)
    func = (3*(2*(2*(1+mu)/(3-6*mu))**1.5 + ((1+mu)/(3-3*mu))**1.5)**-1)**(1.0/3.0)
    Tetta = constants.hbar * ( 6*3.1415**2 * V0**0.5 * get_weight()[1] )**(1.0/3.0) * func * ((BH * 10**9 * constants.N_A) / ( constants.k**2 * get_weight()[0] ))**0.5

    f.write('mu = ' + str(mu) + '\n')
    f.write('E = ' + str(E)  + '\n')
    f.write('Tetta = ' + str(Tetta) + '\n')
    f.write('\n')

    if (C11 > 0 and C33 > 0 and C44 > 0 and C66 > 0 and C11 - C12 > 0 and C11 - 2*C13 + C33 > 0 and 2*C11 + 2*C12 + 4*C13 + C33 > 0):
        f.write('The crystal lattice is stable under distortion')
    else:
        f.write('The crystal lattice is unstable under distortion')

    f.close()

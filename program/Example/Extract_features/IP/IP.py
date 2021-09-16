# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 15:33:51 2021

@author: Tong Zhou
"""

import pandas as pd
import numpy as np
import sys

# input PDB file
file = sys.argv[1].split('/')[-1]
pdbID = file.split('.')[0]
chain = pdbID.split('_')[1]
pdb = open(sys.argv[1], 'r')
content = pdb.read()
pdb.close()
line = content.strip().split("\n")

# table lookup
ATOM = []
IP_i = []
for i in range(0,len(line)):
    if line[i][0:4] == 'ATOM' and line[i][21] == chain and ((line[i][13] == 'C' and line[i][14] == 'A')):
        ATOM_split = line[i].split()
        ATOM.append(ATOM_split)
        pdbname = line[i][17:20]
        if pdbname =='GLY':
            resname = 'G'
            IP = [0.993]
        elif pdbname =='ALA':
            resname = 'A'
            IP = [0.736]
        elif pdbname == 'VAL':
            resname = 'V'
            IP = [0.853]
        elif pdbname == 'LEU':
            resname = 'L'
            IP = [0.567]
        elif pdbname == 'ILE':
            resname = 'I';
            IP = [0.906]
        elif pdbname == 'PRO':  
            resname = 'P'
            IP = [0.655]
        elif pdbname == 'PHE':
            resname = 'F'
            IP = [1.478]
        elif pdbname == 'TRP':
            resname = 'W';
            IP = [1.169]
        elif pdbname == 'TYR':
            resname = 'Y'
            IP = [1.723]
        elif pdbname == 'SER':
            resname = 'S'
            IP = [1.058]
        elif pdbname == 'THR':
            resname = 'T'
            IP = [1.032]
        elif pdbname == 'CYS':
            resname = 'C'
            IP = [1.668]
        elif pdbname == 'MET':
            resname = 'M'
            IP = [1.118]
        elif pdbname == 'ASN':
            resname = 'N'
            IP = [1.282]
        elif pdbname == 'GLN':
            resname = 'Q'
            IP = [1.185]
        elif pdbname == 'ASP':
            resname = 'D'
            IP = [0.523]
        elif pdbname == 'GLU':
            resname = 'E'
            IP = [0.351]
        elif pdbname == 'HIS':
            resname = 'H'
            IP = [1.584]
        elif pdbname == 'LYS':
            resname = 'K'
            IP = [1.393]
        elif pdbname == 'ARG':
            resname = 'R'
            IP = [2.628]
        IP_i.append(IP)
        
IP_i = np.array(IP_i)
pos = [[x[5]] for x in ATOM]
pos = np.array(pos)
data = pd.DataFrame(np.hstack((pos,IP_i)))
data.columns = ['pos','IP']
data.to_csv(pdbID+'_IP.csv',index=False)
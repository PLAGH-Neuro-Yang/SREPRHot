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
F_i = []
for i in range(0,len(line)):
    if line[i][0:4] == 'ATOM' and line[i][21] == chain and ((line[i][13] == 'C' and line[i][14] == 'A')):
        ATOM_split = line[i].split()
        ATOM.append(ATOM_split)
        pdbname = line[i][17:20]
        if pdbname =='GLY':
            resname = 'G'
            F = [0,2,0.005]
        elif pdbname =='ALA':
            resname = 'A'
            F = [0,2,0.0373]
        elif pdbname =='VAL':
            resname = 'V'
            F = [0,2,0.0057]
        elif pdbname == 'LEU':
            resname = 'L'
            F = [0,2,0]
        elif pdbname == 'ILE':
            resname = 'I';
            F = [0,2,0]
        elif pdbname == 'PRO':  
            resname = 'P'
            F = [0,2,0.0198]
        elif pdbname == 'PHE':
            resname = 'F'
            F = [0,2,0.0946]
        elif pdbname == 'TRP':
            resname = 'W';
            F = [0,3,0.0548]
        elif pdbname == 'TYR':
            resname = 'Y'
            F = [0,3,0.0516]
        elif pdbname == 'SER':
            resname = 'S'
            F = [0,4,0.0829]
        elif pdbname == 'THR':
            resname = 'T'
            F = [0,4,0.0941]
        elif pdbname == 'CYS':
            resname = 'C'
            F = [0,2,0.0829]
        elif pdbname == 'MET':
            resname = 'M'
            F = [0,2,0.0823]
        elif pdbname == 'ASN':
            resname = 'N'
            F = [0,4,0.0036]
        elif pdbname == 'GLN':
            resname = 'Q'
            F = [0,4,0.0761]
        elif pdbname == 'ASP':
            resname = 'D'
            F = [-1,4,0.1263]
        elif pdbname == 'GLU':
            resname = 'E'
            F = [-1,4,0.0058]
        elif pdbname == 'HIS':
            resname = 'H'
            F = [1,4,0.0242]
        elif pdbname == 'LYS':
            resname = 'K'
            F = [1,2,0.0371]
        elif pdbname == 'ARG':
            resname = 'R'
            F = [1,4,0.0959]
        F_i.append(F)
F_i = np.array(F_i)
pos = [[x[5]] for x in ATOM]
pos = np.array(pos)
data = pd.DataFrame(np.hstack((pos,F_i)))
data.columns = ['pos','Nec','Nphb','Eiip']
data.to_csv(pdbID+'_PC.csv',index=False)
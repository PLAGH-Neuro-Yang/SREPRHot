# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 15:33:51 2021

@author: Tong Zhou
"""

import pandas as pd
import sys
import os

# input PDB file
outdir = "./Output/"
if not os.path.exists(outdir):
    os.makedirs(outdir)
file = sys.argv[1].split('/')[-1]
pdbID = file.split('.')[0]
chain = pdbID.split('_')[1]
pdb = open(sys.argv[1], 'r')
content = pdb.read()
pdb.close
line = content.strip().split("\n")

# confirm number of 1st residue in FASTA file
SEQADV = []
for i in range(0,len(line)):
    if line[i][0:6] == 'SEQADV':
        SEQADV.append(line[i])
if SEQADV == []:
    Seq1stNum = 1
else:
    Seq1stNum = eval(SEQADV[0].split()[4])

# Physicochemical characteristics
PC = pd.read_csv("./Input/"+pdbID+"_PC.csv",
                 encoding='utf-8', sep=',',quoting=1)

# PSSM
PSSM = pd.read_table("./Input/"+pdbID+".asn_matrix.txt",
                     encoding='utf-8', sep='\t+',quoting=1,engine='python')
PSSM['P']=range(Seq1stNum,len(PSSM)+Seq1stNum)
PSSM = PSSM[['P','L','P.1']]
PSSM.columns = ['pos','PL','PP']

# CX/DPX
tbl_bound = pd.read_table("./Input/"+pdbID+"_bound.tbl",
                          encoding='utf-8', sep=' +',quoting=1,
                          engine='python',skiprows=9,header=None)
tbl_bound.columns = ['chain','pos','AA',
                     'BmDPXa','BsdDPXa','BmDPXs','BsdDPXs','BmaxDPX','BminDPX',
                     'BmCXa','BsdCXa','BmCXs','BsdCXs','BmaxCX','BminCX']
tbl_bound = tbl_bound[['chain','pos','AA','BsdDPXa','BmDPXs','BmCXs']]
tbl_bound = tbl_bound.loc[tbl_bound['chain']==chain]
tbl_unbound = pd.read_table("./Input/"+pdbID+"_unbound.tbl",
                            encoding='utf-8', sep=' +',quoting=1,
                            engine='python',skiprows=9,header=None)
tbl_unbound.columns = ['chain','pos','AA',
                       'UmDPXa','UsdDPXa','UmDPXs','UsdDPXs','UmaxDPX','UminDPX',
                       'UmCXa','UsdCXa','UmCXs','UsdCXs','UmaxCX','UminCX']
tbl_unbound = tbl_unbound.loc[tbl_unbound['chain']==chain]
tbl_unbound = tbl_unbound[['pos','UsdDPXs','UmDPXs']]
tbl = pd.merge(tbl_bound, tbl_unbound, on=['pos'])
tbl['DmDPXs'] = tbl['BmDPXs']-tbl['UmDPXs']
tbl = tbl[['pos','BsdDPXa','DmDPXs','BmCXs','UsdDPXs']]

# ASA
ASA_bound = open("./Input/"+pdbID+"_bound.rsa", 'r')
content = ASA_bound.read()
ASA_bound.close()
line1 = content.strip().split("\n")
ASA_bound = []
for i in range(0,len(line1)):
    if line1[i][0:3] == 'RES' and line1[i][8] == chain:
        ASA_bound_split = line1[i].split()
        ASA_bound.append(ASA_bound_split)
ASA_bound = [[x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13]] for x in ASA_bound]
ASA_bound = pd.DataFrame(ASA_bound)
ASA_bound.columns = ['pos','BaASAaa','BrASAaa','BaASAts','BrASAts','BaASAmc',
                     'BrASAmc','BaASAnp','BrASAnp','BaASAap','BrASAap']
ASA_unbound = open("./Input/"+pdbID + '_unbound.rsa', 'r')
content = ASA_unbound.read()
ASA_unbound.close()
line2 = content.strip().split("\n")
ASA_unbound = []
for i in range(0,len(line2)):
    if line2[i][0:3] == 'RES' and line2[i][8] == chain:
        ASA_unbound_split = line2[i].split()
        ASA_unbound.append(ASA_unbound_split)
ASA_unbound = [[x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13]] for x in ASA_unbound]
ASA_unbound = pd.DataFrame(ASA_unbound)
ASA_unbound.columns = ['pos','UaASAaa','UrASAaa','UaASAts','UrASAts','UaASAmc',
                       'UrASAmc','UaASAnp','UrASAnp','UaASAap','UrASAap']
ASA_bound = ASA_bound[['pos','BaASAnp','BrASAnp']]
ASA_unbound = ASA_unbound[['pos','UrASAts']]
ASA = pd.merge(ASA_bound, ASA_unbound, on=['pos'])
ASA['pos'] = ASA['pos'].apply(int)

# Solvent exposure
HSE = pd.read_table("./Input/"+pdbID+"_HSE.txt",
                    encoding='utf-8', quoting=1,sep=' +',engine='python')
HSE['Position']=range(Seq1stNum,len(HSE)+Seq1stNum)
HSE.columns = ['pos','Residue','HSE_up','HSE_down','CN']

# IP
IP = pd.read_csv("./Input/"+pdbID+"_IP.csv",
                 encoding='utf-8', sep=',',quoting=1)

# ANN
ANN_H = pd.read_csv("./Input/"+pdbID+"_ANN_Hydrophobicity.txt",
                    encoding='utf-8', sep=' ',quoting=1)
ANN_H = ANN_H[['Resid','Kw']][ANN_H['chain']==chain]
ANN_H.columns = ['pos','Kh']
ANN_H = ANN_H.reset_index(drop=True)
ANN_S = pd.read_csv("./Input/"+pdbID+"_ANN_SAS.txt",
                    encoding='utf-8', sep=' ',quoting=1)
ANN_S = ANN_S[['Resid','Kw']][ANN_S['chain']==chain]
ANN_S.columns = ['pos','Ks']
ANN_S = ANN_S.reset_index(drop=True)
ANN = pd.merge(ANN_H, ANN_S, on=['pos'])

# Integration
M1 = pd.merge(PC, tbl, on=['pos'])
M2 = pd.merge(M1,ASA, on=['pos'])
M3 = pd.merge(M2,IP, on=['pos'])
M4 = pd.merge(M3,ANN, on=['pos'])
M5 = pd.merge(M4,PSSM, on=['pos'])
M6 = pd.merge(M5,HSE, on=['pos'])
res = M6[['pos','Residue','DmDPXs','BaASAnp','IP','UrASAts','BsdDPXa','Eiip',
          'BrASAnp','Ks','HSE_down','PL','PP','UsdDPXs','Nphb','CN','Nec',
          'BmCXs','HSE_up','Kh']]
res.to_csv(outdir+pdbID+"_feature.csv",index=False)


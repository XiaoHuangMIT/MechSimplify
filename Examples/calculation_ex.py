import os
import sys
import numpy as np
import pandas as pd
from molSimplify.Classes.mol3D import mol3D
from molSimplify.Classes.atom3D import atom3D



df = pd.read_csv('opt_natoms.csv')
#Spin-splitting filtering
df = df[df.Ehsls > 0]
df = df[df.Ehsls < 30]
#Size: 0 - 74 
df = df[df.natoms < 75]



idxs = df.index.values.tolist()

for idx in idxs:
    
    row = df.loc[idx]
    refcode = row['refcode']
    metal = row['metal']
    charge = row['ox_csd']
    appair = None
    if row['apnum'] > 1:
        appair = [row['ap1'],row['ap2']]
    
    is_co2 = False #co2 doesn`t have IS state
    if metal == 'Co' and charge == 2:
        is_co2 = True
        
    pp1,pp2 = row['pp1'],row['pp2']
    pp1 += 1 #Since molSimplify starts with 0
    pp2 += 1
    mol2ls = row['mol2ls']
    distls = get_atoms_distance(mol2ls,pp1,pp2)
    mol2hs = row['mol2hs']
    disths = get_atoms_distance(mol2hs,pp1,pp2)
    if not is_co2:
        mol2is = row['mol2is']
        distis = get_atoms_distance(mol2is,pp1,pp2)
    
    if metal == 'Fe' and charge == 2:
        prep_tera_input(refcode,'LS_C',1,2,mol2ls,ap_pair=appair,COGEF=[pp1,pp2,distls,10,50])
        prep_tera_input(refcode,'IS_C',3,2,mol2is,ap_pair=appair,COGEF=[pp1,pp2,distis,10,50])
        prep_tera_input(refcode,'HS_C',5,2,mol2hs,ap_pair=appair,COGEF=[pp1,pp2,disths,10,50])
    
    elif metal == 'Fe' and charge == 3:
        prep_tera_input(refcode,'LS_C',2,3,mol2ls,ap_pair=appair,COGEF=[pp1,pp2,distls,10,50])
        prep_tera_input(refcode,'IS_C',4,3,mol2is,ap_pair=appair,COGEF=[pp1,pp2,distis,10,50])
        prep_tera_input(refcode,'HS_C',6,3,mol2hs,ap_pair=appair,COGEF=[pp1,pp2,disths,10,50])
    
    elif metal == 'Co' and charge == 2:
        prep_tera_input(refcode,'LS_C',2,2,mol2ls,ap_pair=appair,COGEF=[pp1,pp2,distls,10,50])
        prep_tera_input(refcode,'HS_C',4,2,mol2hs,ap_pair=appair,COGEF=[pp1,pp2,disths,10,50])
    
    elif metal == 'Co' and charge == 3:
        prep_tera_input(refcode,'LS_C',1,3,mol2ls,ap_pair=appair,COGEF=[pp1,pp2,distls,10,50])
        prep_tera_input(refcode,'IS_C',3,3,mol2is,ap_pair=appair,COGEF=[pp1,pp2,distis,10,50])
        prep_tera_input(refcode,'HS_C',5,3,mol2hs,ap_pair=appair,COGEF=[pp1,pp2,disths,10,50])
    
    elif metal == 'Mn' and charge == 2:
        prep_tera_input(refcode,'LS_C',2,2,mol2ls,ap_pair=appair,COGEF=[pp1,pp2,distls,10,50])
        prep_tera_input(refcode,'IS_C',4,2,mol2is,ap_pair=appair,COGEF=[pp1,pp2,distis,10,50])
        prep_tera_input(refcode,'HS_C',6,2,mol2hs,ap_pair=appair,COGEF=[pp1,pp2,disths,10,50])

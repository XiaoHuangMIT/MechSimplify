#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
import numpy as np
import pandas as pd
from molSimplify.Classes.mol3D import mol3D


# In[2]:


def prep_orca_input(refcode, charge, spin, spinval, mol2, multiPP = None, dist_constraint = None, EFEI = None):


    #Write orca EFEI geometry optimization input, including distance-constraint and force_modified variations, that
    #can be run on Expanse/Slurm and used with job_manager
    #Method by default:
    #B3LYP*(VMN5)-D3BJ: 0.15 a, 0.77b, 0.81c
    #Basis: def2-svp
    #Solvent: gaseous phase

    #Inputs:
    #refcode: 6-letter refcode of complex
    #charge: charge of metal and molecule (neutral ligands only)
    #spin: LS/IS/HS
    #spinval: 1/3/5 (ex: Fe2+)
    #mol2: mol2 string of molecule to be written into xyz

    #multiPP: if the molecule has more than one pair of pulling points
    #format: None/[pp1idx,pp2idx]

    #dist_constraint: if optimization with distance constraint between two atoms (pps) should be performed
    #format: None/[pp1idx,pp2idx,distance(A)]

    #EFEI: if optimization under force should be performed
    #format: None/[pp1idx,pp2idx,force(nN)]


    #Obtaining base name
    basename = refcode + '_' + spin
    if multiPP != None:
        basename = basename + '_' + str(multiPP[0]) + '_' + str(multiPP[1])
    if dist_constraint != None:
        basename = basename + '_' + str(dist_constraint[2])
    elif EFEI != None:
        basename = basename + '_' + str(EFEI[2])

    #Make folder
    if not os.path.exists(basename + '/'):
        os.mkdir(basename + '/')

    #Generating jobscript:
    with open(basename + '/' + basename + '_jobscript', 'w') as f:
        f.write('#!/usr/bin/env bash\n')
        f.write('#SBATCH --job-name=' + basename + '\n')
        f.write('#SBATCH --output=' + basename + '.out' + '\n')
        f.write('#SBATCH --mem=48GB\n')
        f.write('#SBATCH --time=48:00:00\n')
        f.write('#SBATCH --nodes=1\n')
        f.write('#SBATCH --ntasks-per-node=8\n')
        f.write('#SBATCH -p shared\n')
        f.write('#SBATCH -A TG-CHE140073\n')
        f.write('#SBATCH --export=ALL\n')
        f.write('\n')
        f.write('module load cpu/0.15.4  gcc/9.2.0  openmpi/4.1.1 orca/5.0.1\n')
        f.write('\n')
        f.write('# Change to working directory\n')
        f.write('export SCRDIR=/scratch/$USER/job_$SLURM_JOB_ID\n')
        #f.write('cp $SLURM_SUBMIT_DIR/' + basename + '.in ' + '$SCRDIR\n') check again
        f.write('cp $SLURM_SUBMIT_DIR/* $SCRDIR\n')
        f.write('mkdir $SLURM_SUBMIT_DIR/scr\n')
        f.write('cd $SCRDIR\n')
        f.write('\n')
        f.write('ulimit -s unlimited\n')
        f.write('export OMP_NUM_THREADS=8\n')
        f.write('\n')
        f.write('$ORCAHOME/bin/orca ' + basename + '.in > $SLURM_SUBMIT_DIR/' + basename + '.out\n')
        f.write('mv * $SLURM_SUBMIT_DIR/scr\n')
        f.write('\n')
        f.close()

    #Generating inputs
    with open(basename + '/' + basename + '.in', 'w') as f:
        f.write('! uks B3LYP D3BJ Opt def2-SVP\n')
        f.write('%method\n')
        f.write('ScalHFX = 0.15\n') #0.15a, 0.77b, 0.81c, VMN5(default)
        f.write('ScalDFX = 0.77\n')
        f.write('ScalDFC = 0.81\n')
        f.write('ScalLDAC = 1\n')
        f.write('end\n')
        f.write('%pal nprocs 8 end\n')
        if dist_constraint != None: #Geometry optimization with distance constraint (COGEF steps)
            f.write('%geom\n')
            f.write('Constraints\n')
            f.write('{B ' + str(dist_constraint[0]) + ' ' + str(dist_constraint[1]) + ' '+ str(dist_constraint[2]) + ' C}\n')
            f.write('end\n')
            f.write('end\n')
        elif EFEI != None: #Geometry optimization under force (EFEI)
            f.write('%geom POTENTIALS\n')
            f.write('{C ' + str(EFEI[0]) + ' ' + str(EFEI[1]) + ' ' + str(EFEI[2]) + '}\n')
            f.write('end\n')
            f.write('end\n')
        f.write('* xyzfile ' + str(int(charge)) + ' ' + str(spinval) + ' ' + basename + '.xyz\n')
        f.close()

    #Generating xyz file
    molecule = mol3D()
    molecule.readfrommol2(mol2,readstring=True)
    molecule.writexyz(basename + '/' + basename + '.xyz',symbsonly=True)


# In[5]:


df = pd.read_csv('1stround_recorded_corrected.csv')


# In[19]:


for i in np.arange(df.shape[0]):
    
    row = df.iloc[i]
    refcode = row['refcode']
    metal = row['metal']
    charge = row['ox_csd']
    appair = None
    if row['apnum'] > 1:
        appair = [row['ap1'],row['ap2']]
    
    is_co2 = False #co2 doesn`t have IS state
    if metal == 'Co' and charge == 2:
        is_co2 = True
        
    pp1,pp2 = row['pp1'],row['pp2'] #Orca and csv(molSimplify): same indexing system starting from 0
    mol2ls = row['mol2LS_0nN']
    mol2hs = row['mol2HS_0nN']
    if not is_co2:
        mol2is = row['mol2IS_0nN']
        
    #Determine which ones have been performed
    ls1,is1,hs1,ls3,is3,hs3 = False,False,False,False,False,False
    if row['ELS_1nN'] != 'Unperformed' and row['ELS_1nN'] != 'Failed':
        ls1 = True
    if row['EIS_1nN'] != 'Unperformed' and row['EIS_1nN'] != 'Failed':
        is1 = True
    if row['EHS_1nN'] != 'Unperformed' and row['EHS_1nN'] != 'Failed':
        hs1 = True
    if row['ELS_3nN'] != 'Unperformed' and row['ELS_3nN'] != 'Failed':
        ls3 = True
    if row['EIS_3nN'] != 'Unperformed' and row['EIS_3nN'] != 'Failed':
        is3 = True
    if row['EHS_3nN'] != 'Unperformed' and row['EHS_3nN'] != 'Failed':
        hs3 = True

    if metal == 'Fe' and charge == 2:
        #1nN
        if ls1 == False:
            prep_orca_input(refcode,charge,'LS',1,mol2ls,multiPP=appair,EFEI=[pp1,pp2,1])
        if is1 == False:
            prep_orca_input(refcode,charge,'IS',3,mol2ls,multiPP=appair,EFEI=[pp1,pp2,1])
        if hs1 == False:
            prep_orca_input(refcode,charge,'HS',5,mol2ls,multiPP=appair,EFEI=[pp1,pp2,1])
        #3nN
        if ls3 == False:
            prep_orca_input(refcode,charge,'LS',1,mol2ls,multiPP=appair,EFEI=[pp1,pp2,3])
        if is3 == False:
            prep_orca_input(refcode,charge,'IS',3,mol2ls,multiPP=appair,EFEI=[pp1,pp2,3])
        if hs3 == False:
            prep_orca_input(refcode,charge,'HS',5,mol2ls,multiPP=appair,EFEI=[pp1,pp2,3])
    
    elif metal == 'Fe' and charge == 3:
        #1nN
        if ls1 == False:
            prep_orca_input(refcode,charge,'LS',2,mol2ls,multiPP=appair,EFEI=[pp1,pp2,1])
        if is1 == False: 
            prep_orca_input(refcode,charge,'IS',4,mol2ls,multiPP=appair,EFEI=[pp1,pp2,1])
        if hs1 == False:    
            prep_orca_input(refcode,charge,'HS',6,mol2ls,multiPP=appair,EFEI=[pp1,pp2,1])
        #3nN
        if ls3 == False:
            prep_orca_input(refcode,charge,'LS',2,mol2ls,multiPP=appair,EFEI=[pp1,pp2,3])
        if is3 == False:
            prep_orca_input(refcode,charge,'IS',4,mol2ls,multiPP=appair,EFEI=[pp1,pp2,3])
        if hs3 == False:
            prep_orca_input(refcode,charge,'HS',6,mol2ls,multiPP=appair,EFEI=[pp1,pp2,3])
    
    elif metal == 'Co' and charge == 2:
        #1nN
        if ls1 == False:
            prep_orca_input(refcode,charge,'LS',2,mol2ls,multiPP=appair,EFEI=[pp1,pp2,1])
        if hs1 == False:
            prep_orca_input(refcode,charge,'HS',4,mol2ls,multiPP=appair,EFEI=[pp1,pp2,1])
        #3nN
        if ls3 == False:
            prep_orca_input(refcode,charge,'LS',2,mol2ls,multiPP=appair,EFEI=[pp1,pp2,3])
        if hs3 == False:
            prep_orca_input(refcode,charge,'HS',4,mol2ls,multiPP=appair,EFEI=[pp1,pp2,3])
    
    elif metal == 'Co' and charge == 3:
        #1nN
        if ls1 == False:
            prep_orca_input(refcode,charge,'LS',1,mol2ls,multiPP=appair,EFEI=[pp1,pp2,1])
        if is1 == False:
            prep_orca_input(refcode,charge,'IS',3,mol2ls,multiPP=appair,EFEI=[pp1,pp2,1])
        if hs1 == False:
            prep_orca_input(refcode,charge,'HS',5,mol2ls,multiPP=appair,EFEI=[pp1,pp2,1])
        #3nN
        if ls3 == False:
            prep_orca_input(refcode,charge,'LS',1,mol2ls,multiPP=appair,EFEI=[pp1,pp2,3])
        if is3 == False:
            prep_orca_input(refcode,charge,'IS',3,mol2ls,multiPP=appair,EFEI=[pp1,pp2,3])
        if hs3 == False:
            prep_orca_input(refcode,charge,'HS',5,mol2ls,multiPP=appair,EFEI=[pp1,pp2,3])
    
    elif metal == 'Mn' and charge == 2:
        #1nN
        if ls1 == False:
            prep_orca_input(refcode,charge,'LS',2,mol2ls,multiPP=appair,EFEI=[pp1,pp2,1])
        if is1 == False:
            prep_orca_input(refcode,charge,'IS',4,mol2ls,multiPP=appair,EFEI=[pp1,pp2,1])
        if hs1 == False:
            prep_orca_input(refcode,charge,'HS',6,mol2ls,multiPP=appair,EFEI=[pp1,pp2,1])
        #3nN
        if ls3 == False:
            prep_orca_input(refcode,charge,'LS',2,mol2ls,multiPP=appair,EFEI=[pp1,pp2,3])
        if is3 == False:
            prep_orca_input(refcode,charge,'IS',4,mol2ls,multiPP=appair,EFEI=[pp1,pp2,3])
        if hs3 == False:
            prep_orca_input(refcode,charge,'HS',6,mol2ls,multiPP=appair,EFEI=[pp1,pp2,3])







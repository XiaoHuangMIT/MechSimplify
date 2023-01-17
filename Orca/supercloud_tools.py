import os
import sys
import numpy as np
import pandas as pd
from molSimplify.Classes.mol3D import mol3D
from molSimplify.Classes.atom3D import atom3D #distance
from molSimplify.job_manager.manager_io import read_outfile

def prep_orca_input(refcode, charge, spin, spinval, mol2, metal = None, multiPP = None, dist_constraint = None, EFEI = None):


    #Write orca EFEI geometry optimization input, including distance-constraint and force_modified variations, that
    #can be run on Expanse/Slurm and used with job_manager
    #Method by default:
    #B3LYP*(VMN5)-D3BJ: 0.15 a, 0.77b, 0.81c
    #Basis: def2-svp
    #If metal is specified: used def2-tzvp
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
    basename = refcode 
    if multiPP != None:
        basename = basename + '_' + str(multiPP[0]) + '_' + str(multiPP[1])
    if dist_constraint != None:
        basename = basename + '_' + str(dist_constraint[2])
    elif EFEI != None:
        basename = basename + '_' + str(EFEI[2])
    basename = basename + '_' + spin

    #Make folder
    if not os.path.exists(basename + '/'):
        os.mkdir(basename + '/')

    #Generating jobscript:
    with open(basename + '/' + basename + '_jobscript', 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('#SBATCH --job-name=' + basename + '\n')
        f.write('#SBATCH --output=' + basename + '.out' + '\n')
        f.write('#SBATCH --mem-per-cpu=4gb\n')
        f.write('#SBATCH --time=48:00:00\n')
        f.write('#SBATCH --ntasks=8\n')
        f.write('#SBATCH --export=*.xyz\n')
        f.write('\n')
        f.write('module load mpi/openmpi-4.1.1\n')
        f.write('export PATH=/home/gridsan/xiaohuang/orca:$PATH\n')
        f.write('export LD_LIBRARY_PATH=/home/gridsan/xiaohuang/orca:$LD_LIBRARY_PATH\n')
        f.write('\n')
        f.write('# Change to working directory\n')
        f.write('export SCRDIR=/scratch/$USER/job_$SLURM_JOB_ID\n')
        #f.write('cp $SLURM_SUBMIT_DIR/' + basename + '.in ' + '$SCRDIR\n') check again
        f.write('cp $SLURM_SUBMIT_DIR/* $SCRDIR\n')
        f.write('mkdir $SLURM_SUBMIT_DIR/scr\n')
        f.write('cd $SCRDIR\n')
        f.write('\n')
        f.write('export OMP_NUM_THREADS=8\n')
        f.write('\n')
        f.write('/home/gridsan/xiaohuang/orca/orca ' + basename + '.in > $SLURM_SUBMIT_DIR/' + basename + '.out\n')
        f.write('mv $SLURM_SUBMIT_DIR/scr\n')
        f.write('\n')
        f.close()

    #Generating inputs
    with open(basename + '/' + basename + '.in', 'w') as f:
        f.write('! uks B3LYP D3BJ Opt def2-SVP\n')
        f.write('%method\n')
        f.write('ScalHFX = 0.15\n') 
        f.write('ScalDFX = 0.77\n')
        f.write('ScalDFC = 0.81\n')
        f.write('ScalLDAC = 1\n')
        f.write('end\n')
        if metal != None:
            f.write('%basis\n')
            metal_basis = 'newgto ' + metal + ' "def2-TZVP" end\n'
            f.write(metal_basis)
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


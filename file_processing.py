#!/usr/bin/env python
# coding: utf-8

import os
import sys
import numpy as np
import pandas as pd
from molSimplify.Classes.mol3D import *
from molSimplify.job_manager.manager_io import read_outfile



######################################################################################################################
#########################################################Terachem#####################################################
######################################################################################################################


def prep_tera_input(refcode, spin, spinval, charge, mol2string, ap_pair=None):

    #Function to prepare terachem geometry optimization inputs which could be used with job_manager on Gibraltar
    #Method by default:
    #B3LYP*(VMN5)-D3BJ: 0.15 a, 0.77b, 0.81c
    #Basis: LACVP* (6-31g* and LANL2DZ), without spherical basis
    #Solvent: gaseous phase

    #In current directory, make a folder: refcode_(ap1ap2)_ls/is/hs
    #Inside the folder, write foldername_jobscript, foldername.in, foldername.xyz

    #Inputs:
    #refcode: six-letter CSD code of molecule
    #spin: string 'ls/is/hs'
    #spinval: integer value of spin
    #charge: integer value of charge
    #mol2string: mol2 string of the molecule
    #ap_pair: pair of terminal Hs to be replaced by ethyl, in length-2 list

    #Generate name of job
    if ap_pair == None:
        basename = refcode + '_' + spin
    else:
        ap1str,ap2str = str(ap_pair[0]), str(ap_pair[1])
        basename = refcode + '_' + ap1str + '_' + ap2str + '_' + spin

    #Make directory
    dirname = basename+'/'
    if not os.path.exists(dirname):
        os.mkdir(dirname)

    #Write jobscripts
    with open(dirname + basename + '_jobscript', 'w') as f:
        f.write('#$ -S /bin/bash\n')
        f.write('#$ -N ' + basename + '\n')
        f.write('#$ -cwd\n')
        f.write('#$ -R y\n')
        f.write('#$ -cwd\n')
        f.write('#$ -l h_rt=48:00:00\n')
        f.write('#$ -l h_rss=8G\n')
        f.write('#$ -q (gpusnew|gpus)\n')
        f.write('#$ -l gpus=1\n')
        f.write('#$ -pe smp 1\n')
        f.write('# -fin ' + basename + '.in\n')
        f.write('# -fin ' + basename + '.xyz\n')
        #f.write('# -fin pcm_radii\n')
        f.write('# -fout scr\n')
        f.write('#$ -cwd\n')
        f.write('\n')
        f.write('module load terachem\n')
        f.write('export OMP_NUM_THREADS=1\n')
        f.write('terachem ' + basename + '.in > $SGE_O_WORKDIR/' + basename + '.out')
        f.write('\n')
        f.close()

    #Write input files
    with open(dirname + basename + '.in', 'w') as f:
        f.write('run minimize\n')
        f.write('scf diis+a\n')
        f.write('coordinates ' + basename + '.xyz\n')
        f.write('levelshift yes\n')
        f.write('levelshiftvala 0.25\n')
        f.write('levelshiftvalb 0.25\n')
        f.write('spinmult ' + str(spinval) + '\n') #specify spin
        f.write('scrdir ./scr\n')
        f.write('basis lacvps_ecp\n')
        f.write('timings yes\n')
        f.write('charge ' + str(charge) + '\n') #specify charge
        f.write('method ub3lyp5\n')
        f.write('hfx 0.15\n') #optimal for spin-splitting
        f.write('precision dynamic\n')
        f.write('dftgrid 1\n')
        f.write('dynamicgrid yes\n')
        f.write('maxit 500\n')
        f.write('gpus 1\n')
        f.write('dispersion yes\n')
        f.write('new_minimizer yes\n')
        #f.write('pcm cosmo\n')b mostly non solvent or nonpolar solvent system
        #f.write('epsilon 78.39\n')
        #f.write('pcm_radii read\n')
        #f.write('pcm_radii_file /home/x_huang/useful_parameters/pcm_radii\n')
        f.write('nbo no\n')
        f.write('### props ###\n')
        f.write('ml_prop yes\n')
        f.write('poptype mulliken\n')
        f.write('bond_order_list no\n')
        f.write('end\n')
        f.close()

    #Write xyz files
    molecule = mol3D()
    molecule.readfrommol2(mol2string,readstring=True)
    molecule.writexyz(dirname + basename + '.xyz',symbsonly=True)

    
    
def find_opt_xyz(optim):
    
    #Getting last (optimized) geometry structure from terachem input
    #Input: pathway of optim.xyz output file generated by terachem geometry optimization
    #Output: generate a temp.xyz that stores the optimized structure
    
    file = open(optim, 'r')
    lines = file.readlines()
    
    #Count number of atoms of molecule
    num_atoms = lines[0].split()[0]
    num_atoms = int(num_atoms)
    
    #Count number of frames
    #Each frame num lines: num_atoms + 2 header
    num_lines = num_atoms+2
    num_frames = len(lines)/(num_lines)
    num_frames = int(num_frames)
    
    #Number of lines to read: 2 + number of atoms
    first_line = (num_frames - 1)  * num_lines
    lines_output = lines[first_line:]
    file2 = open('temp.xyz', 'w')
    file2.writelines(lines_output)
    file2.close()
    

    
def get_tera_opt_out(filepath):
    
    #Getting optimized energy and geometry structure from a terachem geom opt job
    #filepath: folder named X, containing X_jobscript, X.in, X.xyz, X.out and scr folder
    #dependency: find_opt_xyz; read_outfile
    
    if os.path.exists(filepath + '/scr') == False: #if not finished
        Eout = None
        mol2out = 'Failed Convergence'
    else:
        dict_out = read_outfile(filepath + '/' + filepath + '.out')
        Eout = dict_out['finalenergy']
        if Eout != None: 
            find_opt_xyz(filepath + '/scr/optim.xyz')
            mol3Dout = mol3D()
            mol3Dout.readfromxyz('temp_opt.xyz')
            mol2out = mol3Dout.writemol2('temp.mol',writestring = True)
        else: #if convergence failed
            mol2out = 'Failed Convergence'
    
    return Eout,mol2out    
 
    
    
######################################################################################################################
#########################################################Orca#########################################################
######################################################################################################################



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



def read_orca(filepath):

    #Check if optimization has converged normally (not inconvergence or out of time)
    #If so,return the final energy: energy + external potential

    if exists(filepath) == False:
        return 'Failed'

    lines = open(filepath,'r').readlines()
    cond_1 = False #'Optimization run done' is found in file
    cond_2 = False #'Orca terminated normally' is found (as second-last line)
    cond_3 = True # 'The optimization did not converge' is found in file

    Esp,Eep = [],[]
    for line in lines:

        if 'OPTIMIZATION RUN DONE' in line:
            cond_1 = True

        elif 'ORCA TERMINATED NORMALLY' in line:
            cond_2 = True

        elif 'The optimization did not converge' in line:
            cond_3 = False

        elif 'FINAL SINGLE POINT ENERGY' in line:
            Esp.append(float(line.split()[-1]))

        elif 'External Potential' in line:
            Eep.append(float(line.split()[3]))

    if cond_1 and cond_2 and cond_3:
        if len(Eep) != 0:#might be empty for 0nN
            return Esp[-1] + Eep[-1]
        else:
            return Esp[-1]
    else:
        return 'Failed'



def make_dir(dirname):

    #Helper function for avoiding errors while making a folder

    if not os.path.exists(dirname):
        os.makedirs(dirname)
    else:
        print('directories already exist')
        sys.exit()
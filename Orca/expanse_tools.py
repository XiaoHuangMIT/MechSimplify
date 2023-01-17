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
    #format: None/[pp1idx,pp2idx,force]

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
        f.write('#!/usr/bin/env bash\n')
        f.write('#SBATCH --job-name=' + basename + '\n')
        f.write('#SBATCH --output=' + basename + '.out' + '\n')
        f.write('#SBATCH --mem=16GB\n')
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
        f.write('mv *.xyz $SLURM_SUBMIT_DIR/scr\n')
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



def read_orca(filepath):

    #Check if optimization has converged normally (not inconvergence or out of time)
    #If so,return the final energy: energy + external potential

    if os.path.exists(filepath) == False:
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

    
    
def record_xyz(filepath):
    
    if os.path.exists(filepath) == False:
        return 'Failed'
    
    molecule = mol3D()
    molecule.readfromxyz(filepath)
    mol2 = molecule.writemol2('temp.mol',writestring = True)
    
    return mol2



def make_dir(dirname):

    #Helper function for avoiding errors while making a folder

    if not os.path.exists(dirname):
        os.makedirs(dirname)
    else:
        print('directories already exist')
        sys.exit()

        
        
 def analyze_spin(df,col1,col2,col3,name):
    
    #Analyze the ground, second and highest spin of molecules recorded in dataframe
    #col1,col2,col3: column name of low, intermediate and high spin
    #name: results summarized in column with name Esplitingname
    
    grounds,seconds,thirds,Esps = [],[],[],[]

    for i in np.arange(df.shape[0]):
    
        row = df.iloc[i]
        metal = row['metal']
        charge = row['ox_csd']
    
        if metal == 'Co' and charge == 2.0:
            Els = float(row[col1])
            Ehs = float(row[col3])
            if Els < Ehs:
                ground,second,third = 'LS','HS','N/A'
                Esp = (Ehs - Els) * 627.509
            else:
                ground,second,third = 'HS','LS','N/A'
                Esp = (Els - Ehs) * 627.509
        
    
        else:
            Els = float(row[col1])
            Eis = float(row[col2])
            Ehs = float(row[col3])
            df1 = pd.DataFrame()
            df1['Spin'] = ['LS','IS','HS']
            df1['E'] = [Els,Eis,Ehs]
            df1 = df1.sort_values(by = 'E')
            ground,second,third = df1['Spin'].iloc[0],df1['Spin'].iloc[1],df1['Spin'].iloc[2]
            Esp = (df1['E'].iloc[1] - df1['E'].iloc[0]) * 627.509
    
        grounds.append(ground)
        seconds.append(second)
        thirds.append(third)
        Esps.append(Esp)

    df['ground_spin' + name] = grounds
    df['second_spin' + name] = seconds
    df['third_spin' + name] = thirds
    df['Esplitting' + name] = Esps

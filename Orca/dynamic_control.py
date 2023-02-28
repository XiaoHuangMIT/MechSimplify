import os
import sys
import numpy as np
import pandas as pd
from molSimplify.Classes.mol3D import mol3D
from molSimplify.Classes.ligand import ligand_breakdown
from molSimplify.job_manager.tools import call_bash, list_active_jobs


#Function to find optimized xyzs from orca out file
def find_opt_frames(filepath,natoms):

    if os.path.exists(filepath) == False:
        return 'Failed'

    lines = open(filepath,'r').readlines()
    xyzs = []
    for linenum in np.arange(len(lines)):
        if 'CARTESIAN COORDINATES (ANGSTROEM)' in lines[linenum]: #xyz beginning two lines later
            firstline = linenum + 2
            lastline = firstline + natoms - 1
            xyz = []
            for i in np.arange(firstline,lastline+1):
                xyz.append(lines[i].lstrip())
            xyzs.append(xyz)

    if len(xyzs) > 0:
        xyzs.pop()#n cycles: n+1 structures, so drop last one

    return xyzs


#Function to examine whether the job has reached dissociation by reading the .out file
def check_diss_by_out(basename,threshold = 10,return_list = False):

    natoms = open(basename + '/' + basename + '.xyz','r').readlines()[0].split()[0]
    natoms = int(natoms)

    frames = find_opt_frames(basename + '/' + basename + '.out', natoms)
    if threshold > len(frames):
        return False

    checks = []
    nframes = len(frames)
    first_frame = nframes - threshold + 1

    for i in np.arange(first_frame,nframes+1):
        file = open('temp_n.xyz', 'w')
        file.writelines(str(natoms) + '\n')
        file.writelines('\n')
        file.writelines(frames[i-1])
        file.close()

        mol = mol3D()
        mol.readfromxyz('temp_n.xyz')
        l1,l2,l3 = ligand_breakdown(mol)
        check = 'intact'
        if l2 != [3,3]:
            check = 'diss'
        checks.append(check)

    result = True
    for check in checks:
        if check == 'intact':
            result = False

    if return_list:
        return result,checks
    else:
        return result


#Get name of active jobs and their idxs
names,ids = list_active_jobs(ids=True)

#The name of job should correspond to the folder that job inputs are located
#Double check while using this script
num_canceled = 0
for i in np.arange(len(names)):

    basename = names[i]
    diss = check_diss_by_out(basename)
    #print(basename, diss)

    #If diss = True, kill the job
    if diss == True:
        call_bash('scancel ' + str(ids[i]))
        print('Killed ', basename, diss)
        num_canceled += 1

print(num_canceled, ' jobs killed')

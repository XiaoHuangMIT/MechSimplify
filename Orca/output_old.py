import os
import sys
import numpy as np
import pandas as pd
from molSimplify.Classes.mol3D import mol3D
from molSimplify.Classes.ligand import ligand_breakdown

 def find_nth_xyz(filepath,nth):

    #Find nth structure from a scan/trajectory xyz file
    #Stores as an temp.xyz file

    file = open(filepath, 'r')
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
    first_line = (nth - 1)  * num_lines
    lines_output = lines[first_line:first_line+num_lines]
    file2 = open('temp_nth.xyz', 'w')
    file2.writelines(lines_output)
    file2.close()

    

def count_num_frames(filepath):

    #Count how many structures are in a scan/trajectory xyz file
    #Returned running: if xyz is being written right now

    file = open(filepath, 'r')
    lines = file.readlines()

    #Count number of atoms of molecule
    num_atoms = lines[0].split()[0]
    num_atoms = int(num_atoms)

    #Count number of frames
    #Each frame num lines: num_atoms + 2 header
    num_lines = num_atoms+2
    num_frames = len(lines)/(num_lines)
    if int(num_frames) == num_frames:
        return int(num_frames)
    else:
        return 'running'


    
def check_dissociation(filepath,threshold = 5,return_list=False):

    #Check if we can conclude that an EFEI opt resulted in dissociation
    #If for the last num threshold jobs, the molecule has been dissociated as indicated by molSimplify
    #we can then end the jobd
    #Default threshold is 5

    nframes = count_num_frames(filepath)
    if threshold > nframes:
        return 'Not enough frames'

    frame_nums = np.arange(nframes-threshold+1, nframes+1)
    checks = []
    for i in frame_nums:
        find_nth_xyz(filepath,i)
        mol = mol3D()
        mol.readfromxyz('temp_nth.xyz')
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

import os
import sys
import numpy as np
import pandas as pd
from os.path import exists
from molSimplify.Classes.mol3D import mol3D
from molSimplify.Classes.ligand import ligand_breakdown



def analyze_scan_optim(filepath): #completed
    
    #Read a scan_optim output file of a COGEF run
    #Returns list of distance, list of energy, list of mol2
    
    file = open(filepath, 'r')
    lines = file.readlines()
    
    #Find Energies and Distances
    energies, distances = [],[]
    for line in lines:
        if 'Converged' in line:
            energies.append(float(line.split()[4]))
            distance = line.split()[6]
            distance = distance.replace('(','')
            distance = distance.replace(')','')
            distances.append(float(distance))
    
    
    #Count number of atoms of molecule
    num_atoms = lines[0].split()[0]
    num_atoms = int(num_atoms)
    
    #Count number of frames
    #Each frame num lines: num_atoms + 2 header
    num_lines = num_atoms+2
    num_frames = len(lines)/(num_lines)
    num_frames = int(num_frames)
    
    #Storing mol2s
    mol2s = []
    for i in np.arange(num_frames):
        first_line = i * num_lines
        last_line = first_line + (num_atoms + 2) - 1 #Each xyz has num_atoms+2 lines
        lines_output = lines_output = lines[first_line:last_line+1]
        file2 = open('temp.xyz', 'w')
        file2.writelines(lines_output)
        file2.close()
        moltemp = mol3D()
        moltemp.readfromxyz('temp.xyz')
        mol2 = moltemp.writemol2('temp.mol',writestring = True)
        mol2s.append(mol2)
    
    return distances, energies, mol2s
  
  
  
  def has_dissociated(mol2):
    
    #Reads a mol2 String and determine if it has dissociated
    molecule = mol3D()
    molecule.readfrommol2(mol2,readstring=True)
    idxs,dents,coords = ligand_breakdown(molecule)
    if dents == [3,3]:
        return False
    else:
        return True
      
      
      
def calculate_force(Es,dist=0.2):
    
    #Reads an scan_optim energy list and calculate the cogef force
    #dist: stretching distance
    
    Erels = (np.array(Es) + Es[0]) * 2625.5 #in kJ
    Fs = []
    for i in 1 + np.arange(len(Erels) - 1):
        F = (Erels[i] - Erels[i-1]) / dist #in kJ/A
        F = F*1.66/100 #in nN
        Fs.append(F)
    return max(Fs)



def continue_COGEF(oldpath,new_parent_path,walltime = 'same',stretch_dist = 10, dist_step = 0.2):
    
    #Generate new directory, input, jobscript and xyz with same name in the new parent folder
    #So that job manager could be run in the new parent path without re-running finished jobs
    #Input will be modified as continuing of COGEF
    #Tip: cannot change strech distance or step distance for a job that starts from scratch
    
    #Read scan_optim
    from_scratch = False
    old_scan_path = oldpath + '/scr/scan_optim.xyz'
    if exists(old_scan_path) == False or os.stat(old_scan_path).st_size == 0:
        from_scratch = True
    else:
        old_scan_file = open(old_scan_path, 'r')
        old_scan_lines = old_scan_file.readlines()
        distances = []
        for line in old_scan_lines:
            if 'Converged' in line:
                distance = line.split()[6]
                distance = distance.replace('(','')
                distance = distance.replace(')','')
                distances.append(float(distance))
        d0,dt = distances[0], distances[-1] #beginning frame and last frame finished
    
    #Make directory for continued COGEF run
    new_path = new_parent_path + '/' + oldpath
    if os.path.exists(new_parent_path) == False: #make things easier
        os.mkdir(new_parent_path)
    if os.path.exists(new_path) == False:
        os.mkdir(new_path)
    
    #Write jobscript
    old_js = open(oldpath + '/' + oldpath + '_jobscript')
    js_lines = old_js.readlines()
    if walltime != 'same': #if we want to change walltime
        for idx in np.arange(len(js_lines)):
            if '#$ -l h_rt' in js_lines[idx]:
                js_lines[idx] = '#$ -l h_rt' + ' '+ walltime + '\n'
    new_js = new_path + '/' + oldpath + '_jobscript'
    with open(new_js, 'w') as f:
        for line in js_lines:
            f.write(str(line))
        f.close()
    
    #Write input
    if from_scratch == False:
        dini, dfin = dt + dist_step, d0 + stretch_dist #beginning and end point for new COGEF
        nsteps = ((dfin - dini) / dist_step) + 1
        nsteps = int(nsteps)
    old_in = open(oldpath + '/' + oldpath + '.in')
    in_lines = old_in.readlines()
    if from_scratch == False:
        for idx in np.arange(len(in_lines)):
            if 'bond' in in_lines[idx]:
                pp_pair = in_lines[idx].split()[-1]
                in_lines[idx] = 'bond ' + str(dini) + ' ' + str(dfin) + ' ' + str(nsteps) + ' ' + pp_pair + '\n'
    new_in = new_path + '/' + oldpath + '.in'
    with open(new_in, 'w') as f:
        for line in in_lines:
            f.write(str(line))
        f.close()
    
    #Write(Copy) xyz
    old_xyz = open(oldpath + '/' + oldpath + '.xyz')
    xyz_lines = old_xyz.readlines()
    new_xyz = new_path + '/' + oldpath + '.xyz'
    with open(new_xyz,'w') as f:
        for line in xyz_lines:
            f.write(str(line))
        f.close()

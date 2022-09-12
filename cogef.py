import os
import sys
import numpy as np
import pandas as pd
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

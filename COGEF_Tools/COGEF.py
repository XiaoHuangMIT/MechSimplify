import os
import sys
import numpy as np
import pandas as pd
from os.path import exists
from molSimplify.Classes.mol3D import mol3D
from molSimplify.Classes.ligand import ligand_breakdown


#General principle: don`t add bunch of try/except (shall clean them later)
#Use helper function to first determine/eliminate extreme cases


def count_steps(filepath):
    
    #Analyze number of steps a COGEF run had proceeded
    #Return 0 as well if no file was found
    #Tip: we count 0A optimization as 1st step as well, so 10A stretching would be 51 steps
    
    try:
        file = open(filepath,'r')
    except:
        return 0
    lines = file.readlines()
    
    count = 0
    for line in lines:
        if 'Converged' in line:
            count = count + 1
    return count



def distance3d(coord1,coord2): 
    
    #A helper function to calculate distance between two 3d points (as list of floats)
    
    coord1 = np.array(coord1)
    coord2 = np.array(coord2)
    squared_dist = np.sum((coord1-coord2)**2, axis=0)
    dist = np.sqrt(squared_dist)
    return dist



def distance_from_xyz(lines,atom1,atom2,frontlines=2):
    
    #A helper function to calculate distance between two atoms with
    #lines: lines of xyz file, stored as list
    #atom1, atom2: index of atoms that START WITH ZERO
    #frontlines: number of additional lines in front, default 2 for terachem optimizations
    
    #0th atom is 3rd line/lines[2], so nth atom is lines[n+2]
    line1,line2 = lines[atom1+2].split(),lines[atom2+2].split()
    x1,y1,z1 = float(line1[1]), float(line1[2]), float(line1[3])
    x2,y2,z2 = float(line2[1]), float(line2[2]), float(line2[3])
    
    dist = distance3d([x1,y1,z1],[x2,y2,z2])
    return dist



def analyze_scan_optim(filepath,no_mol2 = False): 
    
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
    if no_mol2 == False:
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
    
    if no_mol2 == False:
        return distances, energies, mol2s
    else:
        return distances, energies
  
  
  
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



def coord_number_analysis(filepath,threshold=1.5):
    
    #Analyze the coordination number (how many coord bonds were still there) of each structure during COGEF
    #If a bond had become longer than threshold * original length: consider it broken
    #To save time: does not involve extensive calling of mol3D
    #(currently: could only work with octahedral/6-coord complexes)
    #Inputs:
    #filepath: path to scan_optim.xyz
    #Outputs:
    #A list of coordination numbers
    
    
    file = open(filepath, 'r')
    lines = file.readlines()
    
    #Count number of atoms of molecule, thus number of lines in a xyz and number of total structures
    num_atoms = lines[0].split()[0]
    num_atoms = int(num_atoms)
    num_lines = num_atoms+2 
    num_frames = len(lines)/(num_lines)
    num_frames = int(num_frames)
    
    #Summary all xyz into a list of lists, each list correspond to lines of single xyz
    xyzs = [] 
    for i in np.arange(num_frames):
        first_line = i * num_lines
        last_line = first_line + (num_atoms + 2) - 1 #Each xyz has num_atoms+2 lines
        lines_each_xyz = lines[first_line:last_line+1]
        xyzs.append(lines_each_xyz)
    
    #Generate a mol3D instance to find index of metal and coordinating atoms
    file2 = open('temp.xyz', 'w')
    file2.writelines(xyzs[0])
    file2.close()
    moltemp = mol3D()
    moltemp.readfromxyz('temp.xyz')
    liglist,ligdent,ligcons = ligand_breakdown(moltemp)
    coords = ligcons[0] + ligcons[1]
    coords.sort() #index of coordinating atoms
    metal = moltemp.findMetal()[0] #index of metal
    if len(coords) != 6: #If ligand breakdown didnt work well
        return 'Failed'
    
    #Find original bond lengths of six bonds
    bonds_0 = []
    for idx in coords:
        l = distance_from_xyz(xyzs[0],metal,idx)
        bonds_0.append(l)
    
    #Evaluate for each bond at each structure: had it become longer than threshold * original length
    bond_orders = []
    for job in np.arange(len(xyzs)):
        order = 0
        bonds = []
        for idx in coords: #coord bond length of each structure
            l = distance_from_xyz(xyzs[job],metal,idx)
            bonds.append(l)
        for i in np.arange(6): #compare each bond length
            if bonds[i] < bonds_0[i] * threshold:
                order += 1
        bond_orders.append(order)
    
    return bond_orders



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

        
        
def iters_each_step(filepath):
    
    #Analyze for each COGEF scan step, how many optimization steps were perfromed
    #Input:
    #filepath: path to optim.xyz
    #Output:
    #list of num frames for each structure
    
    try:
        file = open(filepath,'r')
    except:
        return 'No File' #File doesnt exist
    lines = file.readlines()
    if len(lines) < 2:
        return 'Empty File' #File is empty
    
    nums = []
    frames = []

    for line in lines:
        if 'frame' in line:
            nums.append(int(line.split()[2]))
    
    for i in np.arange(1,len(nums)):
        if nums[i] == 0:
            frames.append(nums[i-1] + 1) #since frame numbers start with 0
    
    return frames



def plot_nsteps_vs_dist(lfs,dist=10,sep=0.2,cutoff=0):
    
    #For each molecule, plot the number of structures generated at optimization at each distance during COGEF
    #Also plot the sum of number of structures durinng each optimization at distances for all molecules
    #Return the sum num steps for all structures for further analysis
    #Inputs:
    #lfs: list of list, every element correspond to a list containing number of frames during each optimization step
    #at each distance for every molecule
    #dist and sep: total distance and distance at each step for COGEF
    #cutoff: only analyze molecules that undergo at least certain num of steps
    #Outputs:
    #Two plots, and a list of the sum num steps for all structures
    
    plt.figure()
    plt.xlabel('Stretching Distance (A)')
    plt.ylabel('# Structures')
    if cutoff > 0:
        title = 'Above ' + str(cutoff*sep) + 'A only'
        plt.title(title)
    for nframes in lfs:
        if type(nframes) != str and len(nframes) > cutoff:
            xs = np.arange(len(nframes)) * 0.2
            nframes = np.array(nframes)
            plt.plot(xs,nframes)
    plt.show()

    plt.figure()
    plt.xlabel('Stretching Distance (A)')
    plt.ylabel('Sum # Structures')
    if cutoff > 0:
        plt.title(title)
    sumnums = []
    for i in np.arange(int(dist/sep) + 1):
        count = 0
        for nframes in lfs:
            if type(nframes) != str and len(nframes) > i and len(nframes) > cutoff:
                count += int(nframes[i])
        sumnums.append(count)
    plt.plot(np.arange(0,dist+sep,sep),sumnums)

    
    
def plot_pnsteps_vs_dist(lfs,dist=10,sep=0.2,cutoff=0):
    
    #For each molecule, plot the number of structures generated at optimization at each distance during COGEF divided
    #by total number of structures generated at all distances (portion)
    #Also plot the sum of number of structures durinng each optimization at distances for all molecules
    #Return the sum portion num steps for all structures for further analysis
    #Inputs:
    #lfs: list of list, every element correspond to a list containing number of frames during each optimization step
    #at each distance for every molecule
    #dist and sep: total distance and distance at each step for COGEF
    #cutoff: only analyze molecules that undergo at least certain num of steps
    #Outputs:
    #Two plots, and a list of the sum portion num steps for all structures
    
    plt.figure()
    plt.xlabel('Stretching Distance (A)')
    plt.ylabel('# Portion of structures')
    if cutoff > 0:
        title = 'Above ' + str(cutoff*sep) + 'A only'
        plt.title(title)
    for nframes in lfs:
        if type(nframes) != str and len(nframes) > cutoff:
            xs = np.arange(len(nframes)) * 0.2
            nframes = np.array(nframes)
            plt.plot(xs,nframes/sum(nframes))
    plt.show()

    plt.figure()
    plt.xlabel('Stretching Distance (A)')
    plt.ylabel('Sum # portion of structures')
    if cutoff > 0:
        plt.title(title)
    sumnums = []
    for i in np.arange(int(dist/sep) + 1):
        count = 0
        for nframes in lfs:
            if type(nframes) != str and len(nframes) > i and len(nframes) > cutoff:
                count += nframes[i]/sum(nframes)
        sumnums.append(count)
    plt.plot(np.arange(0,dist+sep,sep),sumnums)
    



    
 

import os
import numpy as np
from numpy.random import randint
import pandas as pd
import pandas as pd
import matplotlib.pyplot as plt


def prep_aismd_input(name, reps, steering_file, xyz_file, rep_list=False, walltime='120:00:00',nstep=20000):

    #Function to preparae aimsd input for terachem for organic molecules
    #Method by default:
    #B3LYP(VMN5)-D3BJ
    #6-31gs
    #Gaseous phase

    #In current directory, make a folder: name_reps
    #Inside the folder, write name_reps_jobscript, name_reps.in
    #The folder should contain a steering file and a xyz file, to be then copied to every folder of repeated runs

    #Inputs:
    #name: name of the job
    #reps: number of repeated trials to run
    #reps can be specified as a list if reps_list is True
    #steering_file, xyz_file: name of steering txt and xyz file to be copied
    #walltime: max length of simulation time, default 120 hrs/5 days
    #nstep: max number of aismd simulation time, default 20000 steps/5 ps

    if rep_list == False:
        ljobs = np.arange(reps)
    elif rep_list == True:
        ljobs = reps
    ljobs = np.array(ljobs).as_type('int')
    ljobs = ljobs -1
    
    for i in ljobs:

        #Generate name of job
        if i < 10:
            i = '0' + str(i)
        else:
            i = str(i)
        basename = name + '_' + i  
        i += 1
        
        #Make directory
        dirname = basename + '/'
        if not os.path.exists(dirname):
            os.mkdir(dirname)
            
        #Copy xyz and steering file
        shutil.copy(steering_file,dirname)
        shutil.copy(xyz_file,dirname)
        
        #Write jobscript
        rescale = int(nstep) + 1
        with open(dirname + basename + '_jobscript','w') as f:
            f.write('#S -S /bin/bash\n')
            f.write('#$ -N ' + basename + '\n')
            f.write('#$ -cwd\n')
            f.write('#$ -R y\n')
            f.write('#$ -l h_rt=' + walltime + '\n')
            f.write('#$ -l h_rss=8G\n')
            f.write('#$ -q (gpusnew|gpus)\n')
            f.write('#$ -l gpus=1\n')
            f.write('#$ -pe smp 1\n')
            f.write('# -fin ' + basename + '.in\n')
            f.write('# -fin ' + xyz_file + '\n')
            f.write('# -fin ' + steering_file + '\n')
            f.write('# -fout scr/\n')
            f.write('\n')
            f.write('module load cuda\n')
            f.write('module load terachem\n')
            f.write('\n')
            f.write('export OMP_NUM_THREADS=1\n')
            f.write('\n')
            f.write('terachem ' + basename + '.in > $SGE_O_WORKDIR/' + basename + '.out')
            f.close()
        
        with open(dirname + basename + '.in','w') as f:
            f.write('run md\n')
            f.write('coordinates ' + xyz_file + '\n')
            f.write('maxit 1001\n')
            f.write('nstep ' + str(nstep) + '\n')
            f.write('rescalefreq ' + str(rescale) + '\n')
            f.write('tinit 300\n')
            seed = randint(1000,high=100000)
            f.write('seed ' + str(seed) + '\n')
            f.write('integrator reversible_d\n')
            f.write('timestep 0.25\n')
            f.write('\n')
            f.write('steering ' + steering_file + '\n')
            f.write('\n')
            f.write('basis 6-31G*\n')
            f.write('method ub3lyp\n')
            f.write('dispersion yes\n')
            f.write('\n')
            f.write('charge 0\n')
            f.write('spinmult 1\n')
            f.write('\n')
            f.write('scf diis+a \n')
            f.write('levelshift yes\n')
            f.write('levelshiftvala 1.0\n')
            f.write('levelshiftvalb 0.0\n')
            f.write('maxit 500\n')
            f.write('precision dynamic\n')
            f.write('new_minimizer yes\n')
            f.write('\n')
            f.write('scrdir ./scr\n')
            f.write('Timings yes\n')
            f.write('gpus 1\n')
            f.write('end\n')
            f.close()


def find_coord(line,special_case = False):
    #return xyz coordinate of an atom in xyz file
    line_list = line.split()
    x = float(line_list[1])
    y = float(line_list[2])
    z = float(line_list[3])
    return [x,y,z]


def find_distance(a1,a2):
    #a1,a2: xyz coord of two atoms as recorded in length-3 list (in A)
    xdiff = a1[0] - a2[0]
    ydiff = a1[1] - a2[1]
    zdiff = a1[2] - a2[2]
    sum_squares = xdiff**2 + ydiff**2 + zdiff**2
    return sum_squares**0.5


def analyze_aismd_traj(filename,pltname):
    #Returns dataframe containing six bond lengths each frame, and frame number for easy plotting
    #Also plot
    
    with open(filename, "r") as f:
        lines = f.readlines() #97 lines per frame, 1st line number of atoms(95), 2nd line energy and frame number
        num_frames = len(lines)/97
        num_frames = int(num_frames) 
    
    d60s,d69s,d70s,d13s,d22s,d23s = [],[],[],[],[],[]

    for i in np.arange(num_frames):
   
        coord_1 = find_coord(lines[2 + i*97],special_case=True) # +1
        coord_60 = find_coord(lines[61 + i*97])
        coord_69 = find_coord(lines[70 + i*97])
        coord_70 = find_coord(lines[71 + i*97])
        coord_13 = find_coord(lines[14 + i*97])
        coord_22 = find_coord(lines[23 + i*97])
        coord_23 = find_coord(lines[24 + i*97])
    
        d60 = find_distance(coord_1,coord_60)
        d69 = find_distance(coord_1,coord_69)
        d70 = find_distance(coord_1,coord_70)
        d13 = find_distance(coord_1,coord_13)
        d22 = find_distance(coord_1,coord_22)
        d23 = find_distance(coord_1,coord_23)
    
        d60s.append(d60)
        d69s.append(d69)
        d70s.append(d70)
        d13s.append(d13)
        d22s.append(d22)
        d23s.append(d23)
    
    df = pd.DataFrame()
    df['d60'] = d60s
    df['d69'] = d69s
    df['d70'] = d70s
    df['d13'] = d13s
    df['d22'] = d22s
    df['d23'] = d23s
    
    distances = np.arange(0,num_frames)*0.25
    plt.figure(figsize=(8,6))
    plt.plot(distances,d69s,label='N1-1(c)')
    plt.plot(distances,d60s,label='N1-2')
    plt.plot(distances,d70s,label='N1-3')
    plt.plot(distances,d22s,label='N2-1(c)')
    plt.plot(distances,d13s,label='N2-2')
    plt.plot(distances,d23s,label='N2-3')
    plt.xlabel('Time (fs)')
    plt.ylabel('Bond Length (A)')
    plt.legend()
    plt.title(pltname)
    plt.show()
    
    return df,num_frames

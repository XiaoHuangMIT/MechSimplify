import os
import numpy as np
import pandas as pd
import pandas as pd
import matplotlib.pyplot as plt


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

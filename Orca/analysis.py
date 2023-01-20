import os
import sys
import numpy as np
import pandas as pd
from molSimplify.Classes.mol3D import mol3D
from molSimplify.Classes.atom3D import atom3D #distance
from molSimplify.job_manager.manager_io import read_outfile



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

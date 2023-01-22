import os
import sys
import numpy as np
import pandas as pd
from molSimplify.Classes.mol3D import mol3D
from molSimplify.Classes.atom3D import atom3D #distance
from molSimplify.job_manager.manager_io import read_outfile



#########################################################################################################################################################
##########################################################Functions Used#################################################################################
#########################################################################################################################################################



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



def analyze_efei_expanse(basename):

    #Analyze an EFEI job being performed on expanse given the current configuration:
    #basename/basename.out, scr/basename.xyz(being overwritten after optimization), basename_trj.xyz
    #basename: refcode_aps_LS/IS/HS
    #First check: has the job succeeded. If so, record basename.xyz as mol2 of optimized molecule
    #Then check: has the molecule broken.
    #Lastly,if broken and job failed, change  energy and mol2 from 'Failed' to 'Diss'

    energy = read_orca(basename + '/' + basename + '.out')
    mol2 = 'Failed'
    if energy != 'Failed':
        mol2 = record_xyz(basename + '/scr/' + basename + '.xyz')

    has_dis = False
    if path.exists(basename + '/scr/' + basename + '_trj.xyz'):
        has_dis = check_dissociation(basename + '/scr/' + basename + '_trj.xyz')

    if has_dis == True:
        if energy == 'Failed':
            energy = 'Diss'
        if mol2 == 'Failed':
            mol2 = 'Diss'
    if has_dis == 'Not enough frames': #Makes cases less complex
        has_dis = False
    return energy, mol2, has_dis



def analyze_efei_supercloud(basename):

    #Analyze an EFEI job being performed on supercloud given the current configuration:
    #basename/basename.out, basename.xyz(being overwritten after optimization), basename_trj.xyz
    #basename: refcode_aps_LS/IS/HS
    #First check: has the job succeeded. If so, record basename.xyz as mol2 of optimized molecule
    #Then check: has the molecule broken.
    #Lastly,if broken and job failed, change  energy and mol2 from 'Failed' to 'Diss'

    energy = read_orca(basename + '/' + basename + '.out')
    mol2 = 'Failed'
    if energy != 'Failed':
        mol2 = record_xyz(basename + '/' + basename + '.xyz')

    has_dis = False
    if path.exists(basename + '/' + basename + '_trj.xyz'):
        has_dis = check_dissociation(basename + '/' + basename + '_trj.xyz')

    if has_dis == True:
        if energy == 'Failed':
            energy = 'Diss'
        if mol2 == 'Failed':
            mol2 = 'Diss'
    if has_dis == 'Not enough frames': #Makes cases less complex
        has_dis = False
    return energy, mol2, has_dis



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
    
    
    
#########################################################################################################################################################
#########################################################Example Scripts#################################################################################
#########################################################################################################################################################



r1_fs,E_ls,E_is,E_hs,mol2_ls,mol2_is,mol2_hs,dis_ls,dis_is,dis_hs = [],[],[],[],[],[],[],[],[],[]

indexes = df.index.values
for i in indexes:

    row = df.loc[i]
    refcode,ap1,ap2,apnum,metal,charge = row['refcode'],row['ap1'],row['ap2'],row['apnum'],row['metal'],row['ox_csd']
    aps = ''
    if float(row['apnum']) != 1:
        aps = '_' + str(ap1) + '_' + str(ap2)
    refname = refcode + aps + '_1'

    #Record all
    Els,mol2ls,disls = analyze_efei_supercloud(refname + '_LS')
    if metal == 'Co' and charge == 2:
        Eis,mol2is,disis = 'N/A','N/A','N/A'
    else:
        Eis,mol2is,disis = analyze_efei_supercloud(refname + '_IS')
    Ehs,mol2hs,dishs = analyze_efei_supercloud(refname + '_HS')
    
    #Resub analysis
    #if row['Els'] == 'Failed':
        #Els = read_orca(refname + '_LS' + aps + '/' + refname + '_LS' + aps + '.out')
        #mol2ls = record_xyz(refname + '_LS' + aps + '/scr/' + refname + '_LS' + aps + '.xyz')
        #df.at[i,'Els'] = Els
        #df.at[i,'mol2ls'] = mol2ls    

    #Record all recording    
    r1_fs.append(1)
    E_ls.append(Els)
    E_is.append(Eis)
    E_hs.append(Ehs)
    mol2_ls.append(mol2ls)
    mol2_is.append(mol2is)
    mol2_hs.append(mol2hs)
    dis_ls.append(disls)
    dis_is.append(disis)
    dis_hs.append(dishs)
    
df['round1_force'] = r1_fs
df['round1_Els'] = E_ls
df['round1_Eis'] = E_is
df['round1_Ehs'] = E_hs
df['round1_mol2ls'] = mol2_ls
df['round1_mol2is'] = mol2_is
df['round1_mol2hs'] = mol2_hs
df['round1_diss_ls'] = dis_ls
df['round1_diss_is'] = dis_is
df['round1_diss_hs'] = dis_hs

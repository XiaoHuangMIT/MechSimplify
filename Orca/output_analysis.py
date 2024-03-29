import os
import sys
import numpy as np
import pandas as pd
from molSimplify.Classes.mol3D import mol3D
from molSimplify.Classes.atom3D import atom3D #distance
from molSimplify.Classes.ligand import ligand_breakdown
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

    

def read_orca_Gibbs(filepath):

    #Check if optimization has converged normally (not inconvergence or out of time)
    #If so,return the final energy: energy + external potential

    if os.path.exists(filepath) == False:
        return 'Unperformed'

    lines = open(filepath,'r').readlines()
    cond_1 = False #'Optimization run done' is found in file
    cond_2 = False #'Orca terminated normally' is found (as second-last line)
    cond_3 = True # 'The optimization did not converge' is found in file

    for line in lines:

        if 'OPTIMIZATION RUN DONE' in line:
            cond_1 = True

        elif 'ORCA TERMINATED NORMALLY' in line:
            cond_2 = True

        elif 'The optimization did not converge' in line:
            cond_3 = False

        elif 'Final Gibbs free energy' in line:
            Eg = float(line.split()[-2])

    if cond_1 and cond_2 and cond_3:
        return Eg
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



def check_diss_by_out(basename,check_geom_shift = False, threshold = 10,return_list = False):\

    #If check_geom_shift = True:
    #Regard geometry shifting such as to tetrahedral etc also as dissociation
    #If not: only consider cases when only one ligand stay bonded as dissociation

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
        if check_geom_shift == False:
            if len(l2) < 2:
                check = 'diss'
        else:
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



def analyze_efei_expanse(basename, thres=10, record_failed_mol2 = True):

    #Analyze an EFEI job being performed on expanse given the current configuration:
    #basename/basename.out
    #basename: refcode_aps_LS/IS/HS
    #record_failed_mol2: if a job has failed (neither converge or dissociated)
    #we could record the last mol2 as re-run initial structure
    #First check: has the job succeeded. If so, record basename.xyz as mol2 of optimized molecule
    #Then check: has the molecule broken.
    #Lastly,if broken and job failed, change  energy and mol2 from 'Failed' to 'Diss'

    if os.path.exists(basename) == False: #Helper case: if the molecule has not been analyzed
        return 'Failed','Failed','Failed'
    
    #Record energy (from .out)
    energy = read_orca(basename + '/' + basename + '.out')
    
    #Check if the molecule has dissociated (don`t call an diss job 'failed')
    has_dis = check_diss_by_out(basename,threshold=thres)
    if has_dis == True and energy == 'Failed':
        energy = 'Diss'

    #Record mol2 (from .out)
    natoms = open(basename + '/' + basename + '.xyz','r').readlines()[0].split()[0]
    natoms = int(natoms)
    frames = find_opt_frames(basename + '/' + basename + '.out',natoms)
    if len(frames) < 2 and energy == 'Failed':
        mol2 = 'no_progress'
    else:
        file = open('temp_rec.xyz', 'w')
        file.writelines(str(natoms) + '\n')
        file.writelines('\n')
        file.writelines(frames[-1])
        file.close()
        mol_rec = mol3D()
        mol_rec.readfromxyz('temp_rec.xyz')
        mol2 = mol_rec.writemol2('temp_rec.mol',writestring = True)

    #Return results
    if record_failed_mol2 == False:
        if energy == 'Failed':
            return energy, 'Failed', has_dis
        elif energy == 'Diss':
            return energy, 'Diss', has_dis
    else:
        return energy, mol2, has_dis



def analyze_efei_supercloud(basename,thres=5):

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

    has_dis = check_diss_by_out(basename,threshold=thres)

    if has_dis == True:
        if energy == 'Failed':
            energy = 'Diss'
        if mol2 == 'Failed':
            mol2 = 'Diss'

    return energy, mol2, has_dis



def find_spin_delocalization(filename,metal):
    
    #Check if spin if significantly delocalized away from metal
    #Return the amount of spin that is not on metal
    
    if os.path.exists(filename) == False:
        return 'Failed'
    
    lines = open(filename,'r').readlines()
    l = [] #indexes of lines that contain 'MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS'
    ltotal = [] #indexes of lines that contain 'Sum of atomic spin populations'
    
    for i in np.arange(len(lines)):
        if 'MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS' in lines[i]:
            l.append(i)
        elif 'Sum of atomic spin populations' in lines[i]:
            ltotal.append(i)
    
    if len(l) == 0 or len(ltotal) == 0:
        return 'Failed'
    
    lastflag = l[-1]
    totalspin = lines[ltotal[0]].split()[-1]
    
    #Finding the next line that contains the name of metal after the line of last flag
    for j in np.arange(lastflag,len(lines)):
        if metal in lines[j]:
            metalspin = lines[j].split()[-1]
            break
            
    spindiff = float(totalspin) - float(metalspin)
    
    return spindiff



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



E_ls,E_is,E_hs,mol2_ls,mol2_is,mol2_hs,dis_ls,dis_is,dis_hs = [],[],[],[],[],[],[],[],[]
deloc_ls,deloc_is,deloc_hs = [],[],[]

indexes = df.index.values
for i in indexes:

    row = df.loc[i]
    refcode,ap1,ap2,apnum,metal,charge = row['refcode'],row['ap1'],row['ap2'],row['apnum'],row['metal'],row['ox_csd']
    aps = ''
    if float(row['apnum']) != 1:
        aps = '_' + str(ap1) + '_' + str(ap2)
    refname = refcode + aps + '_round2f'

    #Record all
    Els,mol2ls,disls = analyze_efei_expanse(refname + '_LS')
    delocls = find_spin_delocalization(refname + '_LS/' + refname + '_LS.out', metal)
    if metal == 'Co' and charge == 2:
        Eis,mol2is,disis,delocis = 'N/A','N/A','N/A',0
    else:
        Eis,mol2is,disis = analyze_efei_expanse(refname + '_IS')
        delocis = find_spin_delocalization(refname + '_IS/' + refname + '_IS.out', metal)
    Ehs,mol2hs,dishs = analyze_efei_expanse(refname + '_HS')
    delochs = find_spin_delocalization(refname + '_HS/' + refname + '_HS.out', metal)
    
    #Resub analysis
    #if row['round1f_Ehs'] == 'Failed':
        #Ehs1,mol2hs1,dishs1 = analyze_efei_expanse(refname1 + '_HS')
        #delochs1 = find_spin_delocalization(refname1 + '_HS/' + refname1 + '_HS.out', metal)
        #df.at[i,'round1f_Ehs'] = Ehs1
        #df.at[i,'round1f_mol2hs'] = mol2hs1
        #if Ehs1 != 'Failed':
            #df.at[i,'round1f_delochs'] = delochs1

    #Record all recording    
    E_ls.append(Els)
    E_is.append(Eis)
    E_hs.append(Ehs)
    mol2_ls.append(mol2ls)
    mol2_is.append(mol2is)
    mol2_hs.append(mol2hs)
    dis_ls.append(disls)
    dis_is.append(disis)
    dis_hs.append(dishs)
    deloc_ls.append(delocls)
    deloc_is.append(delocis)
    deloc_hs.append(delochs)

df['round2f_Els'] = E_ls
df['round2f_Eis'] = E_is
df['round2f_Ehs'] = E_hs
df['round2f_mol2ls'] = mol2_ls
df['round2f_mol2is'] = mol2_is
df['round2f_mol2hs'] = mol2_hs
df['round2f_diss_ls'] = dis_ls
df['round2f_diss_is'] = dis_is
df['round2f_diss_hs'] = dis_hs

df['round2f_delocls'] = deloc_ls
df['round2f_delocis'] = deloc_is
df['round2f_delochs'] = deloc_hs

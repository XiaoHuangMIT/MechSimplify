import numpy as np
import pandas as pd


def round_result(df,idx,n,is_co2=False):
    
    #Return roundn_spin and roundn_diss for complex recorded at df.loc[idx] row
    #df: shall record calculations results as roundn_Els/is/hs, roundn_mol2ls/is/hs, round2_diss_ls/is/hs
    #Tip: failed--neither fully finished in all spin, or has dissociated in at least one spin
    
    roundn = 'round' + str(n)
    
    #If the complex has not been analyzed at all due to failure in previous round
    if df.loc[idx][roundn + '_force'] == 'failed_unperformed':
        return 'failed_unperformed', 'failed_unperformed'
    
    #Check the ground state
    Els, Eis, Ehs = df.loc[idx][roundn + '_Els'], df.loc[idx][roundn + '_Eis'], df.loc[idx][roundn + '_Ehs']
    
    if Els == 'Failed' or Eis == 'Failed' or Ehs == 'Failed':
        spin = 'Failed'
    
    if is_co2 == True:
        if Els < Ehs:
            spin = 'LS'
        else:
            spin = 'HS' 
    else:
        df1 = pd.DataFrame()
        df1['Spin'] = ['LS','IS','HS']
        df1['E'] = [Els,Eis,Ehs]
        df1 = df1.sort_values(by = 'E')
        spin = df1['Spin'].iloc[0]
           
    #Check if the molecule has dissociated in at least one state
    dissls = df.loc[idx][roundn + 'diss_ls']
    dissis = df.loc[idx][roundn + 'diss_is']
    disshs = df.loc[idx][roundn + 'diss_hs']
    if type(dissls) != bool or type(dissis) != bool or type(disshs) != bool:#sanity check
        diss = False
    diss = dissls or dissis or disshs
    
    #If dissociated in at least one spin: label spin not as failed but as diss
    if diss == True and spin == 'Failed': #dont use 'if variable'--prevent future bugs
        spin = 'diss'
        
    return spin, diss    
  
  
  
  def force_second_round(df,idx):
    
    #Based on spin and dissociation status of round 1(1nN), derive whether or not to perform round 2
    #and if to perform, the magnitude of force, with regarding to complex at df.loc[idx]
    #Tip: failed--neither fully finished in all spin, or has dissociated in at least one spin
    
    spin,diss = df.loc[idx]['round1_spin'],df.loc[idx]['round1_diss']
    
    #Case 1: the molecule has failed in round 1 
    if spin == 'Failed' and diss == False: #Double-insurance
        return 'failed_unperformed'
    
    #Case 2: the molecule has dissociated in round 1
    if diss == True:
        return 0.5
    
    #Case 3: the molecule in IS/HS
    if spin =='IS' or spin == 'HS':
        return 0.5
    
    #Case 4: the molecule still in LS
    if spin == 'LS':
        return 1.5
      
      
      
def force_nth_round(df,idx,n,threshold=5):
    
    #Based on force, spin and dissociation status of round n-1 and round n-2, 
    #derive whether or not to perform round n
    #and if to perform, the magnitude of force, with regarding to complex at df.loc[idx]
    #Refer to notes for full logic flowchart
    #Returns: string （conclusion) or float (force magnitude)
    #Consider edge cases: IS/HS in 0.5nN (SCO range 0-0.5nN) or LS in 5nN(threshold: stop screening and conclude
    #as non SCO)
    
    if n == 2:
        return force_second_round(df,idx)
    
    curr,prev = 'round' + str(n-1), 'round' + str(n-2)
    f_curr,spin_curr,diss_curr = df.loc[idx][curr + '_force'],df.loc[idx][curr + '_spin'],df.loc[idx][curr + '_diss']
    f_prev,spin_prev,diss_prev = df.loc[idx][prev + '_force'],df.loc[idx][prev + '_spin'],df.loc[idx][prev + '_diss']
  
    
    #Case 1: the molecule has failed in round n--do not perform round n+1
    if spin_curr == 'Failed':
        return 'failed_unperformed'
    
    f_curr,f_prev = float(f_curr), float(f_prev)
    #Case 2: the molecule has dissociated in round n-2
    if diss_prev == True:
        #If also dissociated in round n-1: further decrease force
        if diss_curr == True:
            result = f_curr - 0.5
        #If in round n-1, the complex has not dissociated but at IS/HS: also decrease force
        elif spin_curr == 'HS' or spin_curr == 'IS':
            result = f_curr - 0.5
        #If non-diss in round n-1, but in LS: the complex cannot do SCO (while structurally intact)
        elif spin_curr == 'LS':
            result = 'concluded nonSCO'
    
    #Case 3: the molecule has not dissociated and in low spin in round n-2:
    elif spin_prev == 'LS':
        #If dissociated in round n-1: cannot SCO
        if diss_curr == True:
            result = 'concluded nonSCO'
        #If still in LS in round n-1: futher increase force
        elif spin_curr == 'LS':
            result = f_curr + 0.5
        #If changed to IS or HS in round n-1: SCO range found!
        elif spin_curr == 'IS' or spin_curr == 'HS':
            result = 'SCO ' + str(f_prev) + ' ' + str(f_curr) #f curr should be 0.5nN higher than f prev
    
    #Case 4: the molecule has not dissociated and in IS/HS in round n-2:
    elif spin_prev == 'IS' or spin_prev == 'HS':
        #If dissociated in round n-1--an error ocurred, since the complex has not dissociated at larger f_prev
        if diss_curr == True:
            result = 'Diss Error'
        #If changed to LS in round n-1: SCO range found!
        elif spin_curr == 'LS':
            result = 'SCO ' + str(f_curr) + ' ' + str(f_prev)
        #IF still in IS or HS in round n-1: further decrease force
        elif spin_curr == 'IS' or spin_curr == 'HS':
            result = f_curr - 0.5
    
    #Consider edge cases:
    if result == 0:
        if spin_curr == 'IS' or spin_curr == 'HS':
            result = 'SCO 0 0.5'
    if type(result) != str and result > threshold:
        return 'concluded-threshold nonSCO'
    
    return result
  
  
  
def performed_in_prev_rounds(df,idx,n):
    
    #Check if force to be analyzed under during round n has been analyzed previously
    #Prev forces: recorded in performed_forces
    #Forces under which all complexes has been analyzed: 0nN (not recorded), 1nN
    #0nN: all LS. Checking performed in force_nth_round
    #1nN: recorded
    #Only check float values of roundn_force (to be performed next)
    
    forces = df.loc[idx]['performed_forces'].split()
    forces = [float(i) for i in fs]
    forces.sort()
    
    fn = df.loc[idx]['round' + str(n) + '_forces']
    fn = float(fn)
    
    if fn in forces:
        return True
    else:
        return False
  
  

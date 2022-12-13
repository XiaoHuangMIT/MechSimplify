import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt


def separate_metal_and_charge(df, count_only = False):
    
    dffe = df[df.metal == 'Fe']
    dffe2 = dffe[dffe.ox_csd == 2]
    dffe3 = dffe[dffe.ox_csd == 3]
    
    dfco = df[df.metal == 'Co']
    dfco2 = dfco[dfco.ox_csd == 2]
    dfco3 = dfco[dfco.ox_csd == 3]
    
    dfmn = df[df.metal == 'Mn']
    dfmn2 = dfmn[dfmn.ox_csd == 2]
    
    lcount = []
    lcount.append(dffe2.shape[0])
    lcount.append(dffe3.shape[0])
    lcount.append(dfco2.shape[0])
    lcount.append(dfco3.shape[0])
    lcount.append(dfmn2.shape[0])
  
    if count_only == True:
        return lcount
    else: 
        return dffe2, dffe3, dfco2, dfco3, dfmn2, lcount

    

def oct_check(df,lsname,isname,hsname):
    
    #Return number of oct and non-oct geometry of hs , (is) and ls: [lso,lsno] , ([iso,isno]), [hso,hsno] 
    #For Co2+: set isname to None
    lso,lsno,iso,isno,hso,hsno = 0,0,0,0,0,0
    
    for i in np.arange(df.shape[0]):
        
        lsmol2 = df.iloc[i][lsname]
        hsmol2 = df.iloc[i][hsname]
        ls3d,hs3d = mol3D(), mol3D()
        ls3d.readfrommol2(lsmol2,readstring=True)
        hs3d.readfrommol2(hsmol2,readstring=True)
        
        if isname != None:
            ismol2 = df.iloc[i][isname]
            is3d = mol3D()
            is3d.readfrommol2(ismol2,readstring = True)
        lslabel = ls3d.IsOct()[0]
        hslabel = hs3d.IsOct()[0]
        if isname != None:
                islabel = is3d.IsOct()[0]
        
        if lslabel == 1:
            lso += 1
            lsno += 0
        elif lslabel == 0:
            lso += 0
            lsno += 1
        if hslabel == 1:
            hso += 1
            hsno += 0
        elif hslabel == 0:
            hso += 0
            hsno += 1
        if isname != None: 
            if islabel == 1:
                iso += 1
                isno += 0
            elif islabel == 0:
                iso += 0
                isno += 1
    
    if isname == None:
        return [lso,lsno],[hso,hsno]
    else:
        return [lso,lsno],[iso,isno],[hso,hsno]

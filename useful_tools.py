import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt


def separate_metal_and_charge(df):
  
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
  
  return dffe2, dffe3, dfco2, dfco3, dfmn2, lcount

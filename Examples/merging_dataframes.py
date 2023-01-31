import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from molSimplify.Classes.mol3D import mol3D

plt.rc('font',size=15)


df = pd.merge(df_5,df_10,how='outer',on='refcode')
df = pd.merge(df,df_12_5,how='outer',on='refcode')
df = pd.merge(df,df_15,how='outer',on='refcode')
df = pd.merge(df,df_17_5,how='outer',on='refcode')
df = pd.merge(df,df_20,how='outer',on='refcode')
df = pd.merge(df,df_20nd,how='outer',on='refcode')
df = pd.merge(df,df_15t,how='outer',on='refcode')
df = pd.merge(df,df_15vvt,how='outer',on='refcode')


df1 = pd.read_csv('round1_e_finished.csv')
df2 = pd.read_csv('round1_scre_finished.csv')
df = pd.concat([df1,df2],axis=0)
#See: https://pandas.pydata.org/docs/reference/api/pandas.concat.html
#Combine DataFrame objects with overlapping columns and return everything. Columns outside the intersection will be filled with NaN values



#Example: Combining collection from round2 and round 3
#All rows in round3 collection are in round2 collection
#Round 3 collection has new columns

df1 = pd.read_csv('round2_all_recorded.csv')

refnames = []
for i in np.arange(df1.shape[0]):
    
    ref = df1.iloc[i]['refcode'] + '_' + str(df1.iloc[i]['ap1']) + '_' + str(df1.iloc[i]['ap2'])
    refnames.append(ref)
df1['refname'] = refnames
df1 = df1.set_index('refname')

columns = df1.columns
drops = []
for column in columns:
    if 'Unnamed' in column:
        drops.append(column)
        
df1 = df1.drop(columns=drops)


df2 = pd.read_csv('round3_all_recorded.csv')

refnames = []
for i in np.arange(df2.shape[0]):
    
    ref = df2.iloc[i]['refcode'] + '_' + str(df2.iloc[i]['ap1']) + '_' + str(df2.iloc[i]['ap2'])
    refnames.append(ref)
df2['refname'] = refnames
df2 = df2.set_index('refname')

columns = df1.columns
drops = []
for column in columns:
    if 'Unnamed' in column:
        drops.append(column)
        
df2 = df2.drop(columns=drops)


df = df1.join(df2[['round3_Els', 'round3_Eis', 'round3_Ehs', 'round3_mol2ls',
       'round3_mol2is', 'round3_mol2hs', 'round3_diss_ls', 'round3_diss_is',
       'round3_diss_hs', 'round3_spin', 'round3_diss', 'round4_force']])
#df['round4_force'].value_counts(dropna=False) #checking

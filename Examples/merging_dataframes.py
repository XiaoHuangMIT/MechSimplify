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

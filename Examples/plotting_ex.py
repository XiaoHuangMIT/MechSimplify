import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rc('font',size=15)



####Bar####
metals = ['Fe2','Fe3','Co2','Co3','Ni2','Mn2'] #NI3+ only 5/3
totals = np.array([181,5,78,25,100,21])
singles = np.array([120,4,57,14,72,14])
dists = np.array([25, 0, 6, 3, 18, 1])
indists = np.array([34, 1, 15, 6, 10, 6])


x = np.arange(len(metals))  # the label locations
width = 0.35

fig, ax = plt.subplots(figsize=(10,8))
rects1 = ax.bar(x - width/2, totals, width, label='Total',color='tab:brown')
rects2 = ax.bar(x + width/2, singles, width, label='Single AP',color='tab:blue')
rects3 = ax.bar(x + width/2, dists, width, bottom=singles,label='Dist nAPs',color='tab:purple')
rects4 = ax.bar(x + width/2, indists, width, bottom=singles+dists,label='Indist nAPs',color='tab:red')

ax.set_ylabel('Counts')
ax.set_title('Mer candidates')
ax.set_xticks(x, metals)
ax.legend(loc='upper center')

ax.bar_label(rects1, label_type='center',padding=3)
ax.bar_label(rects2, label_type='center',padding=3)
ax.bar_label(rects3, label_type='center',padding=3)
ax.bar_label(rects4, label_type='center',padding=3)

fig.tight_layout()

plt.show()


################################Analayze Spin Splitting#####################################
def plot_Ehsls(Ecolumn, title, size=(8,6), width=0.3, color='tab:blue', return_counts=True):
    
    counts = np.histogram(Ecolumn,bins = [0,5,10,20,30,1000])
    counts = counts[0]
    
    Es = ['0-5','5-10','10-20','20-30','30+']
    x = np.arange(len(Es))
    wid = width
    
    fig, ax = plt.subplots(figsize=size)
    rects1 = ax.bar(x, counts, wid,color=color)
    ax.bar_label(rects1, label_type='center',padding=3)
    
    ax.set_ylabel('Counts')
    ax.set_title(title)
    ax.set_xticks(x, Es)
    ax.set_xlabel('HS-LS Spin Splitting Energy (kcal/mol)')
    
    plt.show()


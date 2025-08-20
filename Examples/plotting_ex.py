import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rc('font',size=15)



def figure_formatting(): # make the plot look xmgrace-esque
    #font = {'family': 'sans-serif', 'weight': 'bold', 'size': 15}
    plt.rcParams['font.sans-serif'] = ['Helvetica']
    #plt.rc('font', **font)
    plt.rcParams['axes.linewidth'] = 2.0
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.minor.size'] = 5
    plt.rcParams['xtick.major.width'] = 2.0
    plt.rcParams['xtick.minor.width'] = 2.0
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.minor.size'] = 5
    plt.rcParams['ytick.major.width'] = 2.0
    plt.rcParams['ytick.minor.width'] = 2.0
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['mathtext.default'] = 'regular'
    plt.rcParams['legend.fancybox'] = False # no rounded legend box
    plt.rcParams['svg.fonttype']='none'
figure_formatting()
    


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



####Histogram####
figure_formatting()
plt.plot(figsize=(8,6))
plt.hist([df1['E_gap'],df2['E_gap'],df3['E_gap']],
           bins=np.arange(0,121,5),
             stacked=True,edgecolor='black',
           color=[(0,0,1,1),(0,0.5,0,1),(1,0.65,0,1)],linewidth=2,label=['0-0.5','0.5-1','1-1.5'])

plt.vlines(70,0,41,linestyles='dashed',linewidth=2,color='darkgreen')
plt.vlines(80,0,41,linestyles='dashed',linewidth=2,color='darkorange')

plt.xticks(np.arange(0,121,20),np.arange(0,121,20))
plt.xlabel('Partial Dissociation Energy (kcal/mol)')

plt.yticks(np.arange(0,41,5))
plt.ylabel('counts')
plt.ylim(0,40)

plt.title('Ortho-Type Candidates',size=16)
leg = plt.legend(prop = {'size' : 13,'family': 'Helvetica'},edgecolor='black')
leg.get_frame().set_linewidth(2)
    plt.show()


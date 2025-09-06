# -*- coding: utf-8 -*-
"""
Created on Tue Aug 26 17:31:19 2025

@author: xinwe
"""

import matplotlib.pyplot as plt
import numpy as np
import os

figfontsize=20

filedir ='./Analysis/' 
date = '2025_08_28_'
loaddir = filedir+'analyzeddata/'
savedir = filedir+'figures/'

#load all data
jGCaMP8s = ['2025_08_20_jGCaMP8s_F_800res_every5sec_1stlymphgland',
            '2025_08_20_jGCaMP8s_F_800res_every5sec_3rdlymphgland',
            '2025_08_20_jGCaMP8s_F_800res_every5sec_5thlymphgland',
            '2025_08_23_jGCaMP8s_F_800res_every5sec_1stlymphgland',
            '2025_08_23_jGCaMP8s_F_800res_every5sec_2ndlymphgland',
            '2025_08_25_jGCaMP8s_F_800res_every5sec_1stlymphgland',
            '2025_08_25_jGCaMP8s_F_800res_every5sec_3rdlymphgland',
            '2025_08_25_jGCaMP8s_F_800res_every5sec_4thlymphgland']

GFP = ['2025_08_21_GFP_F_800res_every5sec_2ndlymphgland', 
       '2025_08_21_GFP_F_800res_every5sec_3rdlymphgland',  
       '2025_08_21_GFP_F_800res_every5sec_4thlymphgland', 
       '2025_08_25_GFP_M_800res_every5sec_2ndlymphgland']

#load duration of spikes for all
#load signal intensities for jGCaMP8s (and all intensities for cell count)
#load all intensities for GFP

jGCaMP8s_duration = [ [] for i in range(len(jGCaMP8s))]
jGCaMP8s_intensities = [ [] for i in range(len(jGCaMP8s))]

for i in range(len(jGCaMP8s)):
    jGCaMP8s_duration[i]=np.load(loaddir + date + jGCaMP8s[i] +'_signalduration.npy')
    jGCaMP8s_intensities[i]=np.load(loaddir + date + jGCaMP8s[i] +'_basalintensitypercell.npy')

GFP_duration = [ [] for i in range(len(GFP))]
GFP_intensities = [ [] for i in range(len(GFP))]

for i in range(len(GFP)):
    GFP_duration[i]=np.load(loaddir + date + GFP[i] +'_signalduration.npy')
    GFP_intensities[i]=np.load(loaddir + date + GFP[i] +'_basalintensitypercell.npy')


#concatenate because we don't need per-cell resolution, remove zero values that can stem from artifacts in looping 
jGCaMP8s_duration=np.concatenate(jGCaMP8s_duration)
jGCaMP8s_duration=jGCaMP8s_duration[jGCaMP8s_duration !=0]

GFP_duration=np.concatenate(GFP_duration)
GFP_duration=GFP_duration[GFP_duration !=0]

#%%
jGCaMP8sSigCount=len(jGCaMP8s_duration)
jGCaMP8sCellCount=len(np.concatenate(jGCaMP8s_intensities))

GFPSigCount=len(GFP_duration)
GFPCellCount=len(np.concatenate(GFP_intensities))

#%%durations as histograms

#important for range to start at 1, the algorithm will return events of duration '0' because of code artefacts, this excludes them

fig=plt.figure(figsize=[16,8],dpi=600)

subplot1=fig.add_subplot(121)
subplot1.hist(jGCaMP8s_duration, bins = 25, range=(1,25), align='left', color='teal', rwidth=0.95) #extreme outliers excluded from plot
plt.ylim(0,375)
plt.yticks(ticks=[0,50,100,150,200,250,300,350], labels =['0','50','100','150','200','250','300','350'],fontsize=figfontsize)
plt.xticks(ticks=[3,6,9,12,18,24], labels =['15','30','45','60','90','120'],fontsize=figfontsize)
plt.xlabel('Signal duration (sec)',fontsize=figfontsize)
plt.title('Detected signals (jGCaMP8s)\n'+str(jGCaMP8sSigCount)+' events/'+str(jGCaMP8sCellCount)+' cells',fontsize=figfontsize+5)

subplot2=fig.add_subplot(122)
subplot2.hist(GFP_duration, bins = 25, range=(1,25), align='left', color='darkgreen', rwidth=0.95) 
plt.ylim(0,375)
plt.yticks(ticks=[0,50,100,150,200,250,300,350], labels =['0','50','100','150','200','250','300','350'],fontsize=figfontsize)
plt.xticks(ticks=[3,6,9,12,18,24], labels =['15','30','45','60','90','120'],fontsize=figfontsize)
plt.xlabel('Signal duration (sec)',fontsize=figfontsize)
plt.title('Detected signals (EGFP)\n'+str(GFPSigCount)+' events/'+str(GFPCellCount)+' cells',fontsize=figfontsize+5)
plt.show()

fig.savefig(os.path.join(savedir, date+'detectedsignals'))

#%%compare variation in all intensities (keep it simple as opposed to comparing signal intensities)

jGCaMP8sintensitydist = np.log(np.concatenate(jGCaMP8s_intensities))
GFPintensitydist = np.log(np.concatenate(GFP_intensities))

#spread on intensities
stdjGCaMP8s = np.std(jGCaMP8sintensitydist)
stdGFP = np.std(GFPintensitydist)

#function for making beeswarmplot
def simple_beeswarm2(y, nbins=None, width=1.):
    """
    Returns x coordinates for the points in ``y``, so that plotting ``x`` and
    ``y`` results in a bee swarm plot.
    Function created by The-Python-Graph-Gallery
    """
    # Convert y to a numpy array to ensure it is compatible with numpy functions
    y = np.asarray(y)
    # If nbins is not provided, calculate a suitable number of bins based on data length
    if nbins is None:
        # nbins = len(y) // 6
        nbins = np.ceil(len(y) / 6).astype(int)
    # Get the histogram of y and the corresponding bin edges
    nn, ybins = np.histogram(y, bins=nbins)
    # Find the maximum count in any bin to be used in calculating the x positions
    nmax = nn.max()
    # Create an array of zeros with the same length as y, to store x-coordinates
    x = np.zeros(len(y))
    # Divide indices of y-values into corresponding bins
    ibs = []
    for ymin, ymax in zip(ybins[:-1], ybins[1:]):
        # Find the indices where y falls within the current bin
        i = np.nonzero((y > ymin) * (y <= ymax))[0]
        ibs.append(i)
    # Assign x-coordinates to the points in each bin
    dx = width / (nmax // 2)
    for i in ibs:
        yy = y[i]
        if len(i) > 1:
            # Determine the starting index (j) based on the number of elements in the bin
            j = len(i) % 2
            # Sort the indices based on their corresponding y-values
            i = i[np.argsort(yy)]
            # Separate the indices into two halves (a and b) for arranging the points
            a = i[j::2]
            b = i[j+1::2]
            # Assign x-coordinates to points in each half of the bin
            x[a] = (0.5 + j / 3 + np.arange(len(b))) * dx
            x[b] = (0.5 + j / 3 + np.arange(len(b))) * -dx
    return x

#plot and compare jGCaMP8s_signalintensities and GFP_intensities[i]
boxprops = dict(linestyle='-', linewidth=4, color='black')
whiskerprops = dict(linestyle='-', linewidth=4, color='black')
capprops = dict(linestyle='-', linewidth=4, color='black')
medianprops = dict(linestyle='-', linewidth=4, color='firebrick')

boxplot_data = [jGCaMP8sintensitydist,GFPintensitydist]
fig1,ax1=plt.subplots(1, 1, figsize=(8, 8),dpi=600)

x = simple_beeswarm2(jGCaMP8sintensitydist, width=0.25)
y = simple_beeswarm2(GFPintensitydist, width=0.25)
ax1.plot(x+1, jGCaMP8sintensitydist, '.', color='teal',markersize='15')
ax1.plot(y+2, GFPintensitydist, '.', color='darkgreen',markersize='15')

ax1.boxplot(boxplot_data, boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, capprops=capprops, widths=0.5, showfliers=False)
ax1.set_xticks([1,2])
ax1.set_xticklabels(['jGCaMP8s\n ($\sigma$ ='+str(np.round(stdjGCaMP8s,decimals=3))+')','EGFP\n ($\sigma$ ='+str(np.round(stdGFP,decimals=3))+')'],fontsize=figfontsize) #convert to minutes
ax1.set_ylim([-1,8.5])
ax1.set_yticks([0,2,4,6,8])
ax1.set_yticklabels(['0','2','4','6','8'],fontsize=figfontsize) 
ax1.set_ylabel('log(fluorescence intensity)',fontsize=figfontsize)
ax1.set_title('Basal cell intensities',fontsize=figfontsize+5)
plt.show()

fig1.savefig(os.path.join(savedir, date+'distributionofsignals'))

#%%useful stats

#median duration
medjGCAMP8s=np.median(jGCaMP8s_duration)
medGFP=np.median(GFP_duration)

#counts 
for i in range (int(max(jGCaMP8s_duration))+1):
    frames=(jGCaMP8s_duration > i).sum()
    print('The number of jGCaMP8s signals with duration greater than '+str(i)+' frames is '+str(frames))
    
for i in range (int(max(GFP_duration))+1):
    frames=(GFP_duration > i).sum()
    print('The number of GFP signals with duration greater than '+str(i)+' frames is '+str(frames))
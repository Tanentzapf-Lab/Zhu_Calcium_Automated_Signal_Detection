# -*- coding: utf-8 -*-
"""
This in-development script file is for processing individual jGCaMP8s and control GFP videos 
Returns traces for every cell and also basal fluorescence values of every cell
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import math
import peakutils

from skimage import io

date = '2025_08_28_'
filedir ='./Analysis/' 
#when code is running ideally date, filedir and filename are the only lines that need to be changed at the beginning
#filename = '2025_08_20_jGCaMP8s_F_800res_every5sec_1stlymphgland' 
#filename = '2025_08_20_jGCaMP8s_F_800res_every5sec_3rdlymphgland' 
#filename = '2025_08_20_jGCaMP8s_F_800res_every5sec_5thlymphgland' 
#filename = '2025_08_23_jGCaMP8s_F_800res_every5sec_1stlymphgland' 
#filename = '2025_08_23_jGCaMP8s_F_800res_every5sec_2ndlymphgland'
#filename = '2025_08_25_jGCaMP8s_F_800res_every5sec_1stlymphgland'
#filename = '2025_08_25_jGCaMP8s_F_800res_every5sec_3rdlymphgland'
#filename = '2025_08_25_jGCaMP8s_F_800res_every5sec_4thlymphgland'

#filename = '2025_08_21_GFP_F_800res_every5sec_2ndlymphgland' 
#filename = '2025_08_21_GFP_F_800res_every5sec_3rdlymphgland'  
#filename = '2025_08_21_GFP_F_800res_every5sec_4thlymphgland' 
filename = '2025_08_25_GFP_M_800res_every5sec_2ndlymphgland'

backgroundfluor = 72.9976080108125 #this value for background fluorescence was determined experimentally as the intensity value 2 standard deviations above the mean for a sample without a green fluorophore

savedir = filedir+'analyzeddata/'

#important for the mask to be a tif, when imreading the equivalent png for some reason it is not read in as integer values (rescaling occurs)
avrgmask = io.imread(filedir+'Mask/LblImg_'+filename+'_redphotobleachexpcorrected_AVG.tif')
    #cell arranged as [t,y,x]
green = io.imread(filedir+'CalciumSignal/'+filename+'_greenphotobleachexpcorrected.tif')


#%%first need to make a mask where avrgmask is just repeated n times (n being the number of frames that were averaged together for tracking)

navrg=5 #5 works well but possible to change this number

nframes=avrgmask.shape[0]

#initialize a mask
mask = np.empty([nframes*navrg,avrgmask.shape[1],avrgmask.shape[2]]) #shape is (time, Ypos, Xpos)

for i in range(nframes):
    for j in range(navrg):
        mask[(i*navrg+j),:,:] = avrgmask[i,:,:]

#%%get intensities for every cell

#initialize meanintensities, first axis track index, second axis time
meanintensities = np.empty([int(np.max(mask)),mask.shape[0]])

#work with each slice individually (i is slice number)
for i in range(mask.shape[0]):
    intensities = [[] for j in range(int(np.max(mask)))]

#mask.shape[0] is time, mask.shape[1] is the vertical direction from top-bottom, mask.shape[2] is horizontal from left to right

#first take mean across a single slice
    for j in range(mask.shape[1]):
     for k in range(mask.shape[2]):
         if mask[i,j,k] > 0: #value of -1 corresponds to cells not in tracks, avoid
             intensities[int(mask[i,j,k]-1)].append(green[i,j,k]) #index is one-off from cell number

#   meanintensities = [[] for j in range(np.max(mask))]
    for j in range(int(np.max(mask))):
        meanintensities[j,i]=np.mean(intensities[j])

#%%okay now need to remove any traces that have nan in them, or rather move traces that are complete to a new array

a=0 #initialize a new counter for where to put traces

containsnan = np.zeros([meanintensities.shape[0]])

for i in range(meanintensities.shape[0]):
    #go through each trace one by one
    for j in range(meanintensities.shape[1]):
        if math.isnan(meanintensities[i,j]) is True:
            containsnan[i-1]=1
            break
        
fulltraces = np.empty([int(len(containsnan)-sum(containsnan)),meanintensities.shape[1]]) #don't have to initialize fulltraces any longer than necessary
for i in range(meanintensities.shape[0]):        
    if containsnan[i-1]==0:
        fulltraces[a,:]=meanintensities[i,:]
        a=a+1

#%%background subtraction - any numbers that dip below 0 after the subtraction should be set to 0

fulltracesbkgdc=np.array(fulltraces)-backgroundfluor 

for i in range (fulltraces.shape[0]):
    for j in range (fulltraces.shape[1]):
        if fulltracesbkgdc[i,j]<0:
            fulltracesbkgdc[i,j]=0

#%%calcium activity in all complete background-corrected traces

figfontsize=25

fig, ax = plt.subplots(figsize=(15,12),dpi=600)
#ax.imshow(fulltraces, cmap='inferno', aspect='auto', interpolation='nearest')
ax.pcolormesh(fulltracesbkgdc, cmap='inferno')
ax.set_xticks([0,36,72,108,144,180])
ax.set_xticklabels(['0','3','6','9','12','15'],fontsize=figfontsize) #convert to minutes
ax.set_yticks([10,20,30,40,50,60,70])
ax.set_yticklabels(['10','20','30','40','50','60','70'],fontsize=figfontsize) 
ax.set_xlabel('Minutes',fontsize=figfontsize+5)
ax.set_ylabel('Cell index',fontsize=figfontsize+5)
#ax.set_title('Calcium activity in progenitor cells',fontsize=figfontsize+10)
plt.show()

fig.savefig(os.path.join(savedir, date + filename+'_alltraces'))

#%%now that traces have been background-corrected determine baseline per cell? alternative way of determining calcium signaling events: count every instance of the signal getting above 25% baseline level as a signaling event, calcium transient
#also populate fullsignal array to see all points at which a signal is happening

frames=np.arange(1,mask.shape[0]+1)

fullsignal = np.zeros([fulltracesbkgdc.shape[0],fulltracesbkgdc.shape[1]]) 
for i in range(fulltracesbkgdc.shape[0]):        
    extrace=fulltracesbkgdc[i,:] 
    baseline_values = peakutils.baseline(extrace)

    for j in range(len(extrace)):
        if extrace[j]>baseline_values[j]+np.mean(baseline_values)*0.25+20:#how much percent above mean(baseline) to count as a signaling event plus a constant to account for stable noise
            fullsignal[i,j]=1

    #plt.plot(frames, extrace) #plotting the trace
    #plt.plot(frames, baseline_values) 

    #for j in range(len(extrace)):
        #if fullsignal[i,j]==1:
            #plt.scatter(frames[j],extrace[j],color='r',marker='x')

    #plt.title('Trace '+str(i))
    #plt.ylim((np.mean(baseline_values)/2, np.mean(baseline_values)*4)) 
    #plt.xticks(ticks=[0,36,72,108,144,180],labels=['0','3','6','9','12','15'])
    #plt.xlabel('Minutes')
    #plt.ylabel('Fluorescence Intensity')
    #plt.show()

# =============================================================================
# #%%plot just one specific trace     
# extrace=fulltracesbkgdc[78,:] #select trace to plot here
# baseline_values = peakutils.baseline(extrace)
# 
# fig, ax = plt.subplots(figsize=(15,12),dpi=600)
# ax.plot(frames, extrace, linewidth=5) #plotting the trace
# ax.plot(frames, baseline_values, linewidth=5) 
# 
# for j in range(len(extrace)):
#     if fullsignal[78,j]==1:
#         plt.scatter(frames[j],extrace[j],color='r',marker='x',s=250,zorder=2,linewidths=5)
# 
# #ax.title('Example trace')
# ax.set_ylim((np.mean(baseline_values)/2, np.mean(baseline_values)*4)) 
# ax.set_xticks([0,36,72,108,144,180])
# ax.set_xticklabels(['0','3','6','9','12','15'],fontsize=figfontsize)
# ax.set_yticks([])
# ax.set_yticklabels([],fontsize=figfontsize)
# ax.set_xlabel('Minutes',fontsize=figfontsize+5)
# ax.set_ylabel('Fluorescence Intensity',fontsize=figfontsize+5)
# plt.show()
# 
# fig.savefig(os.path.join(savedir, date + filename+'_exampletrace'))
# =============================================================================

#%%plot fullsignal

fig, ax = plt.subplots(figsize=(15,12),dpi=600)
#ax.imshow(fulltraces, cmap='inferno', aspect='auto', interpolation='nearest')
ax.pcolormesh(fullsignal, cmap='inferno')
ax.set_xticks([0,36,72,108,144,180])
ax.set_xticklabels(['0','3','6','9','12','15'],fontsize=figfontsize) #convert to minutes
ax.set_yticks([10,20,30,40,50,60,70])
ax.set_yticklabels(['10','20','30','40','50','60','70'],fontsize=figfontsize) 
ax.set_xlabel('Minutes',fontsize=figfontsize+5)
ax.set_ylabel('Cell index',fontsize=figfontsize+5)
#ax.set_title('Detected signaling events in progenitor cells',fontsize=figfontsize+10)
plt.show()

fig.savefig(os.path.join(savedir, date + filename+'_detectedsignals'))

#%%calculate and save fluorescence intensity in signaling portions

intensitybasal = [[] for i in range(fullsignal.shape[0])] #for each cell, make list of intensities at spiking positions
spikelengths = [[] for i in range(fullsignal.shape[0])] #make a list of spike lengths for each cell

for i in range(fullsignal.shape[0]):
    a=0 #initialize counter at the start for each cell (do not carry over incomplete signal from previous cell)
    for j in range(fullsignal.shape[1]):
        if fullsignal[i,j]==1: #signal start or ongoing
            a=a+1 #signal length has increased
        if fullsignal[i,j]==0: #not a signal
            intensitybasal[i].append(fulltracesbkgdc[i,j])
        if fullsignal[i,j]==0 and fullsignal[i,j-1]==1: #signal end
            spikelengths[i].append(a)
            a=0 #reset signal length
                
#flatten lists for further processing
flat_spikelengths = list(np.concatenate(spikelengths).flat) 
flat_spikelengths = np.array(flat_spikelengths) #this can be used to make histogram of signal duration

#for intensityinsignal, take mean and stdev per cell
mean_basalintensity = np.zeros(fullsignal.shape[0])
for i in range(fullsignal.shape[0]):
    mean_basalintensity[i] = np.mean(intensitybasal[i])

basalintensitypercell = mean_basalintensity[~np.isnan(mean_basalintensity)] #this can be used to calculate CV of jGCaMP8s samples (is mean intensity of signal in signaling cells only)

np.save(savedir + date + filename+ '_signalduration.npy', flat_spikelengths)
np.save(savedir + date + filename+ '_basalintensitypercell.npy', basalintensitypercell)

#these arrays will be loaded in separately in downstream scripts for figure creation




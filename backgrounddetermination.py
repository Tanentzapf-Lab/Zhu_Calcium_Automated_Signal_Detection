# -*- coding: utf-8 -*-
"""
This in-development script file is for processing a background video 
"""

import matplotlib.pyplot as plt
import numpy as np
import math
from skimage import io

filedir ='./Analysis/' 
filename = '2025_08_22_mCherryonly' #name of file containing data for background correction

savedir = filedir+'pythonoutput/'

#important for the mask to be a tif, when imreading the equivalent png for some reason it is not read in as integer values (rescaled)
avrgmask = io.imread(filedir+'Mask/LblImg_'+filename+'_redphotobleachexpcorrected_AVG.tif')
    #cell arranged as [t,y,x]
green = io.imread(filedir+'CalciumSignal/'+filename+'_greenphotobleachexpcorrected.tif')


#%%first need to make a mask where avrgmask is just repeated n times (n being the number of frames that were averaged together for tracking)

navrg=5 

nframes=avrgmask.shape[0]

#initialize a mask
mask = np.empty([nframes*navrg,avrgmask.shape[1],avrgmask.shape[2]]) #shape is (time, Ypos, Xpos)

for i in range(nframes):
    for j in range(navrg):
        mask[(i*navrg+j),:,:] = avrgmask[i,:,:]

#%%version for multislice, get intensities for every cell

#initialize meanintensities, first axis track index, second axis time

meanintensities = np.empty([int(np.max(mask)),mask.shape[0]])

#work with each slice individually (i is slice number)
for i in range(mask.shape[0]):

    intensities = [[] for j in range(int(np.max(mask)))]

#mask.shape[0] is time, mask.shape[1] is the vertical direction from top-bottom, mask.shape[2] is horizontal from left to right

#first take mean across a single slice
    for j in range(mask.shape[1]):
     for k in range(mask.shape[2]):
         if mask[i,j,k] != 0:
             intensities[int(mask[i,j,k]-1)].append(green[i,j,k]) #index is one-off from cell number

    for j in range(int(np.max(mask))):
        meanintensities[j,i]=np.mean(intensities[j])
        

#%%remove any traces that have nan in them, or rather move traces that are complete to a new array
#also move the corresponding lines from the same indices in meanX_pos and meanY_pos

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
        
#%%plotting some example traces - this is calcium activity in all complete traces
fig, ax = plt.subplots()
#ax.imshow(fulltraces, cmap='inferno', aspect='auto', interpolation='nearest')
ax.pcolormesh(fulltraces, cmap='inferno')
ax.set_xticks([0,36,72,108,144,180])
ax.set_xticklabels(['0','3','6','9','12','15']) #convert to minutes
ax.set_xlabel('Minutes')
ax.set_ylabel('Cell index')
ax.set_title('Fluorescence in green channel')
plt.show()

#%%

backgroundtosubtract=np.mean(fulltraces)+2*np.std(fulltraces)
print('The intensity value 2 standard deviations above the mean for background fluoresence is '+str(backgroundtosubtract)+'. Subtract this value when background correcting experimental samples.')

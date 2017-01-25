'''
        Script to calculte the change in synaptic strength as seen in Fig. 4 of:
        
        Graupner M and Brunel N (2012). 
        Calcium-based plasticity model explains sensitivity of synaptic changes to spike pattern, rate, and dendritic location. 
        PNAS 109 (10): 3991-3996.
        
        The stimulation protocol consits of pre-post spike-pairs with varying time difference deltaT 
        presented at different frequencies ranging from 1 to 50 Hz. 
        
'''

from scipy import *
from numpy import *
from pylab import *
import os
import sys
import time
import math
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from timeAboveThreshold.timeAboveThreshold import *
from synapticChange import synapticChange

##########################################################
# output directory
outputDir = 'outputFigures/'

##########################################################
# Parameter of the stimulation protocol
frequency   = ([1.,20.,30.,40.,50.])   # frequency of spike-pair presentations in Hz
N           = 75.     # number of spike-pair presentations
DeltaTstart = -0.1    # start time difference between pre- and post-spike, in sec
DeltaTend   =  0.1    # end time difference between pre- and post-spike, in sec
DeltaTsteps =  2000  # steps between start and end value

# chose the 'cortical slices' parameter set 
synChange = synapticChange('cortical slices')

# initialize class which calculates the time the calcium trace spends above threshold
tat = timeAboveThreshold(synChange.tauCa, synChange.Cpre, synChange.Cpost, synChange.thetaD, synChange.thetaP)

nCases = len(frequency)

# initialize empty arrays 
sChange = zeros((nCases,DeltaTsteps+1))
meanUP    = zeros((nCases,DeltaTsteps+1))
meanDOWN  = zeros((nCases,DeltaTsteps+1))
alphaD     = zeros((nCases,DeltaTsteps+1))
alphaP     = zeros((nCases,DeltaTsteps+1))
deltaT    = zeros((nCases,DeltaTsteps+1))

# loop over different parameter-sets
for k in range(nCases):
        interval    = 1./frequency[k]
        if (interval/2. < fabs(DeltaTstart)) or (interval/2. < DeltaTend):
                DeltaTstart = - interval/2.
                DeltaTend   = interval/2.
        deltaT[k] = linspace(DeltaTstart,DeltaTend,DeltaTsteps+1)
        n = 0
        # loop over range of deltaT values
        for i in deltaT[k]:
                # fraction of time spent above threshold
                (alphaD[k,n],alphaP[k,n]) = tat.spikePairFrequency(i-synChange.D,frequency[k])
                # calculate all values for the change in synaptic strength
                synChange.changeInSynapticStrength(N*interval,0.5,alphaD[k,n],alphaP[k,n])
                
                # mean value of the synaptic strength right at the end of the stimulation protocol
                meanUP[k,n]   =  synChange.meanUP
                meanDOWN[k,n] =  synChange.meanDOWN
                
                # change in synaptic strength after/before
                sChange[k,n] = synChange.synChange
                # 
                n+=1
        

##########################################################
# generate figure 

fig_width = 5 # width in inches
fig_height = 4  # height in inches
fig_size =  [fig_width,fig_height]
params = {'axes.labelsize': 14,
          'axes.titlesize': 13,
          'font.size': 11,
          'xtick.labelsize': 11,
          'ytick.labelsize': 11,
          'figure.figsize': fig_size,
          'savefig.dpi' : 600,
          'axes.linewidth' : 1.3,
          'ytick.major.size' : 4,      # major tick size in points
          'xtick.major.size' : 4      # major tick size in points
          #'edgecolor' : None
          #'xtick.major.size' : 2,
          #'ytick.major.size' : 2,
          }
rcParams.update(params)

# create figure instance
fig = plt.figure()

# define sub-panel grid and possibly width and height ratios
gs = gridspec.GridSpec(1, 1)
plt.subplots_adjust(left=0.14, right=0.92, top=0.92, bottom=0.18)

ax0 = plt.subplot(gs[0])

# title
ax0.set_title('cortical slices')

# diplay of data
ax0.axhline(y=1,ls='--',color='0.5',lw=2)
ax0.axvline(x=0,ls='--',color='0.5',lw=2)
for k in range(len(frequency)):
        ax0.plot(deltaT[k],sChange[k],lw=2,label=str(frequency[k]))

# removes upper and right axes 
# and moves left and bottom axes away
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.spines['bottom'].set_position(('outward', 10))
ax0.spines['left'].set_position(('outward', 10))
ax0.yaxis.set_ticks_position('left')
ax0.xaxis.set_ticks_position('bottom')

ax0.set_xlim(-0.06,0.06)

ax0.legend(frameon=False)

ax0.xaxis.set_major_locator(MaxNLocator(5))
#ax0.yaxis.set_major_locator(MaxNLocator(6))

plt.ylabel('change in synaptic strength')
plt.xlabel(r'$\Delta t$ (sec)')

fname = os.path.basename(__file__)

savefig(outputDir+fname[:-2]+'png')
savefig(outputDir+fname[:-2]+'pdf')



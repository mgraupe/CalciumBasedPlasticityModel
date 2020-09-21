'''
        Script to calculte the change in synaptic strength as seen in Fig. 3 of:
        
        Graupner M and Brunel N (2012). 
        Calcium-based plasticity model explains sensitivity of synaptic changes to spike pattern, rate, and dendritic location. 
        PNAS 109 (10): 3991-3996.
        
        The stimulation protocol consits of pre-post spike-pairs and pre-spike post-pairs with varying time 
        difference deltaT presented at 5 Hz. The number of pre-spike post-pair stimuli is furthermore 
        varied.
        
'''

from scipy import *
from numpy import *
from pylab import *
import os
import sys
import time
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.append("..")
from timeAboveThreshold.timeAboveThreshold import *
from synapticChange import synapticChange


##########################################################
# output directory
outputDir = 'outputFigures/'

##########################################################
# Parameter of the stimulation protocol
frequency   = 5.    # frequency of spike-pair presentations in Hz
interval    = 1./frequency
N           = [200,100,30]     # number of spike-pair presentations
DeltaTstart = -0.1    # start time difference between pre- and post-spike, in sec
DeltaTend   =  0.1    # end time difference between pre- and post-spike, in sec
DeltaTsteps =  2000  # steps between start and end value
deltaBurst  = 0.0115

# chose parameter set
plasticityCase = 'hippocampal slices'
synChange = synapticChange('genericCase',plasticityCase)

# initialize class which calculates the time the calcium trace spends above threshold
tat = timeAboveThreshold(synChange.tauCa, synChange.Cpre, synChange.Cpost, synChange.thetaD, synChange.thetaP)

nCases = 3

# array of the paramter to vary 
if (interval/2. < fabs(DeltaTstart)) or (interval/2. < DeltaTend):
        DeltaTstart = - interval/2.
        DeltaTend   = interval/2.
        
deltaT = linspace(DeltaTstart,DeltaTend,DeltaTsteps+1)

# initialize empty arrays 
sChange = zeros((nCases,DeltaTsteps+1))
meanUP    = zeros((nCases,DeltaTsteps+1))
meanDOWN  = zeros((nCases,DeltaTsteps+1))
alphaD     = zeros((nCases,DeltaTsteps+1))
alphaP     = zeros((nCases,DeltaTsteps+1))

# loop over different parameter-sets
for k in range(nCases):
        
        #
        n = 0
        # loop over range of deltaT values
        for i in deltaT:
                # fraction of time spent above threshold
                # spike-pair
                if k == 0:
                        (alphaD[k,n],alphaP[k,n]) = tat.spikePairFrequency(i-synChange.D,frequency)
                # pre-spike and post-pair
                else:
                        (alphaD[k,n],alphaP[k,n]) = tat.preSpikePostPair(i-synChange.D,frequency,deltaBurst)
                
                # calculate all values for the change in synaptic strength
                synChange.changeInSynapticStrength(float(N[k])*interval,0.5,alphaD[k,n],alphaP[k,n])
                
                # mean value of the synaptic strength right at the end of the stimulation protocol
                meanUP[k,n]   =  synChange.meanUP
                meanDOWN[k,n] =  synChange.meanDOWN
                
                # change in synaptic strength after/before
                sChange[k,n] = synChange.synChange
                # 
                n+=1
        

##########################################################
# generate figure 

fig_width = 6 # width in inches
fig_height = 12  # height in inches
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
gs = gridspec.GridSpec(3, 1)

# define vertical and horizontal spacing between panels
gs.update(wspace=0.3,hspace=0.4)

# possibly change outer margins of the figure
plt.subplots_adjust(left=0.15, right=0.92, top=0.92, bottom=0.1)

for k in range(nCases):
        ax0 = plt.subplot(gs[k])

        # title
        #ax0.set_title(plasticityCase)

        # diplay of data
        ax0.axhline(y=1,ls='--',color='0.5',lw=2)
        ax0.axvline(x=0,ls='--',color='0.5',lw=2)
        ax0.plot(deltaT,sChange[k],lw=3,c='m')

        # removes upper and right axes 
        # and moves left and bottom axes away
        ax0.spines['top'].set_visible(False)
        ax0.spines['right'].set_visible(False)
        ax0.spines['bottom'].set_position(('outward', 10))
        ax0.spines['left'].set_position(('outward', 10))
        ax0.yaxis.set_ticks_position('left')
        ax0.xaxis.set_ticks_position('bottom')
        
        ax0.set_xlim(DeltaTstart,DeltaTend)

        ax0.xaxis.set_major_locator(MaxNLocator(5))
        #ax0.yaxis.set_major_locator(MaxNLocator(6))
        
        plt.ylabel('change in synaptic strength')
        if k==2:
                plt.xlabel(r'$\Delta t$ (sec)')

fname = os.path.basename(__file__)

savefig(outputDir+fname[:-2]+'png')
savefig(outputDir+fname[:-2]+'pdf')


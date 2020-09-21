'''
        Script to calculte the change in synaptic strength as seen in Fig. 2 of:
        
        Graupner M and Brunel N (2012). 
        Calcium-based plasticity model explains sensitivity of synaptic changes to spike pattern, rate, and dendritic location. 
        PNAS 109 (10): 3991-3996.
        
        The stimulation protocol consits of pre-post spike-pairs with varying time difference deltaT 
        presented at 1 Hz. 
        The pre- and postsynaptically induced calcium amplitudes and the depression/potentiation
        thresholds are varied to illustrate the diversity of STDP curves emerging from the model. 
        
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

sys.path.append("..") # adds parent directory to python path
from timeAboveThreshold.timeAboveThreshold import *
from synapticChange import synapticChange

##########################################################
# output directory
outputDir = 'outputFigures/'
numSimDir = 'numericalSimulation/'
##########################################################
# Parameter of the stimulation protocol
frequency   = 1.    # frequency of spike-pair presentations in Hz
interval    = 1./frequency
N           = 60.     # number of spike-pair presentations
DeltaTstart = -0.1    # start time difference between pre- and post-spike, in sec
DeltaTend   =  0.1    # end time difference between pre- and post-spike, in sec
DeltaTsteps =  2000  # steps between start and end value

# chose from one of the six different cases : DP, DPD, DPDprime, P, D, Dprime
plasticityCases = ['DP','P','DPD','D','DPDprime','Dprime']
nCases = len(plasticityCases)

#dpSim = np.loadtxt(numSimDir+plasticityCases[0]+'_curve/final_camkII_state_old.dat')
simData = {}
for p in plasticityCases:
    simData.update({p: np.loadtxt(numSimDir+'output/%s_curve/final_camkII_state.dat' % p)})


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
        n = 0
        synChange = synapticChange('genericCase',plasticityCases[k])
        # initialize class which calculates the time the calcium trace spends above threshold
        tat = timeAboveThreshold(synChange.tauCa, synChange.Cpre, synChange.Cpost, synChange.thetaD, synChange.thetaP)
        # loop over range of deltaT values
        for i in deltaT:
                # fraction of time spent above threshold
                (alphaD[k,n],alphaP[k,n]) = tat.spikePairFrequency(i-synChange.D,frequency)
                # calculate all values for the change in synaptic strength
                synChange.changeInSynapticStrength(N*interval,0.5,alphaD[k,n],alphaP[k,n])
                
                # mean value of the synaptic strength right at the end of the stimulation protocol
                meanUP[k,n]   =  synChange.meanUP
                meanDOWN[k,n] =  synChange.meanDOWN
                
                # change in synaptic strength after/before   (self, T_total,rho0,alphaD,alphaP)
                sChange[k,n] = synChange.synChange
                # 
                n+=1
        

##########################################################
# generate figure 

fig_width = 10 # width in inches
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
gs = gridspec.GridSpec(3, 2 )

# define vertical and horizontal spacing between panels
gs.update(wspace=0.3,hspace=0.4)

# possibly change outer margins of the figure
#plt.subplots_adjust(left=0.14, right=0.92, top=0.92, bottom=0.18)



for k in range(nCases):
        ax0 = plt.subplot(gs[k])

        # title
        ax0.set_title(plasticityCases[k])

        # diplay of data
        ax0.axhline(y=1,ls='--',color='0.5',lw=2)
        ax0.axvline(x=0,ls='--',color='0.5',lw=2)
        ax0.plot(deltaT*1000., sChange[k], lw=3, c='m')

        ax0.fill_between(simData[plasticityCases[k]][:,0],simData[plasticityCases[k]][:,4]+simData[plasticityCases[k]][:,5],simData[plasticityCases[k]][:,4]-simData[plasticityCases[k]][:,5],color='C0',alpha=0.3)
        ax0.plot(simData[plasticityCases[k]][:,0],simData[plasticityCases[k]][:,4],'o-',ms=3,c='C0',lw=1)



        # removes upper and right axes 
        # and moves left and bottom axes away
        ax0.spines['top'].set_visible(False)
        ax0.spines['right'].set_visible(False)
        ax0.spines['bottom'].set_position(('outward', 10))
        ax0.spines['left'].set_position(('outward', 10))
        ax0.yaxis.set_ticks_position('left')
        ax0.xaxis.set_ticks_position('bottom')
        
        if plasticityCases[k] == 'DPD':
                ax0.set_xlim(-40,40)
        elif plasticityCases[k] == 'D':
                ax0.set_xlim(-15,15)
        else:
                ax0.set_xlim(DeltaTstart*1000.,DeltaTend*1000.)

        ax0.xaxis.set_major_locator(MaxNLocator(5))
        #ax0.yaxis.set_major_locator(MaxNLocator(6))
        
        if not k%2:
                plt.ylabel('change in synaptic strength')
        if k>=4:
                plt.xlabel(r'$\Delta t$ (ms)')
fname = os.path.basename(__file__)

savefig(outputDir+fname[:-2]+'png')
savefig(outputDir+fname[:-2]+'pdf')


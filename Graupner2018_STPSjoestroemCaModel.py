######################
## Script to calculate cross-correlations in voltage traces
######################

from scipy import *
from numpy import *
from numpy.fft import *
from math import *
from pylab import *
import matplotlib.pyplot as plt
#import brian
import scipy.signal
import scipy.stats
import os
import pdb
import random
import sys
import time
import pickle
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mc
from scipy.interpolate import interp1d
from scipy import interpolate

from timeAboveThreshold.timeAboveThreshold import *
from synapticChange import synapticChange

darkyellow = '#ff9f00'

##############################################################

par = 'stpJesperCaModel'
dS     = 'jesper'

nonlinear = 1.
w0 = 0.5

#
synChange = synapticChange(par,fromFile=True,nonlinear=nonlinear,dataSet=dS)

# initialize class which calculates the time the calcium trace spends above threshold
tat = timeAboveThreshold(synChange.tauCa, synChange.Cpre, synChange.Cpost, synChange.thetaD, synChange.thetaP,nonlinear=nonlinear)
                

##############################################################
# calculate solution
freq = linspace(0.1,50.,500)
sChange = zeros((len(freq),2))
sChangeStoch = zeros((len(freq)))
#synChangeStochOld = zeros((len(freq)))

deltaTs = linspace(-0.05,0.05,1001)
sChange2 = zeros(len(deltaTs))

for i in range(len(freq)):
    #
    # pre-post delta t = 10 ms
    # here alphaD and alphaP are actually the absolute times spent above threshold
    (alphaD, alphaP) = tat.spikePairFrequencySTP(0.01 - synChange.D, freq[i], synChange.Npairs, synChange.tauRec,synChange.U)
    synChange.changeInSynapticStrength(synChange.Npresentations, w0, alphaD, alphaP)
    #(alphaD,alphaP) = tat.spikePairFrequencyNonlinear(0.01-synChange.D,freq[i])
    #synChange.changeInSynapticStrength(Npresentations/freq[i],w0,alphaD,alphaP)
    sChange[i,0] = synChange.mean/w0
    #
    # post-pre delta t = -10 ms
    (alphaD, alphaP) = tat.spikePairFrequencySTP(-0.01 - synChange.D, freq[i], synChange.Npairs, synChange.tauRec,synChange.U)
    synChange.changeInSynapticStrength(synChange.Npresentations, w0, alphaD, alphaP)
    #(alphaD,alphaP) = tat.spikePairFrequencyNonlinear(-0.01-synChange.D,freq[i])
    #synChange.changeInSynapticStrength(Npresentations/freq[i],w0,alphaD,alphaP)
    sChange[i,1] = synChange.mean/w0
    #
    # stochastic sjoestroem protocol
    # spikePairFrequencySTPStochastic(self, DeltaTStart, DeltaTEnd, freq, Npres, tauRec, U):
    (alphaD,alphaP) = tat.spikePairFrequencySTPStochastic((-0.015-synChange.D),(0.015-synChange.D),freq[i],synChange.Npairs, synChange.tauRec,synChange.U)
    synChange.changeInSynapticStrength(synChange.Npresentations,w0,alphaD,alphaP)
    sChangeStoch[i] = synChange.mean/w0
    
for i in range(len(deltaTs)):
    (alphaD, alphaP) = tat.spikePairFrequencySTP(deltaTs[i] - synChange.D, 0.1 , synChange.Npairs, synChange.tauRec,synChange.U)
    synChange.changeInSynapticStrength(synChange.Npresentations, w0, alphaD, alphaP)
    #(alphaD,alphaP) = tat.spikePairFrequencyNonlinear(deltaTs[i]-synChange.D,0.1)
    #synChange.changeInSynapticStrength(Npresentations/0.1,0.5,alphaD,alphaP)
    sChange2[i] = synChange.mean/w0


##############################################################
fig_width = 16  # width in inches
fig_height = 16  # height in inches
fig_size =  [fig_width,fig_height]
params = {'axes.labelsize': 14,
          'axes.titlesize': 14,
          'font.size': 14,
          'xtick.labelsize': 14,
          'ytick.labelsize': 14,
          'figure.figsize': fig_size,
          'savefig.dpi' : 600,
          'axes.linewidth' : 1.3,
          'ytick.major.size' : 6,      # major tick size in points
          'xtick.major.size' : 6      # major tick size in points
          #'xtick.major.size' : 2,
          #'ytick.major.size' : 2,
          }
rcParams.update(params)

# set sans-serif font to Arial
#rcParams['font.sans-serif'] = 'Arial'
#rcParams['text.usetex'] = True
#rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

fig = plt.figure()

#f = plt.figure()

gs = gridspec.GridSpec(3, 1,#,
                       #width_ratios=[1,2],
                       height_ratios=[1,2,1]
                       )

#fig.subplots_adjust(hspace=0.6)
#ax = plt.axes((0.3, 0.12, 0.65, 0.8))

#fig = plt.figure(figsize=(16,4))
fig.subplots_adjust(hspace=0.4,wspace=0.3)
plt.subplots_adjust(left=0.07, right=0.97, top=0.96, bottom=0.07)


#gssub0 = gridspec.GridSpecFromSubplotSpec(5, 1, subplot_spec=gs[0],hspace=0.2)


#plt.figtext(0.01, 0.93, 'A',clip_on=False,color='black', weight='bold',size=22)
#plt.figtext(0.51, 0.93, 'B',clip_on=False,color='black', weight='bold',size=22)
#plt.figtext(0.01, 0.66, 'C',clip_on=False,color='black', weight='bold',size=22)
#plt.figtext(0.5, 0.66, 'D',clip_on=False,color='black', weight='bold',size=22)
#plt.figtext(0.01, 0.335, 'E',clip_on=False,color='black', weight='bold',size=22)
#plt.figtext(0.5, 0.335, 'F',clip_on=False,color='black', weight='bold',size=22)
#plt.figtext(0.08, 0.36, 'E',clip_on=False,color='black', weight='bold',size=22)

gssub0 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[0],hspace=0.2,wspace=0.2)

##############################################################
# regular sjoestroem data 
ax0 = plt.subplot(gssub0[0])
        

ax0.set_title('Sjoestroem : regular pairs (75 pairs)')
ax0.axhline(y=1.,c='0.7')
plusMask = synChange.xDataReg[:,1] > 0.
minusMask = synChange.xDataReg[:,1] < 0.
ax0.errorbar(synChange.xDataReg[:,0][plusMask],synChange.yDataReg[plusMask],yerr=synChange.sigmaDataReg[plusMask],fmt='s',color='red',clip_on=False)
ax0.errorbar(synChange.xDataReg[:,0][minusMask],synChange.yDataReg[minusMask],yerr=synChange.sigmaDataReg[minusMask],fmt='o',color='blue',clip_on=False)
ax0.plot(freq,sChange[:,0],lw=2,color='red')
ax0.plot(freq,sChange[:,1],lw=2,color='blue')


ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.spines['bottom'].set_position(('outward', 10))
ax0.spines['left'].set_position(('outward', 10))
ax0.yaxis.set_ticks_position('left')
ax0.xaxis.set_ticks_position('bottom')

plt.xlim(0,50)

plt.ylabel('change in synaptic strength')

plt.xlabel(r'frequency (Hz)')

#plt.legend(frameon=False,loc=4)


##############################################################
# regular sjoestroem data 
ax4 = plt.subplot(gssub0[1])
ax4.set_title('Sjoestroem : stochastic pairs (75 pairs)')

# diplay of data
ax4.axhline(y=1.,c='0.7')
ax4.errorbar(synChange.xDataStoch,synChange.yDataStoch,yerr=synChange.sigmaDataStoch,fmt='s',color='green',clip_on=False)
ax4.plot(freq,sChangeStoch,lw=2,color='green')
#ax4.plot(freq,synChangeStochOld,lw=2,color='turquoise')

# removes upper and right axes 
# and moves left and bottom axes away
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.spines['bottom'].set_position(('outward', 10))
ax4.spines['left'].set_position(('outward', 10))
ax4.yaxis.set_ticks_position('left')
ax4.xaxis.set_ticks_position('bottom')

# legends and labels
#plt.legend(loc=1,frameon=False)

plt.xlabel('frequency (Hz)')

########################################################
# stdp at low frequency 
ax1 = plt.subplot(gssub0[2])
ax1.set_title(r'Spike-pairs vs. $\Delta t$ (75 pairs)')

ax1.axhline(y=1.,c='0.7')
ax1.plot(deltaTs*1000.,sChange2,lw=2,color='k')

# removes upper and right axes 
# and moves left and bottom axes away
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_position(('outward', 10))
ax1.spines['left'].set_position(('outward', 10))
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')

# legends and labels
#plt.legend(loc=1,frameon=False)

plt.xlabel(r'$\Delta t$ (ms)')
#plt.ylabel('change in synaptic strength')


#plt.legend(frameon=False,loc=2)
#leg = plt.gca().get_legend()
#ltext  = leg.get_texts()
#plt.setp(ltext, fontsize=9)

gssub1 = gridspec.GridSpecFromSubplotSpec(2, 3, subplot_spec=gs[1],hspace=0.2,wspace=0.2)

dir0 = 'simResults/'


irrPairsDT = np.load(dir0+'irregularSpikePairs_vs_deltaT_differentFreqs_%s.npy' % par)
regPairsDT = np.load(dir0+'regularSpikePairs_vs_deltaT_differentFreqs_%s.npy' % par)


argM = argmin((irrPairsDT[:,0]+0.01)**2)
argP = argmin((irrPairsDT[:,0]-0.01)**2)

ul = 1.8
ll = 0.5

ax80 = plt.subplot(gssub1[0])

idx = 1
ax80.set_title('%s Hz' % irrPairsDT[0,idx])
ax80.axhline(y=1,color='0.4',ls='--')
ax80.axvline(x=0,color='0.4',ls='--')

mask = (regPairsDT[:,0] > -1./(2.*irrPairsDT[0,idx])) & (regPairsDT[:,0] < 1./(2.*irrPairsDT[0,idx]))
ax80.plot(irrPairsDT[:,0]*1000.,irrPairsDT[:,14]/w0,lw=2,c=darkyellow)
ax80.plot(regPairsDT[:,0][mask]*1000.,regPairsDT[:,14][mask]/w0,lw=2,c='blue')

ax80.spines['top'].set_visible(False)
ax80.spines['right'].set_visible(False)
ax80.spines['bottom'].set_visible(False)
ax80.spines['left'].set_position(('outward', 10))
ax80.yaxis.set_ticks_position('left')
ax80.axes.get_xaxis().set_visible(False)
plt.xlim(-80,80)
plt.ylim(ll,ul)

ax81 = plt.subplot(gssub1[1])

idx = 2
ax81.set_title(r'irregular pairs vs. $\Delta t$ (10 sec)'+'\n'+ '%s Hz' % irrPairsDT[0,idx])
ax81.axhline(y=1,color='0.4',ls='--')
#ax0.axhline(y=bb/(aa*0.5),color='black')
ax81.axvline(x=0,color='0.4',ls='--')


mask = (regPairsDT[:,0] > -1./(2.*irrPairsDT[0,idx])) & (regPairsDT[:,0] < 1./(2.*irrPairsDT[0,idx]))
ax81.plot(irrPairsDT[:,0]*1000.,irrPairsDT[:,15]/w0,lw=2,c=darkyellow)
ax81.plot(regPairsDT[:,0][mask]*1000.,regPairsDT[:,15][mask]/w0,lw=2,c='blue')

ax81.spines['top'].set_visible(False)
ax81.spines['right'].set_visible(False)
ax81.spines['bottom'].set_visible(False)
ax81.spines['left'].set_visible(False)
ax81.axes.get_xaxis().set_visible(False)
ax81.axes.get_yaxis().set_visible(False)
plt.xlim(-80,80)
plt.ylim(ll,ul)

ax82 = plt.subplot(gssub1[2])

idx = 3
ax82.set_title('%s Hz' % irrPairsDT[0,idx])
ax82.axhline(y=1,color='0.4',ls='--')
#ax0.axhline(y=bb/(aa*0.5),color='black')
ax82.axvline(x=0,color='0.4',ls='--')


mask = (regPairsDT[:,0] > -1./(2.*irrPairsDT[0,idx])) & (regPairsDT[:,0] < 1./(2.*irrPairsDT[0,idx]))
ax82.plot(irrPairsDT[:,0]*1000.,irrPairsDT[:,16]/w0,lw=2,c=darkyellow)
ax82.plot(regPairsDT[:,0][mask]*1000.,regPairsDT[:,16][mask]/w0,lw=2,c='blue')
#ax82.plot(irrPairsDT[argM,0]*1000.,irrPairsDT[:,16][argM]/w0,'o',c='k')
#ax82.plot(irrPairsDT[argP,0]*1000.,irrPairsDT[:,16][argP]/w0,'o',c='k')


ax82.spines['top'].set_visible(False)
ax82.spines['right'].set_visible(False)
ax82.spines['bottom'].set_visible(False)
ax82.spines['left'].set_visible(False)
ax82.axes.get_xaxis().set_visible(False)
ax82.axes.get_yaxis().set_visible(False)
plt.xlim(-80,80)
plt.ylim(ll,ul)

ax83 = plt.subplot(gssub1[3])

idx=4

ax83.set_title('%s Hz' % irrPairsDT[0,idx])
ax83.axhline(y=1,color='0.4',ls='--')
#ax0.axhline(y=bb/(aa*0.5),color='black')
ax83.axvline(x=0,color='0.4',ls='--')


mask = (regPairsDT[:,0] > -1./(2.*irrPairsDT[0,idx])) & (regPairsDT[:,0] < 1./(2.*irrPairsDT[0,idx]))
ax83.plot(irrPairsDT[:,0]*1000.,irrPairsDT[:,17]/w0,lw=2,c=darkyellow)
ax83.plot(regPairsDT[:,0][mask]*1000.,regPairsDT[:,17][mask]/w0,lw=2,c='blue')

ax83.spines['top'].set_visible(False)
ax83.spines['right'].set_visible(False)
ax83.spines['bottom'].set_position(('outward', 10))
ax83.spines['left'].set_position(('outward', 10))
ax83.yaxis.set_ticks_position('left')
ax83.xaxis.set_ticks_position('bottom')
plt.xlim(-80,80)
plt.ylim(ll,ul)
plt.xticks([-50, 0, 50],[-50, 0, 50])

plt.ylabel('change in synaptic strength',position=(1,1.15))

ax84 = plt.subplot(gssub1[4])

idx = 5

ax84.set_title('%s Hz' % irrPairsDT[0,idx])
ax84.axhline(y=1,color='0.4',ls='--')
#ax0.axhline(y=bb/(aa*0.5),color='black')
ax84.axvline(x=0,color='0.4',ls='--')

mask = (regPairsDT[:,0] > -1./(2.*irrPairsDT[0,idx])) & (regPairsDT[:,0] < 1./(2.*irrPairsDT[0,idx]))
ax84.plot(irrPairsDT[:,0]*1000.,irrPairsDT[:,18]/w0,lw=2,c=darkyellow)
ax84.plot(regPairsDT[:,0][mask]*1000.,regPairsDT[:,18][mask]/w0,lw=2,c='blue')
#ax84.plot(irrPairsDT[argM,0]*1000.,irrPairsDT[:,18][argM]/w0,'>',c='k')
#ax84.plot(irrPairsDT[argP,0]*1000.,irrPairsDT[:,18][argP]/w0,'>',c='k')

ax84.spines['top'].set_visible(False)
ax84.spines['right'].set_visible(False)
ax84.spines['bottom'].set_position(('outward', 10))
ax84.spines['left'].set_visible(False)
ax84.xaxis.set_ticks_position('bottom')
ax84.axes.get_yaxis().set_visible(False)
plt.xlim(-80,80)
plt.ylim(ll,ul)
plt.xticks([-50, 0, 50],[-50, 0, 50])

plt.xlabel(r'$\Delta t$ (msec)')

ax85 = plt.subplot(gssub1[5])

idx = 6

ax85.set_title('%s Hz' % irrPairsDT[0,idx])
ax85.axhline(y=1,color='0.4',ls='--')
ax85.axvline(x=0,color='0.4',ls='--')


mask = (regPairsDT[:,0] > -1./(2.*irrPairsDT[0,idx])) & (regPairsDT[:,0] < 1./(2.*irrPairsDT[0,idx]))
ax85.plot(irrPairsDT[:,0]*1000.,irrPairsDT[:,19]/w0,lw=2,c=darkyellow)
ax85.plot(regPairsDT[:,0][mask]*1000.,regPairsDT[:,19][mask]/w0,lw=2,c='blue')

ax85.spines['top'].set_visible(False)
ax85.spines['right'].set_visible(False)
ax85.spines['bottom'].set_position(('outward', 10))
ax85.spines['left'].set_visible(False)
ax85.xaxis.set_ticks_position('bottom')
ax85.axes.get_yaxis().set_visible(False)
plt.xlim(-80,80)
plt.ylim(ll,ul)
plt.xticks([-50, 0, 50],[-50, 0, 50])

#################################################################################
gssub2 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[2],hspace=0.2,wspace=0.3)

irrPairsFreq = np.load(dir0+'irregularSpikePairs_vs_rate_differentDeltaTs_%s.npy' % par)


ax12 = plt.subplot(gssub2[0])
ax12.set_title('irregular pairs vs. rate (10 sec)')

ax12.axhline(y=1,ls='--',color='0.6',lw=2)

ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,16]/w0,lw=2,c='red',label=r'$\Delta t = %s$ ms, $\rho=%s$' % (irrPairsFreq[0,4]*1000.,irrPairsFreq[0,11]))
ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,15]/w0,lw=2,c='red',alpha=0.6,label=r'$\Delta t = %s$ ms, $\rho=%s$' % (irrPairsFreq[0,3]*1000.,irrPairsFreq[0,10]))
ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,14]/w0,lw=2,c='k',label=r'$\rho=%s$' % (irrPairsFreq[0,9]))
ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,13]/w0,lw=2,c='blue',alpha=0.6,label=r'$\Delta t = %s$ ms, $\rho=%s$' % (irrPairsFreq[0,1]*1000.,irrPairsFreq[0,8]))
ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,12]/w0,lw=2,c='blue',label=r'$\Delta t = %s$ ms, $\rho=%s$' % (irrPairsFreq[0,0]*1000.,irrPairsFreq[0,7]))
#ax12.plot(irrPairsDT[:,3][argM],irrPairsDT[:,16][argM]/w0,'o',c='k')
#ax12.plot(irrPairsDT[:,3][argP],irrPairsDT[:,16][argP]/w0,'o',c='k')
#ax12.plot(irrPairsDT[:,5][argM],irrPairsDT[:,18][argM]/w0,'>',c='k')
#ax12.plot(irrPairsDT[:,5][argP],irrPairsDT[:,18][argP]/w0,'>',c='k')

ax12.spines['top'].set_visible(False)
ax12.spines['right'].set_visible(False)
ax12.spines['bottom'].set_position(('outward', 10))
ax12.spines['left'].set_position(('outward', 10))
ax12.yaxis.set_ticks_position('left')
ax12.xaxis.set_ticks_position('bottom')

ax12.set_ylabel('change in synaptic strength')
ax12.set_xlabel('rate (spk/sec)')

####################
########################################
ax1 = plt.subplot(gssub2[1])


ax1.axhline(y=0,color='0.4',ls='--')

ax1.plot(irrPairsFreq[:,5],(irrPairsFreq[:,16]-irrPairsFreq[:,14])/w0,color='r',lw=3,label=r'$\rho=0.4$, $\Delta = 10$ ms')
ax1.plot(irrPairsFreq[:,5],(irrPairsFreq[:,15]-irrPairsFreq[:,14])/w0,c='red',alpha=0.6,lw=3,label=r'$\rho=0.2$, $\Delta = 10$ ms')

ax1.plot(irrPairsFreq[:,5],(irrPairsFreq[:,12]-irrPairsFreq[:,14])/w0,color='b',lw=3,label=r'$\rho=0.2$, $\Delta = -10$ ms')
ax1.plot(irrPairsFreq[:,5],(irrPairsFreq[:,13]-irrPairsFreq[:,14])/w0,c='blue',alpha=0.6,lw=3,label=r'$\rho=0.4$, $\Delta = -10$ ms')

ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_position(('outward', 10))
ax1.spines['left'].set_position(('outward', 10))
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')
ax1.set_xlim(0,80)
#plt.ylim(lls,uls)
#plt.xticks([-50, 0, 50],[-50, 0, 50])
ax1.set_ylabel('sensitivity to correlations')
ax1.set_xlabel('firing rate (spk/s)')

#majorLocator_y = MaxNLocator(4)
ax1.yaxis.set_major_locator(MaxNLocator(5))
ax1.xaxis.set_major_locator(MaxNLocator(4))

############################
#########################################################################################################
# Calcium-model firing rate sensitivity
ax1 = plt.subplot(gssub2[2])

NN = 1001

rateMatrix = zeros((NN, NN))
# pdb.set_trace()
for i in range(NN):
    rateMatrix[:, i] = linspace(0.1, 2. * 40., 2 * NN)[i:(i + NN)]

# pdb.set_trace()

#####New piece of code #####################################################
# rho0_interp=interp1d(stdpanc[:,0],stdpanc[:,19]/wcal0)

# caSurface = rho0_interp(rateMatrix)-rho0_interp(rateMatrix[:,0])

minC = -0.1
maxC = 0.55

rho04_interp = interp1d(irrPairsFreq[:,5], (irrPairsFreq[:,16] - irrPairsFreq[:,14]) / w0)
rho02_interp = interp1d(irrPairsFreq[:,5], (irrPairsFreq[:,15] - irrPairsFreq[:,14]) / w0)
rho0_interp = interp1d(irrPairsFreq[:,5], irrPairsFreq[:,14] / w0)

caSurface = rho0_interp(rateMatrix) - rho0_interp(rateMatrix[:, 0])
corr04_line = rho04_interp(rateMatrix[:,0])
corr02_line = rho02_interp(rateMatrix[:,0])

corr04_match = argmin(abs(caSurface - corr04_line), axis=0)
corr02_match = argmin(abs(caSurface - corr02_line), axis=0)

f04CaMatch = rateMatrix[:, 0][corr04_match]
f02CaMatch = rateMatrix[:, 0][corr02_match]

# ax1.axhline(y=3.,ls='--',c='w',lw=1)
ax1.plot(rateMatrix[:, 0], f04CaMatch, lw=3, c='r')
ax1.plot(rateMatrix[:, 0], f02CaMatch, lw=3, c='red',alpha=0.6)
ax1c = ax1.imshow(caSurface,
                  extent=(rateMatrix[:, 0][0], rateMatrix[:, 0][-1], rateMatrix[:, 0][0], rateMatrix[:, 0][-1]),
                  aspect=1, origin='lower')
cb = plt.colorbar(ax1c, shrink=0.9)
cb.locator = plt.MaxNLocator(nbins=5)
cb.update_ticks()

# ax1.imshow(ax1.plot(stdpanc[:,0][:90],rho1,"m--",label='5+50% rate increase \n no correlations',lw=3)
# ax1.plot(stdpanc[:,0][:90],rho11,"g--",label='50% rate increase \n no correlations',lw=3)
###########################################################################
# plt.xscale('log')
# plt.yscale('log')
# ax1.text(32, 38, "100 %",color='w', va="center", ha="center")
# ax1.text(35, 9, "20 %",color='w', va="center", ha="center")

ax1.text(32, 42., '     change in \nsynaptic strength')

ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_position(('outward', 10))
ax1.spines['left'].set_position(('outward', 10))
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')
plt.xlim(0, 40)
plt.ylim(0, 40)
# plt.xticks([-50, 0, 50],[-50, 0, 50])
plt.ylabel('increase in firing rate (spk/s)', multialignment='center')
plt.xlabel('baseline firing rate (spk/s)')

# majorLocator_y = MaxNLocator(4)
# ax1.yaxis.set_major_locator(majorLocator_y)
majorLocator_y = MaxNLocator(4)
ax1.yaxis.set_major_locator(majorLocator_y)
ax1.xaxis.set_major_locator(majorLocator_y)



#ax12.set_xlim(0,40)

#tauCa = sol[0][0]
#Cpre  = sol[0][1]
#Cpost = sol[0][2]
#thetaD = 1.
#thetaP = 1.3
#gammaD = sol[0][3]
#gammaP = sol[0][4]
#sigma  = 1.
#tau    = sol[0][5]
#rhoStar= 0.5
#D      = sol[0][6]
#beta   = 0.5
#b      = 2.


#ttt = r'$\mathbf{ mse = %s}$' + '\n' + r'$\tau_{\rm Ca} = %s$' + '\n' + r'$C_{\rm pre} = %s$' + '\n' + r'$C_{\rm post} = %s$' + '\n' + r'$\theta_d = %s$' + '\n' + r'$\theta_p = %s$' + '\n' + r'$\gamma_d = %s$' + '\n' + r'$\gamma_p = %s$' + '\n' + r'$\tau = %s$' + '\n' + r'$D = %s$'

#synChange.tauCa, synChange.Cpre, synChange.Cpost, synChange.thetaD, synChange.thetaP
#ax12.text(70,1.1,ttt % (sol[1],sol[0][0],sol[0][1],sol[0][2],thetaD,sol[0][3],sol[0][4],sol[0][5],sol[0][6],sol[0][7]))
#ax12.text(90,1.1,ttt % (synChange.mse,synChange.tauCa,synChange.Cpre,synChange.Cpost,synChange.thetaD,synChange.thetaP,synChange.gammaD,synChange.gammaP,synChange.tau,synChange.D))

#lll = plt.legend(frameon=True,loc=4)
#frame = lll.get_frame()
#frame.set_edgecolor('white')

#plt.xlim(0,100)



###############################################

fname = os.path.basename(__file__)

#savefig('ccgs_constant_coupled_networks.svg')
plt.savefig('outputFigures/'+fname[:-3]+'.pdf')
plt.savefig('outputFigures/'+fname[:-3]+'.png')
#plt.savefig(fname[:-3]+'_%s.png' % par)
#plt.savefig('delta_t.png')

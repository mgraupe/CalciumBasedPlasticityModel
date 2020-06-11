######################
## Script to calculate cross-correlations in voltage traces
######################

from scipy import *
from numpy import *
from numpy.fft import *
from math import *
from pylab import *
import matplotlib.pyplot as plt
# import brian
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

parameterSetVis = 'sJFullNoSTDSim0'

nonlinear = 1.
w0 = 0.5

#
#synChangeVis = synapticChange(parameterSetVis, fromFile=True, nonlinear=nonlinear)

# initialize class which calculates the time the calcium trace spends above threshold
#tatVis = timeAboveThreshold(synChangeVis.tauCa, synChangeVis.Cpre, synChangeVis.Cpost, synChangeVis.thetaD, synChangeVis.thetaP, nonlinear=nonlinear)


parameterSetSom = 'sHFullNoSTDSim1b'
#dS     = 'jesper'

#synChangeSom = synapticChange(parameterSetSom,fromFile=True,nonlinear=nonlinear)

# initialize class which calculates the time the calcium trace spends above threshold
#tatSom = timeAboveThreshold(synChangeSom.tauCa, synChangeSom.Cpre, synChangeSom.Cpost, synChangeSom.thetaD, synChangeSom.thetaP,nonlinear=nonlinear)


#############################################################
fig_width = 10  # width in inches
fig_height = 10  # height in inches
fig_size = [fig_width, fig_height]
params = {'axes.labelsize': 14, 'axes.titlesize': 14, 'font.size': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'figure.figsize': fig_size, 'savefig.dpi': 600, 'axes.linewidth': 1.3,
          'ytick.major.size': 6,  # major tick size in points
          'xtick.major.size': 6  # major tick size in points
          # 'xtick.major.size' : 2,
          # 'ytick.major.size' : 2,
          }
rcParams.update(params)

# set sans-serif font to Arial
# rcParams['font.sans-serif'] = 'Arial'
# rcParams['text.usetex'] = True
# rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

fig = plt.figure()

# f = plt.figure()

gs = gridspec.GridSpec(2, 1,  # ,
                       # width_ratios=[1,2],
                       #height_ratios=[1, 2, 1]
                       )

# fig.subplots_adjust(hspace=0.6)
# ax = plt.axes((0.3, 0.12, 0.65, 0.8))

# fig = plt.figure(figsize=(16,4))
fig.subplots_adjust(hspace=0.35, wspace=0.2)
plt.subplots_adjust(left=0.11, right=0.97, top=0.95, bottom=0.07)

# gssub0 = gridspec.GridSpecFromSubplotSpec(5, 1, subplot_spec=gs[0],hspace=0.2)


plt.figtext(0.01, 0.97, 'A',clip_on=False,color='black', weight='bold',size=24)
plt.figtext(0.01, 0.47, 'B',clip_on=False,color='black', weight='bold',size=24)
# plt.figtext(0.01, 0.66, 'C',clip_on=False,color='black', weight='bold',size=22)
# plt.figtext(0.5, 0.66, 'D',clip_on=False,color='black', weight='bold',size=22)
# plt.figtext(0.01, 0.335, 'E',clip_on=False,color='black', weight='bold',size=22)
# plt.figtext(0.5, 0.335, 'F',clip_on=False,color='black', weight='bold',size=22)
# plt.figtext(0.08, 0.36, 'E',clip_on=False,color='black', weight='bold',size=22)

###################################################################################################
gssub0 = gridspec.GridSpecFromSubplotSpec(2, 3, subplot_spec=gs[0], hspace=0.2, wspace=0.2)

dir0 = 'simResults/'

irrPairsDT = np.load(dir0+'irregularSpikePairs_vs_deltaT_differentFreqs_STDdet_%s.npy' % parameterSetVis)
regPairsDT = np.load(dir0 + 'regularSpikePairs_vs_deltaT_differentFreqs_STDdet_%s.npy' % parameterSetVis)

argM = argmin((irrPairsDT[:, 0] + 0.01) ** 2)
argP = argmin((irrPairsDT[:, 0] - 0.01) ** 2)

ul = 1.8
ll = 0.5

ax80 = plt.subplot(gssub0[0])

idx = 1
ax80.set_title('%s spk/s' % int(irrPairsDT[0, idx]))
ax80.axhline(y=1, color='0.4', ls='--')
ax80.axvline(x=0, color='0.4', ls='--')

mask = (regPairsDT[:, 0] > -1. / (2. * irrPairsDT[0, idx])) & (regPairsDT[:, 0] < 1. / (2. * irrPairsDT[0, idx]))
ax80.plot(irrPairsDT[:, 0] * 1000., irrPairsDT[:, 14] / w0, lw=2, c=darkyellow)
ax80.plot(regPairsDT[:, 0][mask] * 1000., regPairsDT[:, 14][mask] / w0, lw=2, c='blue')

ax80.spines['top'].set_visible(False)
ax80.spines['right'].set_visible(False)
ax80.spines['bottom'].set_visible(False)
ax80.spines['left'].set_position(('outward', 10))
ax80.yaxis.set_ticks_position('left')
ax80.axes.get_xaxis().set_visible(False)
plt.xlim(-80, 80)
plt.ylim(ll, ul)

ax81 = plt.subplot(gssub0[1])

idx = 2
ax81.set_title(r'visual cortex' + '\n' + '%s spk/s' % int(irrPairsDT[0, idx]))
ax81.axhline(y=1, color='0.4', ls='--')
# ax0.axhline(y=bb/(aa*0.5),color='black')
ax81.axvline(x=0, color='0.4', ls='--')

mask = (regPairsDT[:, 0] > -1. / (2. * irrPairsDT[0, idx])) & (regPairsDT[:, 0] < 1. / (2. * irrPairsDT[0, idx]))
ax81.plot(irrPairsDT[:, 0] * 1000., irrPairsDT[:, 15] / w0, lw=2, c=darkyellow)
ax81.plot(regPairsDT[:, 0][mask] * 1000., regPairsDT[:, 15][mask] / w0, lw=2, c='blue')

ax81.spines['top'].set_visible(False)
ax81.spines['right'].set_visible(False)
ax81.spines['bottom'].set_visible(False)
ax81.spines['left'].set_visible(False)
ax81.axes.get_xaxis().set_visible(False)
ax81.axes.get_yaxis().set_visible(False)
plt.xlim(-80, 80)
plt.ylim(ll, ul)

ax82 = plt.subplot(gssub0[2])

idx = 3
ax82.set_title('%s spk/s' % int(irrPairsDT[0, idx]))
ax82.axhline(y=1, color='0.4', ls='--')
# ax0.axhline(y=bb/(aa*0.5),color='black')
ax82.axvline(x=0, color='0.4', ls='--')

mask = (regPairsDT[:, 0] > -1. / (2. * irrPairsDT[0, idx])) & (regPairsDT[:, 0] < 1. / (2. * irrPairsDT[0, idx]))
ax82.plot(irrPairsDT[:, 0] * 1000., irrPairsDT[:, 16] / w0, lw=2, c=darkyellow)
ax82.plot(regPairsDT[:, 0][mask] * 1000., regPairsDT[:, 16][mask] / w0, lw=2, c='blue')
# ax82.plot(irrPairsDT[argM,0]*1000.,irrPairsDT[:,16][argM]/w0,'o',c='k')
# ax82.plot(irrPairsDT[argP,0]*1000.,irrPairsDT[:,16][argP]/w0,'o',c='k')


ax82.spines['top'].set_visible(False)
ax82.spines['right'].set_visible(False)
ax82.spines['bottom'].set_visible(False)
ax82.spines['left'].set_visible(False)
ax82.axes.get_xaxis().set_visible(False)
ax82.axes.get_yaxis().set_visible(False)
plt.xlim(-80, 80)
plt.ylim(ll, ul)

ax83 = plt.subplot(gssub0[3])

idx = 4

ax83.set_title('%s spk/s' % int(irrPairsDT[0, idx]))
ax83.axhline(y=1, color='0.4', ls='--')
# ax0.axhline(y=bb/(aa*0.5),color='black')
ax83.axvline(x=0, color='0.4', ls='--')

mask = (regPairsDT[:, 0] > -1. / (2. * irrPairsDT[0, idx])) & (regPairsDT[:, 0] < 1. / (2. * irrPairsDT[0, idx]))
ax83.plot(irrPairsDT[:, 0] * 1000., irrPairsDT[:, 17] / w0, lw=2, c=darkyellow)
ax83.plot(regPairsDT[:, 0][mask] * 1000., regPairsDT[:, 17][mask] / w0, lw=2, c='blue')

ax83.spines['top'].set_visible(False)
ax83.spines['right'].set_visible(False)
ax83.spines['bottom'].set_position(('outward', 10))
ax83.spines['left'].set_position(('outward', 10))
ax83.yaxis.set_ticks_position('left')
ax83.xaxis.set_ticks_position('bottom')
plt.xlim(-80, 80)
plt.ylim(ll, ul)
plt.xticks([-50, 0, 50], [-50, 0, 50])

plt.ylabel('change in synaptic strength', position=(1, 1.15))

ax84 = plt.subplot(gssub0[4])

idx = 5

ax84.set_title('%s spk/s' % int(irrPairsDT[0, idx]))
ax84.axhline(y=1, color='0.4', ls='--')
# ax0.axhline(y=bb/(aa*0.5),color='black')
ax84.axvline(x=0, color='0.4', ls='--')

mask = (regPairsDT[:, 0] > -1. / (2. * irrPairsDT[0, idx])) & (regPairsDT[:, 0] < 1. / (2. * irrPairsDT[0, idx]))
ax84.plot(irrPairsDT[:, 0] * 1000., irrPairsDT[:, 18] / w0, lw=2, c=darkyellow)
ax84.plot(regPairsDT[:, 0][mask] * 1000., regPairsDT[:, 18][mask] / w0, lw=2, c='blue')
# ax84.plot(irrPairsDT[argM,0]*1000.,irrPairsDT[:,18][argM]/w0,'>',c='k')
# ax84.plot(irrPairsDT[argP,0]*1000.,irrPairsDT[:,18][argP]/w0,'>',c='k')

ax84.spines['top'].set_visible(False)
ax84.spines['right'].set_visible(False)
ax84.spines['bottom'].set_position(('outward', 10))
ax84.spines['left'].set_visible(False)
ax84.xaxis.set_ticks_position('bottom')
ax84.axes.get_yaxis().set_visible(False)
plt.xlim(-80, 80)
plt.ylim(ll, ul)
plt.xticks([-50, 0, 50], [-50, 0, 50])

plt.xlabel(r'$\Delta t$ (msec)')

ax85 = plt.subplot(gssub0[5])

idx = 6

ax85.set_title('%s spk/s' % int(irrPairsDT[0, idx]))
ax85.axhline(y=1, color='0.4', ls='--')
ax85.axvline(x=0, color='0.4', ls='--')

mask = (regPairsDT[:, 0] > -1. / (2. * irrPairsDT[0, idx])) & (regPairsDT[:, 0] < 1. / (2. * irrPairsDT[0, idx]))
ax85.plot(irrPairsDT[:, 0] * 1000., irrPairsDT[:, 19] / w0, lw=2, c=darkyellow, label='Poisson')
ax85.plot(regPairsDT[:, 0][mask] * 1000., regPairsDT[:, 19][mask] / w0, lw=2, c='blue',label='regular')

ax85.spines['top'].set_visible(False)
ax85.spines['right'].set_visible(False)
ax85.spines['bottom'].set_position(('outward', 10))
ax85.spines['left'].set_visible(False)
ax85.xaxis.set_ticks_position('bottom')
ax85.axes.get_yaxis().set_visible(False)
plt.xlim(-80, 80)
plt.ylim(ll, ul)
plt.xticks([-50, 0, 50], [-50, 0, 50])

lll = plt.legend(loc=(0.35,-0.02),frameon=True)

# change legend text size 
frame = lll.get_frame()
frame.set_edgecolor('white')
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=13)

#############################################################

gssub1 = gridspec.GridSpecFromSubplotSpec(2, 3, subplot_spec=gs[1],hspace=0.2,wspace=0.2)

dir0 = 'simResults/'


irrPairsDT = np.load(dir0+'irregularSpikePairs_vs_deltaT_differentFreqs_STDdet_%s.npy' % parameterSetSom)
regPairsDT = np.load(dir0+'regularSpikePairs_vs_deltaT_differentFreqs_STDdet_%s.npy' % parameterSetSom)
#irrPairsDTStoch = np.load(dir0+'irregularSpikePairs_vs_deltaT_differentFreqs_STDstoch_%s.npy' % parameterSetSom)
#regPairsDTStoch = np.load(dir0+'regularSpikePairs_vs_deltaT_differentFreqs_STDstoch_%s.npy' % parameterSetSom)


argM = argmin((irrPairsDT[:,0]+0.01)**2)
argP = argmin((irrPairsDT[:,0]-0.01)**2)

ul = 1.8
ll = 0.1

ax80 = plt.subplot(gssub1[0])

idx = 1
ax80.set_title('%s spk/s' % int(irrPairsDT[0,idx]))
ax80.axhline(y=1,color='0.4',ls='--')
ax80.axvline(x=0,color='0.4',ls='--')

mask = (regPairsDT[:,0] > -1./(2.*irrPairsDT[0,idx])) & (regPairsDT[:,0] < 1./(2.*irrPairsDT[0,idx]))
#maskStoch = (regPairsDTStoch[:,0] > -1./(2.*irrPairsDTStoch[0,idx])) & (regPairsDTStoch[:,0] < 1./(2.*irrPairsDTStoch[0,idx]))
ax80.plot(irrPairsDT[:,0]*1000.,irrPairsDT[:,(13+idx)]/w0,lw=2,c=darkyellow)
ax80.plot(regPairsDT[:,0][mask]*1000.,regPairsDT[:,(13+idx)][mask]/w0,lw=2,c='blue')
#ax80.plot(irrPairsDTStoch[:,0]*1000.,irrPairsDTStoch[:,(13+idx)]/w0,lw=2,ls='--',c=darkyellow)
#ax80.plot(regPairsDTStoch[:,0][maskStoch]*1000.,regPairsDTStoch[:,(13+idx)][maskStoch]/w0,lw=2,ls='--',c='blue')

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
ax81.set_title(r'somatosensory cortex'+'\n'+ '%s spk/s' % int(irrPairsDT[0,idx]))
ax81.axhline(y=1,color='0.4',ls='--')
#ax0.axhline(y=bb/(aa*0.5),color='black')
ax81.axvline(x=0,color='0.4',ls='--')

mask = (regPairsDT[:,0] > -1./(2.*irrPairsDT[0,idx])) & (regPairsDT[:,0] < 1./(2.*irrPairsDT[0,idx]))
#maskStoch = (regPairsDTStoch[:,0] > -1./(2.*irrPairsDTStoch[0,idx])) & (regPairsDTStoch[:,0] < 1./(2.*irrPairsDTStoch[0,idx]))
ax81.plot(irrPairsDT[:,0]*1000.,irrPairsDT[:,(13+idx)]/w0,lw=2,c=darkyellow)
ax81.plot(regPairsDT[:,0][mask]*1000.,regPairsDT[:,(13+idx)][mask]/w0,lw=2,c='blue')
#ax81.plot(irrPairsDTStoch[:,0]*1000.,irrPairsDTStoch[:,(13+idx)]/w0,lw=2,ls='--',c=darkyellow)
#ax81.plot(regPairsDTStoch[:,0][maskStoch]*1000.,regPairsDTStoch[:,(13+idx)][maskStoch]/w0,lw=2,ls='--',c='blue')

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
ax82.set_title('%s spk/s' % int(irrPairsDT[0,idx]))
ax82.axhline(y=1,color='0.4',ls='--')
#ax0.axhline(y=bb/(aa*0.5),color='black')
ax82.axvline(x=0,color='0.4',ls='--')

mask = (regPairsDT[:,0] > -1./(2.*irrPairsDT[0,idx])) & (regPairsDT[:,0] < 1./(2.*irrPairsDT[0,idx]))
#maskStoch = (regPairsDTStoch[:,0] > -1./(2.*irrPairsDTStoch[0,idx])) & (regPairsDTStoch[:,0] < 1./(2.*irrPairsDTStoch[0,idx]))
ax82.plot(irrPairsDT[:,0]*1000.,irrPairsDT[:,(13+idx)]/w0,lw=2,c=darkyellow)
ax82.plot(regPairsDT[:,0][mask]*1000.,regPairsDT[:,(13+idx)][mask]/w0,lw=2,c='blue')
#ax82.plot(irrPairsDTStoch[:,0]*1000.,irrPairsDTStoch[:,(13+idx)]/w0,lw=2,ls='--',c=darkyellow)
#ax82.plot(regPairsDTStoch[:,0][maskStoch]*1000.,regPairsDTStoch[:,(13+idx)][maskStoch]/w0,lw=2,ls='--',c='blue')

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

ax83.set_title('%s spk/s' % int(irrPairsDT[0,idx]))
ax83.axhline(y=1,color='0.4',ls='--')
#ax0.axhline(y=bb/(aa*0.5),color='black')
ax83.axvline(x=0,color='0.4',ls='--')

mask = (regPairsDT[:,0] > -1./(2.*irrPairsDT[0,idx])) & (regPairsDT[:,0] < 1./(2.*irrPairsDT[0,idx]))
#maskStoch = (regPairsDTStoch[:,0] > -1./(2.*irrPairsDTStoch[0,idx])) & (regPairsDTStoch[:,0] < 1./(2.*irrPairsDTStoch[0,idx]))
ax83.plot(irrPairsDT[:,0]*1000.,irrPairsDT[:,(13+idx)]/w0,lw=2,c=darkyellow)
ax83.plot(regPairsDT[:,0][mask]*1000.,regPairsDT[:,(13+idx)][mask]/w0,lw=2,c='blue')
#ax83.plot(irrPairsDTStoch[:,0]*1000.,irrPairsDTStoch[:,(13+idx)]/w0,lw=2,ls='--',c=darkyellow)
#ax83.plot(regPairsDTStoch[:,0][maskStoch]*1000.,regPairsDTStoch[:,(13+idx)][maskStoch]/w0,lw=2,ls='--',c='blue')

#mask = (regPairsDT[:,0] > -1./(2.*irrPairsDT[0,idx])) & (regPairsDT[:,0] < 1./(2.*irrPairsDT[0,idx]))
#ax83.plot(irrPairsDT[:,0]*1000.,irrPairsDT[:,17]/w0,lw=2,c=darkyellow)
#ax83.plot(regPairsDT[:,0][mask]*1000.,regPairsDT[:,17][mask]/w0,lw=2,c='blue')

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

ax84.set_title('%s spk/s' % int(irrPairsDT[0,idx]))
ax84.axhline(y=1,color='0.4',ls='--')
#ax0.axhline(y=bb/(aa*0.5),color='black')
ax84.axvline(x=0,color='0.4',ls='--')

mask = (regPairsDT[:,0] > -1./(2.*irrPairsDT[0,idx])) & (regPairsDT[:,0] < 1./(2.*irrPairsDT[0,idx]))
#maskStoch = (regPairsDTStoch[:,0] > -1./(2.*irrPairsDTStoch[0,idx])) & (regPairsDTStoch[:,0] < 1./(2.*irrPairsDTStoch[0,idx]))
ax84.plot(irrPairsDT[:,0]*1000.,irrPairsDT[:,(13+idx)]/w0,lw=2,c=darkyellow)
ax84.plot(regPairsDT[:,0][mask]*1000.,regPairsDT[:,(13+idx)][mask]/w0,lw=2,c='blue')
#ax84.plot(irrPairsDTStoch[:,0]*1000.,irrPairsDTStoch[:,(13+idx)]/w0,lw=2,ls='--',c=darkyellow)
#ax84.plot(regPairsDTStoch[:,0][maskStoch]*1000.,regPairsDTStoch[:,(13+idx)][maskStoch]/w0,lw=2,ls='--',c='blue')

#mask = (regPairsDT[:,0] > -1./(2.*irrPairsDT[0,idx])) & (regPairsDT[:,0] < 1./(2.*irrPairsDT[0,idx]))
#ax84.plot(irrPairsDT[:,0]*1000.,irrPairsDT[:,18]/w0,lw=2,c=darkyellow)
#ax84.plot(regPairsDT[:,0][mask]*1000.,regPairsDT[:,18][mask]/w0,lw=2,c='blue')
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

ax85.set_title('%s spk/s' % int(irrPairsDT[0,idx]))
ax85.axhline(y=1,color='0.4',ls='--')
ax85.axvline(x=0,color='0.4',ls='--')

mask = (regPairsDT[:,0] > -1./(2.*irrPairsDT[0,idx])) & (regPairsDT[:,0] < 1./(2.*irrPairsDT[0,idx]))
#maskStoch = (regPairsDTStoch[:,0] > -1./(2.*irrPairsDTStoch[0,idx])) & (regPairsDTStoch[:,0] < 1./(2.*irrPairsDTStoch[0,idx]))
ax85.plot(irrPairsDT[:,0]*1000.,irrPairsDT[:,(13+idx)]/w0,lw=2,c=darkyellow)
ax85.plot(regPairsDT[:,0][mask]*1000.,regPairsDT[:,(13+idx)][mask]/w0,lw=2,c='blue')
#ax85.plot(irrPairsDTStoch[:,0]*1000.,irrPairsDTStoch[:,(13+idx)]/w0,lw=2,ls='--',c=darkyellow)
#ax85.plot(regPairsDTStoch[:,0][maskStoch]*1000.,regPairsDTStoch[:,(13+idx)][maskStoch]/w0,lw=2,ls='--',c='blue')

#mask = (regPairsDT[:,0] > -1./(2.*irrPairsDT[0,idx])) & (regPairsDT[:,0] < 1./(2.*irrPairsDT[0,idx]))
#ax85.plot(irrPairsDT[:,0]*1000.,irrPairsDT[:,19]/w0,lw=2,c=darkyellow)
#ax85.plot(regPairsDT[:,0][mask]*1000.,regPairsDT[:,19][mask]/w0,lw=2,c='blue')

ax85.spines['top'].set_visible(False)
ax85.spines['right'].set_visible(False)
ax85.spines['bottom'].set_position(('outward', 10))
ax85.spines['left'].set_visible(False)
ax85.xaxis.set_ticks_position('bottom')
ax85.axes.get_yaxis().set_visible(False)
plt.xlim(-80,80)
plt.ylim(ll,ul)
plt.xticks([-50, 0, 50],[-50, 0, 50])

###########################################################



fname = os.path.basename(__file__)

# savefig('ccgs_constant_coupled_networks.svg')
plt.savefig('outputFigures/' + fname[:-3] + '_v2.pdf')
plt.savefig('outputFigures/' + fname[:-3] + '_v2.png')  # plt.savefig(fname[:-3]+'_%s.png' % par)
# plt.savefig('delta_t.png')

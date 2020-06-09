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


parameterSetSom = 'sHFullNoSTDSim1'
#dS     = 'jesper'

dir0 = 'simResults/'

#synChangeSom = synapticChange(parameterSetSom,fromFile=True,nonlinear=nonlinear)

# initialize class which calculates the time the calcium trace spends above threshold
#tatSom = timeAboveThreshold(synChangeSom.tauCa, synChangeSom.Cpre, synChangeSom.Cpost, synChangeSom.thetaD, synChangeSom.thetaP,nonlinear=nonlinear)


#############################################################
fig_width = 10  # width in inches
fig_height = 12  # height in inches
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

gs = gridspec.GridSpec(2, 2,  # ,
                       # width_ratios=[1,2],
                       height_ratios=[2, 1]
                       )

# fig.subplots_adjust(hspace=0.6)
# ax = plt.axes((0.3, 0.12, 0.65, 0.8))

# fig = plt.figure(figsize=(16,4))
fig.subplots_adjust(hspace=0.3, wspace=0.2)
plt.subplots_adjust(left=0.11, right=0.97, top=0.95, bottom=0.07)

# gssub0 = gridspec.GridSpecFromSubplotSpec(5, 1, subplot_spec=gs[0],hspace=0.2)


plt.figtext(0.0, 0.95, 'A',clip_on=False,color='black', weight='bold',size=24)
plt.figtext(0.52, 0.95, 'B',clip_on=False,color='black', weight='bold',size=24)
plt.figtext(0.0, 0.67, 'C',clip_on=False,color='black', weight='bold',size=22)
plt.figtext(0.52, 0.67, 'D',clip_on=False,color='black', weight='bold',size=22)
plt.figtext(0.0, 0.35, 'E',clip_on=False,color='black', weight='bold',size=22)
plt.figtext(0.52, 0.35, 'F',clip_on=False,color='black', weight='bold',size=22)
# plt.figtext(0.08, 0.36, 'E',clip_on=False,color='black', weight='bold',size=22)

###################################################################################################
gssub2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[0],hspace=0.2,wspace=0.3)

irrPairsFreq = np.load(dir0+'irregularSpikePairs_vs_rate_differentDeltaTs_STDdet_%s.npy' % parameterSetVis)


ax12 = plt.subplot(gssub2[0])
ax12.set_title('visual cortex')

ax12.axhline(y=1,ls='--',color='0.6',lw=2)
ax12.fill_between([0,25],0.95,1.7,facecolor='0.9')

ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,16]/w0,lw=2,c='red',label=r'$\Delta t = %d$ ms, $p=%s$' % (irrPairsFreq[0,4]*1000.,irrPairsFreq[0,11]))
ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,15]/w0,lw=2,c='red',alpha=0.6,label=r'$\Delta t = %d$ ms, $p=%s$' % (irrPairsFreq[0,3]*1000.,irrPairsFreq[0,10]))
ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,14]/w0,lw=2,c='k',label=r'$p=%s$' % (irrPairsFreq[0,9]))
ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,13]/w0,lw=2,c='blue',alpha=0.6,label=r'$\Delta t = %d$ ms, $p=%s$' % (irrPairsFreq[0,1]*1000.,irrPairsFreq[0,8]))
ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,12]/w0,lw=2,c='blue',label=r'$\Delta t = %d$ ms, $p=%s$' % (irrPairsFreq[0,0]*1000.,irrPairsFreq[0,7]))
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
ax12.set_ylim(0.9,1.7)
ax12.set_xlim(0,80)

plt.legend(loc=(0.35,0.1),frameon=False)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=12)

ax12.set_ylabel('change in synaptic strength')
#ax12.set_xlabel('rate (spk/sec)')

####################
########################################
ax1 = plt.subplot(gssub2[1])


ax1.axhline(y=0,color='0.4',ls='--')
ax1.fill_between([0,25],-0.06,0.16,facecolor='0.9')

ax1.plot(irrPairsFreq[:,5],(irrPairsFreq[:,16]-irrPairsFreq[:,14])/w0,color='r',lw=3,label=r'$p=0.4$, $\Delta = 10$ ms')
ax1.plot(irrPairsFreq[:,5],(irrPairsFreq[:,15]-irrPairsFreq[:,14])/w0,c='red',alpha=0.6,lw=3,label=r'$p=0.2$, $\Delta = 10$ ms')

ax1.plot(irrPairsFreq[:,5],(irrPairsFreq[:,12]-irrPairsFreq[:,14])/w0,color='b',lw=3,label=r'$p=0.2$, $\Delta = -10$ ms')
ax1.plot(irrPairsFreq[:,5],(irrPairsFreq[:,13]-irrPairsFreq[:,14])/w0,c='blue',alpha=0.6,lw=3,label=r'$p=0.4$, $\Delta = -10$ ms')

ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_position(('outward', 10))
ax1.spines['left'].set_position(('outward', 10))
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')
ax1.set_xlim(0,80)
ax1.set_ylim(-0.06,0.17)
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
ax1 = plt.subplot(gs[2])

NN = 1001

rateMatrix = zeros((NN, NN))
# pdb.set_trace()
for i in range(NN):
    rateMatrix[:, i] = linspace(1., 2. * 40., 2 * NN)[i:(i + NN)]

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
ax1.axvline(x=25, lw=2, c='0.5',ls='--')
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

ax1.text(33, 41, '     change in \nsynaptic strength',size=12)

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


###########################################################################################################
gssub2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1],hspace=0.2,wspace=0.3)

irrPairsFreq = np.load(dir0+'irregularSpikePairs_vs_rate_differentDeltaTs_STDdet_%s.npy' % parameterSetSom)
#irrPairsFreqS = np.load(dir0+'irregularSpikePairs_vs_rate_differentDeltaTs_STDstoch_%s.npy' % parameterSetSom)


ax12 = plt.subplot(gssub2[0])
ax12.set_title('somatosensory cortex')

ax12.fill_between([0,8],0.85,1.45,facecolor='0.9')
ax12.axhline(y=1,ls='--',color='0.6',lw=2)

ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,16]/w0,lw=2,c='red',label=r'$\Delta t = %s$ ms, $p=%s$' % (irrPairsFreq[0,4]*1000.,irrPairsFreq[0,11]))
#ax12.plot(irrPairsFreqS[:,5],irrPairsFreqS[:,16]/w0,lw=2,c='red',label=r'$\Delta t = %s$ ms, $\rho=%s$' % (irrPairsFreq[0,4]*1000.,irrPairsFreqS[0,11]))
ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,15]/w0,lw=2,c='red',alpha=0.6,label=r'$\Delta t = %s$ ms, $p=%s$' % (irrPairsFreq[0,3]*1000.,irrPairsFreq[0,10]))
ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,14]/w0,lw=2,c='k',label=r'$\rho=%s$' % (irrPairsFreq[0,9]))
#ax12.plot(irrPairsFreqS[:,5],irrPairsFreqS[:,14]/w0,lw=2,c='k',label=r'$\rho=%s$' % (irrPairsFreqS[0,9]))
ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,13]/w0,lw=2,c='blue',alpha=0.6,label=r'$\Delta t = %s$ ms, $p=%s$' % (irrPairsFreq[0,1]*1000.,irrPairsFreq[0,8]))
ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,12]/w0,lw=2,c='blue',label=r'$\Delta t = %s$ ms, $p=%s$' % (irrPairsFreq[0,0]*1000.,irrPairsFreq[0,7]))
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

ax12.set_xlim(0,20)
ax12.set_ylim(0.85,1.45)
#ax12.set_ylabel('change in synaptic strength')
#ax12.set_xlabel('rate (spk/sec)')

########################################
ax1 = plt.subplot(gssub2[1])

ax1.fill_between([0,8],-0.03,0.1,facecolor='0.9')
ax1.axhline(y=0,color='0.4',ls='--')

ax1.plot(irrPairsFreq[:,5],(irrPairsFreq[:,16]-irrPairsFreq[:,14])/w0,color='r',lw=3,label=r'$p=0.4$, $\Delta = 10$ ms')
ax1.plot(irrPairsFreq[:,5],(irrPairsFreq[:,15]-irrPairsFreq[:,14])/w0,c='red',alpha=0.6,lw=3,label=r'$p=0.2$, $\Delta = 10$ ms')

ax1.plot(irrPairsFreq[:,5],(irrPairsFreq[:,12]-irrPairsFreq[:,14])/w0,color='b',lw=3,label=r'$p=0.2$, $\Delta = -10$ ms')
ax1.plot(irrPairsFreq[:,5],(irrPairsFreq[:,13]-irrPairsFreq[:,14])/w0,c='blue',alpha=0.6,lw=3,label=r'$p=0.4$, $\Delta = -10$ ms')

ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_position(('outward', 10))
ax1.spines['left'].set_position(('outward', 10))
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')
ax1.set_xlim(0,20)
ax1.set_ylim(-0.06,0.16)
#plt.ylim(lls,uls)
#plt.xticks([-50, 0, 50],[-50, 0, 50])
#ax1.set_ylabel('sensitivity to correlations')
ax1.set_xlabel('firing rate (spk/s)')

#majorLocator_y = MaxNLocator(4)
ax1.yaxis.set_major_locator(MaxNLocator(5))
ax1.xaxis.set_major_locator(MaxNLocator(4))

############################
#########################################################################################################
# Calcium-model firing rate sensitivity
ax1 = plt.subplot(gs[3])

NN = 1001

rateMatrix = zeros((NN, NN))
# pdb.set_trace()
for i in range(NN):
    rateMatrix[:, i] = linspace(1., 2. * 10., 2 * NN)[i:(i + NN)]

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
ax1.axvline(x=8, lw=2, c='0.5',ls='--')
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

ax1.text(8.5, 10.3, '     change in \nsynaptic strength',size=12)

ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_position(('outward', 10))
ax1.spines['left'].set_position(('outward', 10))
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')
plt.xlim(0, 10)
plt.ylim(0, 10)
# plt.xticks([-50, 0, 50],[-50, 0, 50])
#plt.ylabel('increase in firing rate (spk/s)', multialignment='center')
plt.xlabel('baseline firing rate (spk/s)')

# majorLocator_y = MaxNLocator(4)
# ax1.yaxis.set_major_locator(majorLocator_y)
majorLocator_y = MaxNLocator(4)
ax1.yaxis.set_major_locator(majorLocator_y)
ax1.xaxis.set_major_locator(majorLocator_y)


###################################################################

fname = os.path.basename(__file__)

# savefig('ccgs_constant_coupled_networks.svg')
plt.savefig('outputFigures/' + fname[:-3] + '_v1.pdf')
plt.savefig('outputFigures/' + fname[:-3] + '_v1.png')  # plt.savefig(fname[:-3]+'_%s.png' % par)
# plt.savefig('delta_t.png')

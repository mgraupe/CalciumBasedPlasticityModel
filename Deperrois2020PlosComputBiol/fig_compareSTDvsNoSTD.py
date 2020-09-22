import numpy as  np
import scipy as sci
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
import pdb
import os
import matplotlib.ticker as ticker
import sys

sys.path.append("../")
#import synUtils as synUtils
import synapticChange as synChange

import params as par
import parameter_fit_solutions as goodSol


########################################################
# generate calcium trace from amplitudes
def generateCalciumTrace(tt,tListSorted,tauCa):
    catemp = np.zeros(len(tt))
    catempNoSTP = np.zeros(len(tt))
    ca = []
    # CaTotal, CaPre, CaPost, time
    ca.append([0., 0., 0., 0., 0., -1.])
    pre = 0
    post = 0
    preamps = []
    preampsNoSTD = []
    for i in tListSorted:
        #
        tMask = (tt >= ca[-1][5]) & (tt < i[0])
        caTotOld = ca[-1][0]
        caTotOldNoSTP = ca[-1][1]
        caPreOld = ca[-1][2]
        caPreOldNoSTP = ca[-1][3]
        caPostOld = ca[-1][4]
        tOld = ca[-1][5]
        # caTotTemp  = caTotOld*exp(-(i[0]-tOld)/self.tauCa)
        caPreTemp = caPreOld * np.exp(-(i[0] - tOld) / tauCa)
        caPreTempNoSTP = caPreOldNoSTP * np.exp(-(i[0] - tOld) / tauCa)
        caPostTemp = caPostOld * np.exp(-(i[0] - tOld) / tauCa)
        catemp[tMask] += caTotOld * np.exp(-(tt[tMask] - tt[tMask][0]) / tauCa)
        catempNoSTP[tMask] += caTotOldNoSTP * np.exp(-(tt[tMask] - tt[tMask][0]) / tauCa)
        # print caTotOld, tt[tMask],
        # postsynaptic spike
        if i[1] == 1:
            caPostTemp += i[2]
            post += 1
        # presynaptic spike
        if i[1] == 0:
            caPreTemp += i[2]
            caPreTempNoSTP += i[3]
            pre += 1
            preamps.append(caPreTemp + caPostTemp)
            preampsNoSTD.append(caPreTempNoSTP + caPostTemp)
        caTotTemp = caPreTemp + caPostTemp
        caTotTempNoSTP = caPreTempNoSTP + caPostTemp
        ca.append([caTotTemp, caTotTempNoSTP, caPreTemp, caPreTempNoSTP, caPostTemp, i[0]])

    return (catemp,catempNoSTP,preamps,preampsNoSTD)


####################################
def fitfunc(p,x,frequency):
    y   = np.zeros(len(x))
    y[0]= 1.
    for i in range(1,len(x)):
        y[i] = 1. - (1.-(y[i-1]-p[0]*y[i-1]))*np.exp(-1./(frequency*p[1]))
    return y

#####################################
parameterSetVis = 'sJFullSim0'
parameterSetVisNoSTD = 'sJFullNoSTDSim0'
parameterSetSom = 'sHFullSim1'
parameterSetSomNoSTD = 'sHFullNoSTDSim1b'

#def __init__(self, thetaD, nonlinear, Npairs, Npresentations, presentationInterval, w0, tauRec, U, dataOrigin = 'jesper', thetaP = None, Nves = None):
synUVis2 = synChange.synapticChange('sjoestroem',parameterSetVis,fromFile=True,nonlinear=1.,par=par)
synUSom2 = synChange.synapticChange('markram',parameterSetSom,fromFile=True,nonlinear=1.,par=par)
synUVisNoSTD2 = synChange.synapticChange('sjoestroem',parameterSetVisNoSTD,fromFile=True,USTD=0.,nonlinear=1.,par=par,thetaP=par.thetaPJ)
synUSomNoSTD2 = synChange.synapticChange('markram',parameterSetSomNoSTD,fromFile=True,USTD=0.,nonlinear=1.,par=par,thetaP=par.thetaPH)

figDir = 'outputFigures/'



exec('jesperSol = goodSol.%s' % parameterSetVis)
exec('jesperNoSTDSol = goodSol.%s' % parameterSetVisNoSTD) #goodSol.sJFullNoSTDSim0
exec('henrySol = goodSol.%s' % parameterSetSom) #goodSol.sHFullSim1
exec('henryNoSTDSol = goodSol.%s' % parameterSetSomNoSTD) #goodSol.sHFullNoSTDSim1

print('jesperSol[0] : ', jesperSol[0])
print('henrySol[0] :', henrySol[0])
print('jesperNoSTDSol[0] : ', jesperNoSTDSol[0])
print('henryNoSTDSol[0] :', henryNoSTDSol[0])

#pdb.set_trace()

####################################
# calculate solution
freq = np.linspace(1.5,50. ,500)
synChangeVis2 = np.zeros((len(freq) ,2))
synChangeVisNoSTD2 = np.zeros((len(freq) ,2))
synChangeSom2 = np.zeros((len(freq) ,2))
synChangeSomNoSTD2 = np.zeros((len(freq) ,2))

deltaTs = np.linspace(-0.05 ,0.05 ,501)
synChange2Vis = np.zeros((len(deltaTs),3))
synChange2Som = np.zeros((len(deltaTs),3))
freqVis = [0.1,10.,20.]
freqSom = [0.1,10.,20.]

for i in range(len(freq)):
    # calculateChangeInSynapticStrength(frequency,deltaT,params):
    # calculateChangeInSynapticStrengthSTP(self, frequency, deltaT, params, Npairs, stoch=False,DeltaTRange=None):
    # synUVis.calculateChangeInSynapticStrengthSTP(freq[i] ,0.01 ,jesperSol[0] ,synUVis.Npresentations ,pairCase='fullSimulation')
    synChangeVis2[i, 0] = synUVis2.calculateChangeInSynapticStrengthSTP(freq[i], 0.01, jesperSol[0], synUVis2.Npresentations, pairCase='fullSimulation')
    synChangeVis2[i, 1] = synUVis2.calculateChangeInSynapticStrengthSTP(freq[i], -0.01, jesperSol[0], synUVis2.Npresentations, pairCase='fullSimulation')
    synChangeVisNoSTD2[i, 0] = synUVisNoSTD2.calculateChangeInSynapticStrengthSTP(freq[i], 0.01, jesperNoSTDSol[0], synUVisNoSTD2.Npresentations, pairCase='fullSimulation')
    synChangeVisNoSTD2[i, 1] = synUVisNoSTD2.calculateChangeInSynapticStrengthSTP(freq[i], -0.01, jesperNoSTDSol[0], synUVisNoSTD2.Npresentations, pairCase='fullSimulation')

    synChangeSom2[i, 0] = synUSom2.calculateChangeInSynapticStrengthSTP(freq[i], 0.005, henrySol[0], synUSom2.Npresentations, pairCase='fullSimulation')
    synChangeSom2[i, 1] = synUSom2.calculateChangeInSynapticStrengthSTP(freq[i], -0.01, henrySol[0], synUSom2.Npresentations, pairCase='fullSimulation')
    synChangeSomNoSTD2[i, 0] = synUSomNoSTD2.calculateChangeInSynapticStrengthSTP(freq[i], 0.005, henryNoSTDSol[0], synUSomNoSTD2.Npresentations, pairCase='fullSimulation')
    synChangeSomNoSTD2[i, 1] = synUSomNoSTD2.calculateChangeInSynapticStrengthSTP(freq[i], -0.01, henryNoSTDSol[0], synUSomNoSTD2.Npresentations, pairCase='fullSimulation')

for i in range(len(deltaTs)):
    synChange2Vis[i,0] = synUVis2.calculateChangeInSynapticStrengthSTP(freqVis[0] ,deltaTs[i] ,jesperSol[0] ,synUVis2.Npresentations ,pairCase='fullSimulation')
    synChange2Vis[i,1] = synUVis2.calculateChangeInSynapticStrengthSTP(freqVis[1] ,deltaTs[i] ,jesperSol[0] ,synUVis2.Npresentations ,pairCase='fullSimulation')
    synChange2Vis[i,2] = synUVis2.calculateChangeInSynapticStrengthSTP(freqVis[2] ,deltaTs[i] ,jesperSol[0] ,synUVis2.Npresentations ,pairCase='fullSimulation')
    synChange2Som[i,0] = synUSom2.calculateChangeInSynapticStrengthSTP(freqSom[0], deltaTs[i], henrySol[0], synUSom2.Npresentations, pairCase='fullSimulation')
    synChange2Som[i,1] = synUSom2.calculateChangeInSynapticStrengthSTP(freqSom[1], deltaTs[i], henrySol[0], synUSom2.Npresentations, pairCase='fullSimulation')
    synChange2Som[i,2] = synUSom2.calculateChangeInSynapticStrengthSTP(freqSom[2], deltaTs[i], henrySol[0], synUSom2.Npresentations, pairCase='fullSimulation')

##################################################################
# location of irregular spike-pair simulations

nonlinear = 1.
w0 = 0.5



#dS     = 'jesper'

dir0 = 'simResults/'


###################################################################
# figure generation
fig_width = 10  # width in inches
fig_height = 11  # height in inches
fig_size = [fig_width, fig_height]
params = {'axes.labelsize': 14, 'axes.titlesize': 13, 'font.size': 13, 'xtick.labelsize': 13, 'ytick.labelsize': 13, 'figure.figsize': fig_size,  'axes.linewidth': 1.3,
          'ytick.major.size': 4,  # major tick size in points
          'xtick.major.size': 4  # major tick size in points
          # 'edgecolor' : None
          # 'xtick.major.size' : 2,
          # 'ytick.major.size' : 2,
          }
rcParams.update(params)

# set sans-serif font to Arial
#rcParams['font.sans-serif'] = 'Arial'

# create figure instance
fig = plt.figure()

# define sub-panel grid and possibly width and height ratios
gs = gridspec.GridSpec(3, 2,  # width_ratios=[1,1.2],
                       # height_ratios=[1,1]
                       )

# define vertical and horizontal spacing between panels
gs.update(wspace=0.3, hspace=0.3)

# possibly change outer margins of the figure
plt.subplots_adjust(left=0.12, right=0.96, top=0.95, bottom=0.06)

plt.figtext(0.018, 0.96, 'A',clip_on=False,color='black', weight='bold',size=22)
plt.figtext(0.52, 0.96, 'B',clip_on=False,color='black', weight='bold',size=22)
plt.figtext(0.018, 0.63, 'C',clip_on=False,color='black', weight='bold',size=22)
plt.figtext(0.52, 0.63, 'D',clip_on=False,color='black', weight='bold',size=22)
plt.figtext(0.018, 0.31, 'E',clip_on=False,color='black', weight='bold',size=22)
plt.figtext(0.52, 0.31, 'F',clip_on=False,color='black', weight='bold',size=22)


# first sub-plot #######################################################
ax0 = plt.subplot(gs[0])
ax0.set_title('visual cortex')

# diplay of data
ax0.axhline(y=1., c='0.7',ls='--')
plusMask = synUVis2.xDataReg[:, 1] > 0.
minusMask = synUVis2.xDataReg[:, 1] < 0.
ax0.errorbar(synUVis2.xDataReg[:, 0][plusMask], synUVis2.yDataReg[plusMask], yerr=synUVis2.sigmaDataReg[plusMask], fmt='s', color='C4', clip_on=False,label='data : $\Delta t = 10$ ms')
ax0.errorbar(synUVis2.xDataReg[:, 0][minusMask], synUVis2.yDataReg[minusMask], yerr=synUVis2.sigmaDataReg[minusMask], fmt='o', color='C5', clip_on=False,label='data : $\Delta t = -10$ ms')
ax0.plot(freq, synChangeVis2[:, 0], color='0.3')
ax0.plot(freq, synChangeVis2[:, 1], color='0.7')
ax0.plot(freq, synChangeVisNoSTD2[:, 0], color='C4',label='model : $\Delta t = 10$ ms')
ax0.plot(freq, synChangeVisNoSTD2[:, 1], color='C5',label='model : $\Delta t = -10$ ms')

# removes upper and right axes
# and moves left and bottom axes away
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.spines['bottom'].set_position(('outward', 10))
ax0.spines['left'].set_position(('outward', 10))
ax0.yaxis.set_ticks_position('left')
ax0.xaxis.set_ticks_position('bottom')

ax0.yaxis.set_major_locator(ticker.MaxNLocator(4))

# legends and labels
plt.legend(loc=(0.5,0.), frameon=False)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=10)

plt.xlabel('frequency (Hz)')
plt.ylabel('change in synaptic strength',position=(0.5,-0.17))

# first sub-plot #######################################################
ax12 = plt.subplot(gs[2])


irrPairsFreq = np.load(dir0+'irregularSpikePairs_vs_rate_differentDeltaTs_STDdet_%s.npy' % parameterSetVis)
irrPairsFreqNoSTD = np.load(dir0+'irregularSpikePairs_vs_rate_differentDeltaTs_STDdet_%s.npy' % parameterSetVisNoSTD)

#ax12 = plt.subplot(gssub2[0])
#ax12.set_title('visual cortex')

ax12.axhline(y=1,ls='--',color='0.6',lw=2)
ax12.fill_between([0,25],0.85,1.7,facecolor='0.9')

ax12.plot(irrPairsFreqNoSTD[:,5],irrPairsFreqNoSTD[:,16]/w0,lw=2,c='red',label=r'$\Delta t = %d$ ms, $p=%s$, w/o STD' % (irrPairsFreqNoSTD[0,4]*1000.,irrPairsFreqNoSTD[0,11]))
#ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,15]/w0,lw=2,c='red',alpha=0.6,label=r'$\Delta t = %d$ ms, $p=%s$' % (irrPairsFreq[0,3]*1000.,irrPairsFreq[0,10]))
ax12.plot(irrPairsFreqNoSTD[:,5],irrPairsFreqNoSTD[:,14]/w0,lw=2,c='darkred',label=r'$p=%s$, w/o STD' % (irrPairsFreqNoSTD[0,9]))

ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,16]/w0,lw=2,c='0.4',label=r'$\Delta t = %d$ ms, $p=%s$, with STD' % (irrPairsFreqNoSTD[0,4]*1000.,irrPairsFreqNoSTD[0,11]))#,label=r'$\Delta t = %d$ ms, $p=%s$' % (irrPairsFreq[0,4]*1000.,irrPairsFreq[0,11]))
#ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,15]/w0,lw=2,c='red',alpha=0.6,label=r'$\Delta t = %d$ ms, $p=%s$' % (irrPairsFreq[0,3]*1000.,irrPairsFreq[0,10]))
ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,14]/w0,lw=2,c='black',label=r'$p=0$, with STD')#,label=r'$p=%s$' % (irrPairsFreq[0,9]))

#ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,13]/w0,lw=2,c='blue',alpha=0.6,label=r'$\Delta t = %d$ ms, $p=%s$' % (irrPairsFreq[0,1]*1000.,irrPairsFreq[0,8]))
#ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,12]/w0,lw=2,c='blue',label=r'$\Delta t = %d$ ms, $p=%s$' % (irrPairsFreq[0,0]*1000.,irrPairsFreq[0,7]))

# removes upper and right axes
# and moves left and bottom axes away
ax12.spines['top'].set_visible(False)
ax12.spines['right'].set_visible(False)
ax12.spines['bottom'].set_position(('outward', 10))
ax12.spines['left'].set_position(('outward', 10))
ax12.yaxis.set_ticks_position('left')
ax12.xaxis.set_ticks_position('bottom')

ax12.set_ylim(0.85,1.7)
ax12.set_xlim(0,80)
# legends and labels
plt.legend(loc=(0.3,0.18), frameon=False)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=11)
#plt.xlabel(r'$\Delta t$ (ms)')
plt.xlabel('firing rate (spk/s)')

# first sub-plot #######################################################
ax14 = plt.subplot(gs[4])

ax14.fill_between([0,25],-0.03,0.24,facecolor='0.9')

# title
# ax1.set_title('regular Sjoestroem')
ax14.axhline(y=0,ls='--',color='0.4',lw=2)
ax14.plot(irrPairsFreqNoSTD[:,5],(irrPairsFreqNoSTD[:,16]-irrPairsFreqNoSTD[:,14])/w0,color='r',lw=3,label='w/o STD')#,label=r'$p=0.4$, $\Delta = 10$ ms')
ax14.plot(irrPairsFreq[:,5],(irrPairsFreq[:,16]-irrPairsFreq[:,14])/w0,color='0.4',lw=3,label='with STD')#,label=r'$p=0.4$, $\Delta = 10$ ms')
#ax14.plot(irrPairsFreq[:,5],(irrPairsFreq[:,15]-irrPairsFreq[:,14])/w0,c='red',alpha=0.6,lw=3,label=r'$p=0.2$, $\Delta = 10$ ms')

#ax14.plot(irrPairsFreq[:,5],(irrPairsFreq[:,12]-irrPairsFreq[:,14])/w0,color='b',lw=3,label=r'$p=0.2$, $\Delta = -10$ ms')
#ax14.plot(irrPairsFreq[:,5],(irrPairsFreq[:,13]-irrPairsFreq[:,14])/w0,c='blue',alpha=0.6,lw=3,label=r'$p=0.4$, $\Delta = -10$ ms')

# removes upper and right axes
# and moves left and bottom axes away
ax14.spines['top'].set_visible(False)
ax14.spines['right'].set_visible(False)
ax14.spines['bottom'].set_position(('outward', 10))
ax14.spines['left'].set_position(('outward', 10))
ax14.yaxis.set_ticks_position('left')
ax14.xaxis.set_ticks_position('bottom')

ax14.set_xlim(0,80)
ax14.set_ylim(-0.03,0.24)
ax14.yaxis.set_major_locator(ticker.MaxNLocator(4))
#ax1.set_xlim(-100,600)
# legends and labels
plt.legend(loc=1, frameon=False)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=11)

plt.xlabel(r'time (ms)')
#plt.ylabel('calcium')
plt.ylabel('sensitivity to correlations')

# first sub-plot #######################################################
ax0 = plt.subplot(gs[1])
ax0.set_title('somatosensory cortex')

# diplay of data
ax0.axhline(y=1., c='0.7',ls='--')
plusMask = synUSom2.xDataReg[:, 1] > 0.
minusMask = synUSom2.xDataReg[:, 1] < 0.
ax0.errorbar(synUSom2.xDataReg[:, 0][plusMask], synUSom2.yDataReg[plusMask], yerr=synUSom2.sigmaDataReg[plusMask], fmt='s', color='C4', clip_on=False,label='data : $\Delta t = 5$ ms')
ax0.errorbar(synUSom2.xDataReg[:, 0][minusMask], synUSom2.yDataReg[minusMask], yerr=synUSom2.sigmaDataReg[minusMask], fmt='o', color='C5', clip_on=False,label='data : $\Delta t = -10$ ms')
ax0.plot(freq, synChangeSom2[:, 0], color='0.3')
ax0.plot(freq, synChangeSom2[:, 1], color='0.7')
#ax0.plot(freq, synChangeSom2[:, 0],)
#ax0.plot(freq, synChangeSom2[:, 1],)
ax0.plot(freq, synChangeSomNoSTD2[:, 0], color='C4',label='model : $\Delta t = 5$ ms')
ax0.plot(freq, synChangeSomNoSTD2[:, 1], color='C5',label='model : $\Delta t = -10$ ms')
#ax0.plot(freq, synChangeSomNoSTD2[:, 0], color='C6',label='model : $\Delta t = 5$ ms')
#ax0.plot(freq, synChangeSomNoSTD2[:, 1], color='C7',label='model : $\Delta t = -10$ ms')

# removes upper and right axes
# and moves left and bottom axes away
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.spines['bottom'].set_position(('outward', 10))
ax0.spines['left'].set_position(('outward', 10))
ax0.yaxis.set_ticks_position('left')
ax0.xaxis.set_ticks_position('bottom')

ax0.yaxis.set_major_locator(ticker.MaxNLocator(4))
# legends and labels
plt.legend(loc=4, frameon=False)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=10)

plt.xlabel('frequency (Hz)')
#plt.ylabel('change in synaptic strength')

# first sub-plot #######################################################
ax13 = plt.subplot(gs[3])

irrPairsFreqSC = np.load(dir0+'irregularSpikePairs_vs_rate_differentDeltaTs_STDdet_%s.npy' % parameterSetSom)
irrPairsFreqNoSTDSC = np.load(dir0+'irregularSpikePairs_vs_rate_differentDeltaTs_STDdet_%s.npy' % parameterSetSomNoSTD)

# title
# ax1.set_title('regular Sjoestroem')
ax13.fill_between([0,8],0.81,1.45,facecolor='0.9')
ax13.axhline(y=1,color='0.4',ls='--')

ax13.plot(irrPairsFreqNoSTDSC[:,5],irrPairsFreqNoSTDSC[:,16]/w0,lw=2,c='red',label=r'$\Delta t = %d$ ms, $p=%s$, w/o STD' % (irrPairsFreqNoSTDSC[0,4]*1000.,irrPairsFreqNoSTD[0,11]))
#ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,15]/w0,lw=2,c='red',alpha=0.6,label=r'$\Delta t = %d$ ms, $p=%s$' % (irrPairsFreq[0,3]*1000.,irrPairsFreq[0,10]))
ax13.plot(irrPairsFreqNoSTDSC[:,5],irrPairsFreqNoSTDSC[:,14]/w0,lw=2,c='darkred',label=r'$p=%s$, w/o STD' % (irrPairsFreqNoSTDSC[0,9]))

ax13.plot(irrPairsFreqSC[:,5],irrPairsFreqSC[:,16]/w0,lw=2,c='0.4',label=r'$\Delta t = %d$ ms, $p=%s$, with STD' % (irrPairsFreqNoSTDSC[0,4]*1000.,irrPairsFreqNoSTD[0,11]))#,label=r'$\Delta t = %d$ ms, $p=%s$' % (irrPairsFreq[0,4]*1000.,irrPairsFreq[0,11]))
#ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,15]/w0,lw=2,c='red',alpha=0.6,label=r'$\Delta t = %d$ ms, $p=%s$' % (irrPairsFreq[0,3]*1000.,irrPairsFreq[0,10]))
ax13.plot(irrPairsFreqSC[:,5],irrPairsFreqSC[:,14]/w0,lw=2,c='black',label=r'$p=%s$, with STD' % (irrPairsFreqNoSTDSC[0,9]))#,label=r'$p=%s$' % (irrPairsFreq[0,9]))

#ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,13]/w0,lw=2,c='blue',alpha=0.6,label=r'$\Delta t = %d$ ms, $p=%s$' % (irrPairsFreq[0,1]*1000.,irrPairsFreq[0,8]))
#ax12.plot(irrPairsFreq[:,5],irrPairsFreq[:,12]/w0,lw=2,c='blue',label=r'$\Delta t = %d$ ms, $p=%s$' % (irrPairsFreq[0,0]*1000.,irrPairsFreq[0,7]))


# removes upper and right axes
# and moves left and bottom axes away
ax13.spines['top'].set_visible(False)
ax13.spines['right'].set_visible(False)
ax13.spines['bottom'].set_position(('outward', 10))
ax13.spines['left'].set_position(('outward', 10))
ax13.yaxis.set_ticks_position('left')
ax13.xaxis.set_ticks_position('bottom')
ax13.set_xlim(0,20)
ax13.set_ylim(0.81,1.45)
ax13.yaxis.set_major_locator(ticker.MaxNLocator(4))
# legends and labels
plt.legend(loc=(0.35,0.18), frameon=False)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=11)

plt.xlabel(r'firing rate (spk/s)')
#plt.ylabel('change in synaptic strength')

# first sub-plot #######################################################
ax15 = plt.subplot(gs[5])

ax15.fill_between([0,8],-0.015,0.13,facecolor='0.9')
ax15.axhline(y=0,color='0.4',ls='--')
# title
# ax1.set_title('regular Sjoestroem')
ax15.plot(irrPairsFreqNoSTDSC[:,5],(irrPairsFreqNoSTDSC[:,16]-irrPairsFreqNoSTDSC[:,14])/w0,color='r',lw=3,label='w/o STD')#,label=r'$p=0.4$, $\Delta = 10$ ms')
ax15.plot(irrPairsFreqSC[:,5],(irrPairsFreqSC[:,16]-irrPairsFreqSC[:,14])/w0,color='0.4',lw=3,label=r'with STD')


# removes upper and right axes
# and moves left and bottom axes away
ax15.spines['top'].set_visible(False)
ax15.spines['right'].set_visible(False)
ax15.spines['bottom'].set_position(('outward', 10))
ax15.spines['left'].set_position(('outward', 10))
ax15.yaxis.set_ticks_position('left')
ax15.xaxis.set_ticks_position('bottom')

ax15.yaxis.set_major_locator(ticker.MaxNLocator(4))
ax15.set_xlim(0,20)
ax15.set_ylim(-0.015,0.13)
# legends and labels
plt.legend(loc=1, frameon=False)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=11)

plt.xlabel(r'firing rate (spk/s)')


## save figure ############################################################
## save figure ############################################################
fname = os.path.basename(__file__)


plt.savefig(figDir + fname[:-3] + '_v1.png')
plt.savefig(figDir + fname[:-3] + '_v1.pdf')

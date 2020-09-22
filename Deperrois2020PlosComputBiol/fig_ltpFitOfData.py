import numpy as  np
import scipy as sci
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
import pdb
import os
import sys

sys.path.append("../")
#from synUtils import *
import synapticChange as synChange

import params as par
import parameter_fit_solutions as goodSol
import matplotlib.ticker as ticker


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
        #pdb.set_trace()
        tMask = (tt >= ca[-1][5]) & (tt < i[0])
        #pdb.set_trace()
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
#def __init__(self, thetaD, nonlinear, Npairs, Npresentations, presentationInterval, w0, tauRec, U, dataOrigin = 'jesper', thetaP = None, Nves = None):
parameterSetVis = 'sJFullSim0'
parameterSetSom = 'sHFullSim1'


synUVis2 = synChange.synapticChange('sjoestroem',parameterSetVis,fromFile=True,nonlinear=1.,par=par)
synUSom2 = synChange.synapticChange('markram',parameterSetSom,fromFile=True,nonlinear=1.,par=par)

figDir = 'outputFigures/'

exec('jesperSol = goodSol.%s' % parameterSetVis)
exec('henrySol = goodSol.%s' % parameterSetSom) #goodSol.sHFullSim1

#jesperSol = goodSol.sJFullSim0
#henrySol = goodSol.sHFullSim1

# choice of parameter set


####################################
# calculate solution
freq = np.linspace(1.3 ,50. ,500)
synChangeVis2 = np.zeros((len(freq) ,2))
synChangeSom2 = np.zeros((len(freq) ,2))

deltaTs = np.linspace(-0.05 ,0.05 ,501)
synChange2Vis2 = np.zeros((len(deltaTs),3))
synChange2Som2 = np.zeros((len(deltaTs),3))
freqVis = [2.,10.,20.]
freqSom = [2.,10.,20.]

for i in range(len(freq)):
    # calculateChangeInSynapticStrength(frequency,deltaT,params):
    # calculateChangeInSynapticStrengthSTP(self, frequency, deltaT, params, Npairs, stoch=False,DeltaTRange=None):
    # self.calculateChangeInSynapticStrengthSTP(freq[i],self.deltaTPlus,paraOpt[0],Npresentations,pairCase='fullSimulation')
    synChangeVis2[i, 0] = synUVis2.calculateChangeInSynapticStrengthSTP(freq[i], 0.01, jesperSol[0], synUVis2.Npresentations, pairCase='fullSimulation')
    synChangeVis2[i, 1] = synUVis2.calculateChangeInSynapticStrengthSTP(freq[i], -0.01, jesperSol[0], synUVis2.Npresentations, pairCase='fullSimulation')

    synChangeSom2[i, 0] = synUSom2.calculateChangeInSynapticStrengthSTP(freq[i], 0.005, henrySol[0], synUVis2.Npresentations, pairCase='fullSimulation')
    synChangeSom2[i, 1] = synUSom2.calculateChangeInSynapticStrengthSTP(freq[i], -0.01, henrySol[0], synUVis2.Npresentations, pairCase='fullSimulation')

for i in range(len(deltaTs)):
    synChange2Vis2[i, 0] = synUVis2.calculateChangeInSynapticStrengthSTP(freqVis[0], deltaTs[i], jesperSol[0], synUVis2.Npresentations, pairCase='fullSimulation')
    synChange2Vis2[i, 1] = synUVis2.calculateChangeInSynapticStrengthSTP(freqVis[1], deltaTs[i], jesperSol[0], synUVis2.Npresentations, pairCase='fullSimulation')
    synChange2Vis2[i, 2] = synUVis2.calculateChangeInSynapticStrengthSTP(freqVis[2], deltaTs[i], jesperSol[0], synUVis2.Npresentations, pairCase='fullSimulation')
    synChange2Som2[i,0] = synUSom2.calculateChangeInSynapticStrengthSTP(freqSom[0], deltaTs[i], henrySol[0] ,synUVis2.Npresentations, pairCase='fullSimulation')
    synChange2Som2[i,1] = synUSom2.calculateChangeInSynapticStrengthSTP(freqSom[1], deltaTs[i], henrySol[0] ,synUVis2.Npresentations, pairCase='fullSimulation')
    synChange2Som2[i,2] = synUSom2.calculateChangeInSynapticStrengthSTP(freqSom[2], deltaTs[i], henrySol[0], synUVis2.Npresentations, pairCase='fullSimulation')
##################################################
# for different frequencies generate example calcium traces for vis. and som. STD and no STD
# to get amplitudes for spike-pair stimulation

#############################################################################

# STD parameters for somatosensory cortex
Usom = 0.46
tauRecSom = 0.525
# STD parameters for visual cortex : from fit of deterministic model to the data
Uvis = 0.38375319
tauRecVis = 0.14891922

tauCaSom = henrySol[0][0]
CpreSom0  = henrySol[0][1]
CpostSom = henrySol[0][2]
thetaDSom = 1.
thetaPSom = henrySol[0][3]
DSom = henrySol[0][7]

tauCaVis = jesperSol[0][0]
CpreVis0  = jesperSol[0][1]
CpostVis = jesperSol[0][2]
thetaDVis = 1.
thetaPVis = jesperSol[0][3]
DVis = jesperSol[0][7]


freqs = [10.,22.]

calciumVis = []
calciumSom = []
#ampsNoSTD = []

for i in range(len(freqs)) :
    Npairs = 5
    tStart = 0.1  # start time at 100 ms
    deltaT = -0.01
    dt = 0.00001

    frequency = freqs[i]
    print(frequency)
    tPre = np.arange(Npairs) / frequency + tStart
    tPost = tPre + deltaT

    # somatosensory cortex
    tAll = np.zeros((2 * Npairs, 4))

    cpreSom = fitfunc([Usom,tauRecSom],np.arange(Npairs)/frequency,frequency)
    cpreSom *= Usom*CpreSom0
    tAll[:, 0] = np.hstack((tPre-DSom, tPost))
    tAll[:, 1] = np.hstack((np.zeros(Npairs), np.ones(Npairs)))
    tAll[:, 2] = np.hstack((cpreSom, np.repeat(CpostSom, Npairs)))
    tAll[:, 3] = np.hstack((np.repeat(CpreSom0 * Usom, Npairs), np.repeat(CpostSom, Npairs)))

    tList = tAll.tolist()
    tListSortedSom = sorted(tList, key=lambda tList: tList[0])

    tListSortedSom.append([tListSortedSom[-1][0] + 0.5, -1, 0., 0.])
    tEnd = tListSortedSom[-1][0] + 10. / frequency
    tt = np.linspace(0, tEnd, int(tEnd / dt) + 1)

    ( catempSomExamp,_ ,_,_) = generateCalciumTrace(tt, tListSortedSom, tauCaSom)
    calciumSom.append([tt,catempSomExamp])


    # visual cortex #######################
    tAll = np.zeros((2 * Npairs, 4))

    cpreVis = fitfunc([Uvis,tauRecVis],np.arange(Npairs)/frequency,frequency)
    cpreVis *= Uvis*CpreVis0
    tAll[:, 0] = np.hstack((tPre-DVis+dt, tPost))
    tAll[:, 1] = np.hstack((np.zeros(Npairs), np.ones(Npairs)))
    tAll[:, 2] = np.hstack((cpreVis, np.repeat(CpostVis, Npairs)))
    tAll[:, 3] = np.hstack((np.repeat(CpreVis0 * Uvis, Npairs), np.repeat(CpostVis, Npairs)))

    tList = tAll.tolist()
    tListSortedVis = sorted(tList, key=lambda tList: tList[0])

    tListSortedVis.append([tListSortedVis[-1][0] + 0.5, -1, 0., 0.])
    tEnd = tListSortedVis[-1][0] + 0.5
    tt = np.linspace(0, tEnd, int(tEnd / dt) + 1)

    (catempVisExamp,_, _,_) = generateCalciumTrace(tt, tListSortedVis, tauCaVis)
    calciumVis.append([tt,catempVisExamp])




###################################################################
# figure generation
fig_width = 10  # width in inches
fig_height = 11  # height in inches
fig_size = [fig_width, fig_height]
params = {'axes.labelsize': 14, 'axes.titlesize': 13, 'font.size': 13, 'xtick.labelsize': 13, 'ytick.labelsize': 13, 'figure.figsize': fig_size, 'axes.linewidth': 1.3,
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

plt.figtext(0.01, 0.96, 'A',clip_on=False,color='black', weight='bold',size=22)
plt.figtext(0.52, 0.96, 'B',clip_on=False,color='black', weight='bold',size=22)
plt.figtext(0.01, 0.63, 'C',clip_on=False,color='black', weight='bold',size=22)
plt.figtext(0.52, 0.63, 'D',clip_on=False,color='black', weight='bold',size=22)
plt.figtext(0.01, 0.31, 'E',clip_on=False,color='black', weight='bold',size=22)
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
ax0.plot(freq, synChangeVis2[:, 0], color='C4',label='model : $\Delta t = 10$ ms')
ax0.plot(freq, synChangeVis2[:, 1], color='C5',label='model : $\Delta t = -10$ ms')

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
ax1 = plt.subplot(gs[2])

# title
# ax1.set_title('regular Sjoestroem')

ax1.axhline(y=1., c='0.7',ls='--')
ax1.axvline(x=0., c='0.7',ls='--')
dMask1 = (deltaTs> -1./(2.*freqVis[1]))&(deltaTs < 1./(2.*freqVis[1]))
dMask2 = (deltaTs> -1./(2.*freqVis[2]))&(deltaTs < 1./(2.*freqVis[2]))
ax1.plot(deltaTs * 1000., synChange2Vis2[:,0],label=r'$f=%s$ Hz'%int(freqVis[0]),lw=2)
ax1.plot(deltaTs[dMask1] * 1000., synChange2Vis2[:,1][dMask1],label=r'$f=%s$ Hz'% int(freqVis[1]),lw=2)
ax1.plot(deltaTs[dMask2] * 1000., synChange2Vis2[:,2][dMask2],label=r'$f=%s$ Hz'% int(freqVis[2]),lw=2)

# removes upper and right axes
# and moves left and bottom axes away
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_position(('outward', 10))
ax1.spines['left'].set_position(('outward', 10))
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')
#ax1.yaxis.set_major_locator(ticker.MaxNLocator(4))

# legends and labels
plt.legend(loc=(0.65,0.78), frameon=False)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=10)
plt.xlabel(r'$\Delta t$ (ms)')

# first sub-plot #######################################################
ax1 = plt.subplot(gs[4])

# title
# ax1.set_title('regular Sjoestroem')

ax1.axhline(y=thetaDVis, c='0.7',ls='--')
ax1.axhline(y=thetaPVis, c='0.7',ls='-')
#dMask2 = (deltaTs> -1./(2.*freqVis[2]))&(deltaTs < 1./(2.*freqVis[2]))
#ax1.plot(deltaTs * 1000., synChange2Som[:,0],label=r'$f=%s$ Hz' % freqSom[0],lw=2)
#ax1.plot(deltaTs * 1000., synChange2Som[:,1],label=r'$f=%s$ Hz' % int(freqSom[1]),lw=2)
ax1.plot(calciumVis[0][0] * 1000.-tStart*1000., calciumVis[0][1], label=r'$f=%s$ Hz' % (int(freqs[0])),lw=2,c='0.3')
ax1.plot(calciumVis[1][0] * 1000.-tStart*1000., calciumVis[1][1], label=r'$f=%s$ Hz' % (int(freqs[1])),lw=2,c='0.8')

# removes upper and right axes
# and moves left and bottom axes away
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_position(('outward', 10))
ax1.spines['left'].set_position(('outward', 10))
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')

ax1.yaxis.set_major_locator(ticker.MaxNLocator(4))

ax1.set_xlim(-100,600)
# legends and labels
plt.legend(loc=1, frameon=False)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=10)

plt.xlabel(r'time (ms)')
plt.ylabel('calcium')


# first sub-plot #######################################################
ax0 = plt.subplot(gs[1])
ax0.set_title('somatosensory cortex')

# diplay of data
ax0.axhline(y=1., c='0.7',ls='--')
plusMask = synUSom2.xDataReg[:, 1] > 0.
minusMask = synUSom2.xDataReg[:, 1] < 0.
ax0.errorbar(synUSom2.xDataReg[:, 0][plusMask], synUSom2.yDataReg[plusMask], yerr=synUSom2.sigmaDataReg[plusMask], fmt='s', color='C4', clip_on=False,label='data : $\Delta t = 5$ ms')
ax0.errorbar(synUSom2.xDataReg[:, 0][minusMask], synUSom2.yDataReg[minusMask], yerr=synUSom2.sigmaDataReg[minusMask], fmt='o', color='C5', clip_on=False,label='data : $\Delta t = -10$ ms')
ax0.plot(freq, synChangeSom2[:, 0], color='C4',label='model : $\Delta t = 5$ ms')
ax0.plot(freq, synChangeSom2[:, 1], color='C5',label='model : $\Delta t = -10$ ms')

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
ax1 = plt.subplot(gs[3])

# title
# ax1.set_title('regular Sjoestroem')

ax1.axhline(y=1., c='0.7',ls='--')
ax1.axvline(x=0., c='0.7',ls='--')
dMask2 = (deltaTs> -1./(2.*freqVis[2]))&(deltaTs < 1./(2.*freqVis[2]))
ax1.plot(deltaTs * 1000., synChange2Som2[:,0],label=r'$f=%s$ Hz' % int(freqSom[0]),lw=2)
ax1.plot(deltaTs * 1000., synChange2Som2[:,1],label=r'$f=%s$ Hz' % int(freqSom[1]),lw=2)
ax1.plot(deltaTs[dMask2] * 1000., synChange2Som2[:,2][dMask2],label=r'$f=%s$ Hz' % int(freqSom[2]),lw=2)

# removes upper and right axes
# and moves left and bottom axes away
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_position(('outward', 10))
ax1.spines['left'].set_position(('outward', 10))
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')

ax1.yaxis.set_major_locator(ticker.MaxNLocator(4))
# legends and labels
plt.legend(loc=4, frameon=False)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=10)

plt.xlabel(r'$\Delta t$ (ms)')
#plt.ylabel('change in synaptic strength')

# first sub-plot #######################################################
ax1 = plt.subplot(gs[5])

# title
# ax1.set_title('regular Sjoestroem')

ax1.axhline(y=thetaDSom, c='0.7',ls='--')
ax1.axhline(y=thetaPSom, c='0.7',ls='-')
#dMask2 = (deltaTs> -1./(2.*freqVis[2]))&(deltaTs < 1./(2.*freqVis[2]))
#ax1.plot(deltaTs * 1000., synChange2Som[:,0],label=r'$f=%s$ Hz' % freqSom[0],lw=2)
#ax1.plot(deltaTs * 1000., synChange2Som[:,1],label=r'$f=%s$ Hz' % int(freqSom[1]),lw=2)
ax1.plot(calciumSom[0][0] * 1000.-tStart*1000., calciumSom[0][1], label=r'$f=%s$ Hz' % (int(freqs[0])),lw=2,c='0.3')
ax1.plot(calciumSom[1][0] * 1000.-tStart*1000., calciumSom[1][1], label=r'$f=%s$ Hz' % (int(freqs[1])),lw=2,c='0.8')

# removes upper and right axes
# and moves left and bottom axes away
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_position(('outward', 10))
ax1.spines['left'].set_position(('outward', 10))
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')

ax1.yaxis.set_major_locator(ticker.MaxNLocator(4))
ax1.set_xlim(-100,600)
# legends and labels
plt.legend(loc=1, frameon=False)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=10)

plt.xlabel(r'$\Delta t$ (ms)')
#plt.ylabel('change in synaptic strength')

## save figure ############################################################
## save figure ############################################################
fname = os.path.basename(__file__)


plt.savefig(figDir + fname[:-3] + '_v2.png')
plt.savefig(figDir + fname[:-3] + '_v2.pdf')

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
import multiprocessing
import pdb

from timeAboveThreshold.timeAboveThreshold import timeAboveThreshold
from synapticChange import synapticChange


def runIrregularPairSimulations(args):
    dT       = args[0]
    preRate  = args[1]
    postRate = args[2]
    p        = args[3]
    
    (alphaD,alphaP) = tat.irregularSpikePairs(dT,preRate,postRate,p)
    synChange.changeInSynapticStrength(T_total,rho0,alphaD,alphaP)
    
    return synChange.mean

##########################################################
# output directory
outputDir = 'simResults/'

##########################################################
# synaptic change vs Delta T for irregular Pairs
##########################################################
# Parameter of the stimulation protocol
frequencies   = array([1.,5.,10.,20.,30.,50.])   # frequency of spike-pair presentations in pairs/sec
T_total     = 10.     # total time of stimulation in sec
DeltaTstart = -0.1    # start time difference between pre- and post-spike, in sec
DeltaTend   =  0.1    # end time difference between pre- and post-spike, in sec
DeltaTsteps =  101.  # steps between start and end value
ppp         = 1.
rho0        = 0.5

# nonlinearity factor
nl = 1.

# parameter-set to use
params = 'solOld'

nCases = len(frequencies)


###########################################################
# initiate synaptic change class and chose parameter set from file
synChange = synapticChange(params,fromFile=True,nonlinear=nl)
# initiate class to calculate fraction of time above threshold
tat = timeAboveThreshold(synChange.tauCa, synChange.Cpre, synChange.Cpost, synChange.thetaD, synChange.thetaP, nonlinear=nl)


###########################################################
# initialize arrays 

deltaT = linspace(DeltaTstart,DeltaTend,DeltaTsteps)

pool = multiprocessing.Pool()

resultsIrr = zeros(len(frequencies)*3+2)

###########################################################
# simulation loop over range of deltaT values
for i in range(len(deltaT)):
    #
    print 'deltaT : ', deltaT[i]
    
    args = column_stack((ones(nCases)*deltaT[i],frequencies,frequencies,ones(nCases)*ppp))

    rrr = pool.map(runIrregularPairSimulations,args)
    #for n in range(len(frequencies)):
    #    (synC[i,n],meanU[i,n],meanD[i,n],tD[i,n],tP[i,n]) = rrr[n]
    #pdb.set_trace()
    res1 = hstack((deltaT[i],frequencies,frequencies,ppp,rrr))
    resultsIrr = vstack((resultsIrr,res1))

resultsIrr = resultsIrr[1:]

# the array has to be flipped if Cpre>Cpost 
if synChange.Cpre>synChange.Cpost:
    # invert time axis
    resultsIrr[:,0] = resultsIrr[:,0][::-1]
    # reorder entire array to end up with increasing time
    resultsIrr = resultsIrr[::-1]
    

np.save(outputDir+'irregularSpikePairs_vs_deltaT_%s.npy' % params,resultsIrr)
np.savetxt(outputDir+'irregularSpikePairs_vs_deltaT_%s.dat' % params,resultsIrr)


##########################################################
# synaptic change vs frequency
##########################################################

deltaTs   = array([-0.01,0.,0.01])   # frequency of spike-pair presentations in pairs/sec
T_total     = 10.     # total time of stimulation in sec
Freqstart = 1.    # start time difference between pre- and post-spike, in sec
FreqTend   =  60.    # end time difference between pre- and post-spike, in sec
FreqSteps =  60  # steps between start and end value
ppp         = array([0.4,0.,0.4])
rho0        = 0.5

# nonlinearity factor
nl = 1.

# parameter-set to use
params = 'solOld'

nCases = len(deltaTs)


###########################################################
# initialize arrays 

frequencies = linspace(Freqstart,FreqTend,FreqSteps)

results = zeros(len(deltaTs)*3+2)

###########################################################
# simulation loop over range of deltaT values
for i in range(len(frequencies)):
    #
    print 'rate : ', frequencies[i]
    
    args = column_stack((deltaTs,ones(nCases)*frequencies[i],ones(nCases)*frequencies[i],ppp))

    rrr = pool.map(runIrregularPairSimulations,args)
    #for n in range(len(frequencies)):
    #    (synC[i,n],meanU[i,n],meanD[i,n],tD[i,n],tP[i,n]) = rrr[n]
    #pdb.set_trace()
    res1 = hstack((deltaTs,frequencies[i],frequencies[i],ppp,rrr))
    results = vstack((results,res1))

results = results[1:]

# the array has to be flipped if Cpre>Cpost 
if synChange.Cpre>synChange.Cpost:
    # invert time axis
    results[:,0] = results[:,0][::-1]
    # reorder entire array to end up with increasing time
    results = results[::-1]
    

np.save(outputDir+'irregularSpikePairs_vs_rate_%s.npy' % params,results)
np.savetxt(outputDir+'irregularSpikePairs_vs_rate_%s.dat' % params,results)

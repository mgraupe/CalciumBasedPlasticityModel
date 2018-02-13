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


##########################################################
def runIrregularPairSimulations(args):
    dT = args[0]
    preRate = args[1]
    postRate = args[2]
    p = args[3]

    if synChange.Cpre > synChange.Cpost:
        (alphaD, alphaP) = tat.irregularSpikePairs(dT + synChange.D, preRate, postRate, p, deltaCa)
    else:
        (alphaD, alphaP) = tat.irregularSpikePairs(dT - synChange.D, preRate, postRate, p, deltaCa)
    synChange.changeInSynapticStrength(T_total, rho0, alphaD, alphaP)

    return synChange.mean


##########################################################
def runIrregularPairSTPDeterministicSimulations(args):
    dT = args[0]
    preRate = args[1]
    postRate = args[2]
    p = args[3]

    (alphaD, alphaP) = tat.irregularSpikePairsSTPDeterministic(dT - synChange.D, preRate, postRate, p, synChange.tauRec, synChange.U)

    synChange.changeInSynapticStrength(T_total, rho0, alphaD, alphaP)

    return synChange.mean


##########################################################
def runIrregularPairSTPStochasticSimulations(args):
    dT = args[0]
    preRate = args[1]
    postRate = args[2]
    p = args[3]

    (alphaD, alphaP) = tat.irregularSpikePairsSTPStochastic(dT - synChange.D, preRate, postRate, p, synChange.tauRec,
                                                            synChange.U, synChange.Nvesicles)

    synChange.changeInSynapticStrength(T_total, rho0, alphaD, alphaP)

    return synChange.mean


##########################################################
def runRegularPairSimulations(args):
    dT = args[0]
    preRate = args[1]
    postRate = args[2]
    p = args[3]

    (alphaD, alphaP) = tat.spikePairFrequency(dT - synChange.D, preRate)
    synChange.changeInSynapticStrength(T_total, rho0, alphaD, alphaP)

    return synChange.mean


##########################################################
def runRegularPairSTPDeterministicSimulations(args):
    dT = args[0]
    preRate = args[1]
    postRate = args[2]
    p = args[3]

    # here alphaD and alphaP are actually the absolute times spent above threshold
    (alphaD, alphaP) = tat.spikePairFrequencySTPDeterministic(dT - synChange.D, preRate, synChange.Npairs,
                                                              synChange.tauRec, synChange.U)
    # in turn T_total turns into the number of the stimulation pattern presentation
    synChange.changeInSynapticStrength(synChange.Npresentations, rho0, alphaD, alphaP)

    return synChange.mean


##########################################################
def runRegularPairSTPStochasticSimulations(args):
    dT = args[0]
    preRate = args[1]
    postRate = args[2]
    p = args[3]

    # here alphaD and alphaP are actually the absolute times spent above threshold
    (alphaD, alphaP) = tat.spikePairFrequencySTPStochastic(dT - synChange.D, preRate, synChange.Npairs,
                                                           synChange.tauRec, synChange.U, synChange.Nvesicles)
    # in turn T_total turns into the number of the stimulation pattern presentation
    synChange.changeInSynapticStrength(synChange.Npresentations, rho0, alphaD, alphaP)

    return synChange.mean


##########################################################
##########################################################
# output directory
outputDir = 'simResults/'

# numerical integration step width
deltaCa     = 0.0001 #0.01 #  0.0001 
T_total     = 10.     # total time of stimulation in sec
rho0        = 0.5
nl = 1.  # nonlinearity factor

###########################################################
# initiate synaptic change class and chose parameter set from file
parameterSet = 'stpJesperCaModel'
#dS     = 'jesper'

synChange = synapticChange(parameterSet,fromFile=True,nonlinear=nl)
# initiate class to calculate fraction of time above threshold
tat = timeAboveThreshold(synChange.tauCa, synChange.Cpre, synChange.Cpost, synChange.thetaD, synChange.thetaP, nonlinear=nl)

pool = multiprocessing.Pool()

##################################################################################################
# synaptic change vs Delta T for irregular Pairs
#################################################################################################
print 'irregular pairs : synaptic change vs Delta T for frequencies p\'s, stochastic STD'

# Parameter of the stimulation protocol
frequencies   = array([1.,5.,10.,20.,40.,80.])   # frequency of spike-pair presentations in pairs/sec
DeltaTstart = -0.1    # start time difference between pre- and post-spike, in sec
DeltaTend   =  0.1    # end time difference between pre- and post-spike, in sec
DeltaTsteps =  101.  # steps between start and end value
ppp         = 1.

nCases = len(frequencies)


###########################################################
# initialize arrays 
deltaT = linspace(DeltaTstart,DeltaTend,DeltaTsteps)
resultsIrr = zeros(len(frequencies)*3+2)

###########################################################
# simulation loop over range of deltaT values
for i in range(len(deltaT)):
    #
    print 'deltaT : ', deltaT[i]
    
    args = column_stack((ones(nCases)*deltaT[i],frequencies,frequencies,ones(nCases)*ppp))

    rrr = pool.map(runIrregularPairSTPStochasticSimulations,args)
    #for n in range(len(frequencies)):
    #    (synC[i,n],meanU[i,n],meanD[i,n],tD[i,n],tP[i,n]) = rrr[n]
    #pdb.set_trace()
    res1 = hstack((deltaT[i],frequencies,frequencies,ppp,rrr))
    resultsIrr = vstack((resultsIrr,res1))

resultsIrr = resultsIrr[1:]

if not os.path.exists(outputDir):
    os.makedirs(outputDir)
    
np.save(outputDir+'irregularSpikePairs_vs_deltaT_differentFreqs_STDStoch_%s.npy' % parameterSet,resultsIrr)
np.savetxt(outputDir+'irregularSpikePairs_vs_deltaT_differentFreqs_STDStoch_%s.dat' % parameterSet,resultsIrr)

#pdb.set_trace()

##################################################################################################
# synaptic change vs Delta T for irregular Pairs
##################################################################################################
print 'irregular pairs : synaptic change vs Delta T for different p\'s, stochastic STD'

# Parameter of the stimulation protocol
frequency   = 10.   # frequency of spike-pair presentations in pairs/sec
DeltaTstart = -0.1    # start time difference between pre- and post-spike, in sec
DeltaTend   =  0.1    # end time difference between pre- and post-spike, in sec
DeltaTsteps =  101.  # steps between start and end value
ppp         = array([0.1,0.25,0.5,0.75,1.])

nCases = len(ppp)


###########################################################
# initialize arrays 
deltaT = linspace(DeltaTstart,DeltaTend,DeltaTsteps)
resultsIrr = zeros(nCases*2+3)

###########################################################
# simulation loop over range of deltaT values
for i in range(len(deltaT)):
    #
    print 'deltaT : ', deltaT[i]
    
    args = column_stack((ones(nCases)*deltaT[i],ones(nCases)*frequency,ones(nCases)*frequency,ppp))

    rrr = pool.map(runIrregularPairSTPStochasticSimulations,args)
    #for n in range(len(frequencies)):
    #    (synC[i,n],meanU[i,n],meanD[i,n],tD[i,n],tP[i,n]) = rrr[n]
    #pdb.set_trace()
    res1 = hstack((deltaT[i],frequency,frequency,ppp,rrr))
    resultsIrr = vstack((resultsIrr,res1))

resultsIrr = resultsIrr[1:]


if not os.path.exists(outputDir):
    os.makedirs(outputDir)
    
np.save(outputDir+'irregularSpikePairs_vs_deltaT_differentPs_STDStoch_%s.npy' % parameterSet,resultsIrr)
np.savetxt(outputDir+'irregularSpikePairs_vs_deltaT_differentPs_STDStoch_%s.dat' % parameterSet,resultsIrr)

##########################################################
# synaptic change vs Delta T for regular Pairs
##########################################################
print 'regular pairs : synaptic change vs Delta T for frequencies p\'s, stochastic STD'

# Parameter of the stimulation protocol
frequencies   = array([1.,5.,10.,20.,40.,80.])   # frequency of spike-pair presentations in pairs/sec
DeltaTstart = -0.1    # start time difference between pre- and post-spike, in sec
DeltaTend   =  0.1    # end time difference between pre- and post-spike, in sec
DeltaTsteps =  101.  # steps between start and end value
ppp         = 1.

nCases = len(frequencies)

###########################################################
# initialize arrays 
deltaT = linspace(DeltaTstart,DeltaTend,DeltaTsteps)
resultsReg = zeros(len(frequencies)*3+2)

###########################################################
# simulation loop over range of deltaT values
for i in range(len(deltaT)):
    #
    print 'deltaT : ', deltaT[i]
    
    args = column_stack((ones(nCases)*deltaT[i],frequencies,frequencies,ones(nCases)*ppp))

    rrr = pool.map(runRegularPairSTPStochasticSimulations,args)
    #for n in range(len(frequencies)):
    #    (synC[i,n],meanU[i,n],meanD[i,n],tD[i,n],tP[i,n]) = rrr[n]
    #pdb.set_trace()
    res1 = hstack((deltaT[i],frequencies,frequencies,ppp,rrr))
    resultsReg = vstack((resultsReg,res1))

resultsReg = resultsReg[1:]

np.save(outputDir+'regularSpikePairs_vs_deltaT_differentFreqs_STDStoch_%s.npy' % parameterSet,resultsReg)
np.savetxt(outputDir+'regularSpikePairs_vs_deltaT_differentFreqs_STDStoch_%s.dat' % parameterSet,resultsReg)

##################################################################################################
# synaptic change vs frequency
##################################################################################################
print 'irregular pairs : synaptic change vs rate for different deltaT\'s and p\'s, stochastic STD'

deltaTs   = array([-0.01,-0.01,0.,0.01,0.01])   # frequency of spike-pair presentations in pairs/sec
Freqstart = 0.1    # start time difference between pre- and post-spike, in sec
FreqTend   =  80.    # end time difference between pre- and post-spike, in sec
FreqSteps =  120  # steps between start and end value
ppp         = array([0.4,0.2,0.,0.2,0.4])

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

    rrr = pool.map(runIrregularPairSTPStochasticSimulations,args)
    #for n in range(len(frequencies)):
    #    (synC[i,n],meanU[i,n],meanD[i,n],tD[i,n],tP[i,n]) = rrr[n]
    #pdb.set_trace()
    res1 = hstack((deltaTs,frequencies[i],frequencies[i],ppp,rrr))
    results = vstack((results,res1))

results = results[1:]

np.save(outputDir+'irregularSpikePairs_vs_rate_differentDeltaTs_STDStoch_%s.npy' % parameterSet,results)
np.savetxt(outputDir+'irregularSpikePairs_vs_rate_differentDeltaTs_STDStoch_%s.dat' % parameterSet,results)


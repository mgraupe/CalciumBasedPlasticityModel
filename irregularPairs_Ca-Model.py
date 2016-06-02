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
# Parameter sets for the calcium and synaptic weight dynamics
class calculateChange():
    ###############################################################################
    def __init__(self, params):
        
        
            
        # output directory
        self.outputDir = 'simResults/'

        # numerical integration step width
        self.deltaCa = 0.005 #0.01 #  0.0001 

        # nonlinearity factor
        self.nl     = 1.
        self.rho0   = 0.5
        T_total     = 10.     # total time of stimulation in sec
        
        ###########################################################
        # initiate synaptic change class and chose parameter set from file
        self.synChange = synapticChange(params,fromFile=True,nonlinear=self.nl)
        # initiate class to calculate fraction of time above threshold
        self.tat = timeAboveThreshold(self.synChange.tauCa, self.synChange.Cpre, self.synChange.Cpost, self.synChange.thetaD, self.synChange.thetaP, nonlinear=self.nl)


    ##########################################################
    def runIrregularPairSimulations(self,args):
        dT       = args[0]
        preRate  = args[1]
        postRate = args[2]
        p        = args[3]
        
        if self.synChange.Cpre>self.synChange.Cpost:
            (alphaD,alphaP) = tat.irregularSpikePairs(dT+self.synChange.D,preRate,postRate,p,self.deltaCa)
        else:
            (alphaD,alphaP) = tat.irregularSpikePairs(dT-self.synChange.D,preRate,postRate,p,self.deltaCa)
        self.synChange.changeInSynapticStrength(self.T_total,self.rho0,alphaD,alphaP)
        
        return self.synChange.mean

    ##########################################################
    def runRegularPairSimulations(self,args):
        dT       = args[0]
        preRate  = args[1]
        postRate = args[2]
        p        = args[3]
        
        (alphaD,alphaP) = tat.spikePairFrequency(dT-self.synChange.D,preRate)
        self.synChange.changeInSynapticStrength(self.T_total,self.rho0,alphaD,alphaP)
        
        return self.synChange.mean

    ##################################################################################################
    # synaptic change vs Delta T for irregular Pairs
    ##################################################################################################
    def changeVsDeltaTirregular_differentFrequencies(self):
    
        # Parameter of the stimulation protocol
        frequencies   = array([1.,5.,10.,20.,30.,50.])   # frequency of spike-pair presentations in pairs/sec
        DeltaTstart = -0.1    # start time difference between pre- and post-spike, in sec
        DeltaTend   =  0.1    # end time difference between pre- and post-spike, in sec
        DeltaTsteps =  101.  # steps between start and end value
        ppp         = 1.

        nCases = len(frequencies)
        
        ##########################################################
        pool = multiprocessing.Pool(nCases)
        
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

            rrr = pool.map(self.runIrregularPairSimulations,args)
            #for n in range(len(frequencies)):
            #    (synC[i,n],meanU[i,n],meanD[i,n],tD[i,n],tP[i,n]) = rrr[n]
            #pdb.set_trace()
            res1 = hstack((deltaT[i],frequencies,frequencies,ppp,rrr))
            resultsIrr = vstack((resultsIrr,res1))

        resultsIrr = resultsIrr[1:]

        # the array has to be flipped if Cpre>Cpost 
        if self.synChange.Cpre>self.synChange.Cpost:
            # invert time axis
            resultsIrr[:,0] = resultsIrr[:,0][::-1] 
            # reorder entire array to end up with increasing time
            resultsIrr = resultsIrr[::-1]
            
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)
            
        np.save(outputDir+'irregularSpikePairs_vs_deltaT_differentFreq_%s.npy' % params,resultsIrr)
        np.savetxt(outputDir+'irregularSpikePairs_vs_deltaT_differentFreq_%s.dat' % params,resultsIrr)

    ##################################################################################################
    # synaptic change vs Delta T for regular Pairs
    ##################################################################################################
    def changeVsDeltaTregular_differentFrequencies(self):

        # Parameter of the stimulation protocol
        frequencies   = array([1.,5.,10.,20.,30.,50.])   # frequency of spike-pair presentations in pairs/sec
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

            rrr = pool.map(runRegularPairSimulations,args)
            #for n in range(len(frequencies)):
            #    (synC[i,n],meanU[i,n],meanD[i,n],tD[i,n],tP[i,n]) = rrr[n]
            #pdb.set_trace()
            res1 = hstack((deltaT[i],frequencies,frequencies,ppp,rrr))
            resultsReg = vstack((resultsReg,res1))

        resultsReg = resultsReg[1:]

        np.save(outputDir+'regularSpikePairs_vs_deltaT_%s.npy' % params,resultsReg)
        np.savetxt(outputDir+'regularSpikePairs_vs_deltaT_%s.dat' % params,resultsReg)

    ##################################################################################################
    # synaptic change vs Delta T for irregular Pairs
    ##################################################################################################
    def changeVsDeltaTirregular_differentPs(self):
    
        # Parameter of the stimulation protocol
        frequency   = 10.     # frequency of spike-pair presentations in pairs/sec
        DeltaTstart = -0.1    # start time difference between pre- and post-spike, in sec
        DeltaTend   =  0.1    # end time difference between pre- and post-spike, in sec
        DeltaTsteps =  101.  # steps between start and end value
        ppp         = array([0.1,0.25,0.5,0.75,1.])

        nCases = len(ppp)

        ###########################################################
        # initialize arrays 
        deltaT = linspace(DeltaTstart,DeltaTend,DeltaTsteps)
        resultsIrr = zeros(len(ppp)*2+3)

        ###########################################################
        # simulation loop over range of deltaT values
        for i in range(len(deltaT)):
            #
            print 'deltaT : ', deltaT[i]
            
            args = column_stack((ones(nCases)*deltaT[i],ones(nCases)*frequency,ones(nCases)*frequency,ppp))
            #pdb.set_trace()
            rrr = self.pool.map(runIrregularPairSimulations,args)
            #for n in range(len(frequencies)):
            #    (synC[i,n],meanU[i,n],meanD[i,n],tD[i,n],tP[i,n]) = rrr[n]
            #pdb.set_trace()
            res1 = hstack((deltaT[i],frequency,frequency,ppp,rrr))
            resultsIrr = vstack((resultsIrr,res1))

        resultsIrr = resultsIrr[1:]

        # the array has to be flipped if Cpre>Cpost 
        if self.synChange.Cpre>self.synChange.Cpost:
            # invert time axis
            resultsIrr[:,0] = resultsIrr[:,0][::-1] 
            # reorder entire array to end up with increasing time
            resultsIrr = resultsIrr[::-1]
            
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)
            
        np.save(outputDir+'irregularSpikePairs_vs_deltaT_differentP_%s.npy' % params,resultsIrr)
        np.savetxt(outputDir+'irregularSpikePairs_vs_deltaT_differenP_%s.dat' % params,resultsIrr)


    ##################################################################################################
    # synaptic change vs rate
    ##################################################################################################
    def changeVsRateirregular_differentDeltaTs(self):
        deltaTs   = array([-0.01,-0.01,0.,0.01,0.01])   # frequency of spike-pair presentations in pairs/sec
        Freqstart = 1.    # start time difference between pre- and post-spike, in sec
        FreqTend   =  60.    # end time difference between pre- and post-spike, in sec
        FreqSteps =  60  # steps between start and end value
        ppp         = array([1.,0.4,0.,0.4,1.])

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
        if self.synChange.Cpre>self.synChange.Cpost:
            # invert time axis
            results[:,0] = results[:,0][::-1]
            # reorder entire array to end up with increasing time
            results = results[::-1]
            

        np.save(outputDir+'irregularSpikePairs_vs_rate_%s.npy' % params,results)
        np.savetxt(outputDir+'irregularSpikePairs_vs_rate_%s.dat' % params,results)


##########################################################
# main program

if __name__ == "__main__":

    if len(sys.argv[1]) > 0:
        params = sys.argv[1]
    else:
        print "No parameter set specified!"
        sys.exit(0)

    if len(sys.argv[2]) > 0:
        case = sys.argv[2]
    else:
        print 'No case specified!'
        print 'Chose one : deltaTFreq, deltaTP, freqDeltaT'
        sys.exit(0)

    calc = calculateChange(params)
    
    if case == 'deltaTFreq':
        calc.changeVsDeltaTirregular_differentFrequencies()
        calc.changeVsDeltaTregular_differentFrequencies()
        
    elif case == 'deltaTP':
        calc.changeVsDeltaTirregular_differentPs()
        
    elif case == 'freqDeltaT':
        calc.changeVsRateirregular_differentDeltaTs()

    else:
        print 'Case does not exist!'
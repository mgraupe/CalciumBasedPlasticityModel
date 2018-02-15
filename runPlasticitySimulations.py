import numpy as np
import math
import sys
import pdb

import parameter_fit_solutions


##########################################################
# Parameter sets for the calcium and synaptic weight dynamics
class runPlasticitySimulations():
    '''
        class to calculate the change in synaptic strenght
    '''

    ###############################################################################
    def __init__(self, synC,timeAboveThres,T_tot,r0):
        # read in experimental data
        self.synChange = synC
        self.tat = timeAboveThres
        self.T_total = T_tot
        self.rho0 = r0


    ##########################################################
    def runIrregularPairSimulations(self, args):
        dT = args[0]
        preRate = args[1]
        postRate = args[2]
        p = args[3]

        if self.synChange.Cpre > self.synChange.Cpost:
            (alphaD, alphaP) = self.tat.irregularSpikePairs(dT + self.synChange.D, preRate, postRate, p, deltaCa)
        else:
            (alphaD, alphaP) = self.tat.irregularSpikePairs(dT - self.synChange.D, preRate, postRate, p, deltaCa)
            self.synChange.changeInSynapticStrength(self.T_total, self.rho0, alphaD, alphaP)

        return self.synChange.mean

    ##########################################################
    def runIrregularPairSTPDeterministicSimulations(self,args):
        dT = args[0]
        preRate = args[1]
        postRate = args[2]
        p = args[3]

        (alphaD, alphaP) = self.tat.irregularSpikePairsSTPDeterministic(dT - self.synChange.D, preRate, postRate, p,
                                                                   self.synChange.tauRec, self.synChange.U)

        self.synChange.changeInSynapticStrength(self.T_total, self.rho0, alphaD, alphaP)

        return self.synChange.mean

    ##########################################################
    def runIrregularPairSTPStochasticSimulations(self,args):
        dT = args[0]
        preRate = args[1]
        postRate = args[2]
        p = args[3]

        (alphaD, alphaP) = self.tat.irregularSpikePairsSTPStochastic(dT - self.synChange.D, preRate, postRate, p,
                                                                self.synChange.tauRec,
                                                                self.synChange.U, self.synChange.Nvesicles)

        self.synChange.changeInSynapticStrength(self.T_total, self.rho0, alphaD, alphaP)

        return self.synChange.mean

    ##########################################################
    def runRegularPairSimulations(self,args):
        dT = args[0]
        preRate = args[1]
        postRate = args[2]
        p = args[3]

        (alphaD, alphaP) = self.tat.spikePairFrequency(dT - self.synChange.D, preRate)
        self.synChange.changeInSynapticStrength(self.T_total, self.rho0, alphaD, alphaP)

        return self.synChange.mean

    ##########################################################
    def runRegularPairSTPDeterministicSimulations(self,args):
        dT = args[0]
        preRate = args[1]
        postRate = args[2]
        p = args[3]

        # here alphaD and alphaP are actually the absolute times spent above threshold
        (alphaD, alphaP) = self.tat.spikePairFrequencySTPDeterministic(dT - self.synChange.D, preRate, self.synChange.Npairs,
                                                                  self.synChange.tauRec, self.synChange.U)
        # in turn T_total turns into the number of the stimulation pattern presentation
        self.synChange.changeInSynapticStrength(self.synChange.Npresentations, self.rho0, alphaD, alphaP)

        return self.synChange.mean

    ##########################################################
    def runRegularPairSTPStochasticSimulations(self,args):
        dT = args[0]
        preRate = args[1]
        postRate = args[2]
        p = args[3]

        # here alphaD and alphaP are actually the absolute times spent above threshold
        (alphaD, alphaP) = self.tat.spikePairFrequencySTPStochastic(dT - self.synChange.D, preRate, self.synChange.Npairs,
                                                               self.synChange.tauRec, self.synChange.U, self.synChange.Nvesicles)
        # in turn T_total turns into the number of the stimulation pattern presentation
        self.synChange.changeInSynapticStrength(self.synChange.Npresentations, self.rho0, alphaD, alphaP)

        return self.synChange.mean


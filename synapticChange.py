import numpy as np
import math
import sys
import pdb

import parameter_fit_solutions
from timeAboveThreshold.timeAboveThreshold import timeAboveThreshold

##########################################################
# Parameter sets for the calcium and synaptic weight dynamics
class synapticChange():
    ''' 
        class to calculate the change in synaptic strenght 
    '''
    ###############################################################################
    # synapticChange(dataCase,parameterSetName,fromFile=True,nonlinear=nl)
    def __init__(self, dataCase, parameterSetName ,fromFile=False,nonlinear=1.,USTD = None,thetaP=None,w0=0.5,par=None):
        # read in experimental data
        dataDir = '../experimental_data/'

        if thetaP is not None:
            self.thetaPfixed = True
            self.thetaP  = thetaP
        if par is not None:
            self.w0 = par.w0

        #par.thetaD,par.nonlinear,par.Npairs,par.Npresentations,par.presentationInterval,par.w0,par.tauRec

        # chose parameters from predefined set or from file
        self.choseParameterSet(parameterSetName, fromFile,tP=thetaP)
        print('parameters : ',self.tauCa, self.Cpre, self.Cpost, self.thetaD, self.thetaP, self.gammaD, self.gammaP, self.sigma, self.tau, self.rhoStar, self.D, self.beta, self.b)
        self.Npairs = 5

        if dataCase == 'sjoestroem':
            self.tauRec = par.tauRec[1]
            if USTD is None:
                self.U = par.U[1] # 0.383753
            else:
                self.U = USTD
            self.Npresentations = par.Npresentations[1] #15
            self.Nvesicles = 15
            self.presentationInterval = par.presentationInterval[1] #10
            self.Npairs = par.Npairs #5
            jesperReg = np.loadtxt(dataDir+'sjoestroem_regular_all.dat')
            jesperStoch = np.loadtxt(dataDir+'sjoestroem_stochastic.dat')
        elif dataCase == 'sjoestroemNoSTD':
            jesperReg = np.loadtxt(dataDir + 'sjoestroem_regular_all.dat')
            jesperStoch = np.loadtxt(dataDir + 'sjoestroem_stochastic.dat')
        elif dataCase == 'markram':
            self.tauRec = par.tauRec[0] #0.525
            if USTD is None:
                self.U = par.U[0] #0.46
            else:
                self.U = USTD
            self.Npresentations = par.Npresentations[0] #10
            self.Nvesicles = 15
            self.presentationInterval = par.presentationInterval[0] #4
            self.Npairs = par.Npairs # 5
            jesperReg = np.loadtxt(dataDir+'henry_regular.dat')
            jesperStoch = np.loadtxt(dataDir + 'sjoestroem_stochastic.dat')

        elif dataCase == 'tmmFacilitationJesperCaModel':
            self.U = 0.15
            self.tauFac = 1.
            self.tauDep = 0.05
            self.Npresentations = 15
            jesperReg = np.loadtxt(dataDir+'sjoestroem_regular_all.dat')
            jesperStoch = np.loadtxt(dataDir+'sjoestroem_stochastic.dat')
        elif dataCase == 'tmmfacilitationDepressionJesperCaModel':
            self.U = 0.15
            self.tauFac = 1.5
            self.tauDep = 0.3
            self.Npresentations = 15
            jesperReg = np.loadtxt(dataDir + 'sjoestroem_regular_all.dat')
            jesperStoch = np.loadtxt(dataDir + 'sjoestroem_stochastic.dat')


        if dataCase != 'genericCase':
            self.xDataReg = jesperReg[:,[0,1]]
            self.xDataReg[:,1] = self.xDataReg[:,1]/1000. # everything in sec
            self.yDataReg = jesperReg[:,2]+1. # Sjoestroem's data is normalized to 0
            self.sigmaDataReg = jesperReg[:,3]


            self.xDataStoch = jesperStoch[:,0]
            self.yDataStoch = jesperStoch[:,1]+1. # Sjoestroem's data is normalized to 0
            self.sigmaDataStoch = jesperStoch[:,2]

    ##########################################################
    # calculate UP and DOWN transition probabilities
    def transitionProbability(self, T_total, rho0, rhoBar,sigmaRhoSquared,tauEff):
        
        # argument for the Error Function
        x1 =  -(self.rhoStar - rhoBar + (rhoBar - rho0)*np.exp(-T_total/tauEff))/(np.sqrt(sigmaRhoSquared*(1.-np.exp(-2.*T_total/tauEff))))
        # transition probability
        if rho0 == 0.:
            return (1. + math.erf(x1))/2.
        else:
            return (1. - math.erf(x1))/2.


    ##########################################################
    # calculate all values for the change in synaptic strength
    def changeInSynapticStrength(self, T_total,rho0,alphaD,alphaP):
        
        # fraction of time spent above threshold
        (self.alphaD,self.alphaP) = (alphaD, alphaP)
        
        # average potentiation and depression rates
        self.GammaP = self.gammaP*self.alphaP
        self.GammaD = self.gammaD*self.alphaD
        
        # rhoBar: average value of rho in the limit of a very long protocol equivalent to the minimum of the quadratic potentia
        self.rhoBar = self.GammaP/(self.GammaP + self.GammaD)
        
        # sigmaRhoSquared L standard deviation of rho in the same limit
        self.sigmaRhoSquared = (self.alphaP + self.alphaD)*(self.sigma**2)/(self.GammaP + self.GammaD)
        # tauEff : characteristic time scale of the temporal evolution of the pdf of rho
        self.tauEff = self.tau/(self.GammaP + self.GammaD)
        #
        # UP and the DOWN transition probabilities
        self.UP   = self.transitionProbability(T_total,0.,self.rhoBar,self.sigmaRhoSquared,self.tauEff)
        self.DOWN = self.transitionProbability(T_total,1.,self.rhoBar,self.sigmaRhoSquared,self.tauEff)
        
        # mean value of the synaptic strength right at the end of the stimulation protocol
        self.meanUP   =  self.rhoBar - (self.rhoBar - 0.)*np.exp(-T_total/self.tauEff)
        self.meanDOWN =  self.rhoBar - (self.rhoBar - 1.)*np.exp(-T_total/self.tauEff)
        self.mean     =  self.rhoBar - (self.rhoBar - rho0)*np.exp(-T_total/self.tauEff)
        
        # change in synaptic strength after/before
        self.synChange = ((self.beta*(1.-self.UP) + (1.-self.beta)*self.DOWN) + (self.beta*self.UP+ (1.-self.beta)*(1.-self.DOWN))*self.b)/(self.beta + (1.-self.beta)*self.b)
        
        #synChange = synapticChange.changeSynapticStrength(synapticChange.beta,UP,DOWN,synapticChange.b)
        

    ##########################################################################################
    # calculate change for regular spike-pair vs frequency protocol
    #(xDataReg[n, 0], xDataReg[n, 1], params, par.Npairs, stoch=False)
    def calculateChangeInSynapticStrengthSTP(self, frequency, deltaT, params, Npresentations, pairCase = 'regular', DeltaTRange = None, nonlin=1):
        #####
        if len(params) == 7:
            tauCa = params[0]
            Cpre = params[1]
            Cpost = params[2]
            # thetaD = params[3]
            thetaP = self.thetaP
            gammaD = params[3]
            gammaP = params[4]
            tau = params[5]
            D = params[6]
        elif len(params) == 8:
            tauCa = params[0]
            Cpre = params[1]
            Cpost = params[2]
            # thetaD = params[3]
            thetaP = params[3]
            gammaD = params[4]
            gammaP = params[5]
            tau = params[6]
            D = params[7]

        if pairCase == 'jitter':
            DeltaTStart = DeltaTRange[0]
            DeltaTEnd = DeltaTRange[1]
        #interval = 1. / frequency
        #### (self, tauCa, Cpre, Cpost, thetaD, thetaP, nonlinear=1.,Nves=0,U=None,w0=1.):
        #print()
        #if self.U == 0:
        #    print(tauCa, Cpre, Cpost, self.thetaD, thetaP,nonlin,self.U,self.w0,deltaT - D, frequency, self.Npairs, Npresentations, self.tauRec, self.U, gammaD, gammaP,tau,self.w0,self.presentationInterval)
        #    #pdb.set_trace()
        tat = timeAboveThreshold(tauCa, Cpre, Cpost, self.thetaD, thetaP, nonlinear=nonlin,U=self.U,w0=self.w0)
        if pairCase == 'jitter':
            (self.timeD, self.timeP) = tat.spikePairFrequencySTPJitter(DeltaTStart - D, DeltaTEnd - D, frequency, self.Npairs, self.tauRec, self.U)
        elif pairCase == 'regular':
            (self.timeD, self.timeP) = tat.spikePairFrequencySTP(deltaT - D, frequency, self.Npairs, self.tauRec, self.U)
        elif pairCase == 'stochastic':
            (self.timeD, self.timeP) = tat.spikePairFrequencySTPStochastic(deltaT - D, frequency, self.Npairs, self.tauRec, self.U, self.Nvesicles)
        elif pairCase == 'fullSimulation':
            dyn = tat.spikePairFrequencySTPFullSimulation(deltaT - D, frequency, self.Npairs, Npresentations, self.tauRec, self.U, gammaD, gammaP,tau,self.w0,self.presentationInterval)
            self.timeD = 0.
            self.timeP = 1.
            #if self.U == 0:
            #    print(dyn[-1][5],dyn[-1][5]/self.w0)
            #    pdb.set_trace()
            return (dyn[-1][5]/self.w0)

        # average potentiation and depression rates
        #if not stoch:
        #    print self.timeD, self.timeP, deltaT, deltaT - D, frequency
        GammaP = gammaP * self.timeP
        GammaD = gammaD * self.timeD
        # rhoBar: average value of rho in the limit of a very long protocol equivalent to the minimum of the quadratic potentia
        try:
            rhoBar = GammaP / (GammaP + GammaD)
        except RuntimeWarning:
            print(GammaP, GammaD)
        # tauEff : characteristic time scale of the temporal evolution of the pdf of rho
        tauEff = tau / (GammaP + GammaD)
        #
        # mean value of the synaptic strength right at the end of the stimulation protocol
        mean = rhoBar - (rhoBar - self.w0) * exp(-self.Npresentations / tauEff)
        # change in synaptic strength after/before
        return (mean / self.w0)

    ##########################################################
    # chose parameter set or read file
    def choseParameterSet(self, plasticityCase,fromFile=False,tP=None):
        if plasticityCase == 'DP':
            print('DP')
            self.tauCa = 0.02 # in sec
            self.Cpre = 1.
            self.Cpost = 2.
            self.thetaD  = 1.
            self.thetaP  = 1.3
            self.gammaD  = 200.
            self.gammaP  = 321.808
            self.sigma   = 2.8284
            self.tau     = 150.
            self.rhoStar = 0.5
            self.D       = 0.0137 # in sec
            self.beta    = 0.5
            self.b       = 5.
        elif plasticityCase == 'DPD':
            print('DPD')
            self.tauCa = 0.02 # in sec
            self.Cpre = 0.9
            self.Cpost = 0.9
            self.thetaD  = 1.
            self.thetaP  = 1.3
            self.gammaD  = 250.
            self.gammaP  = 550.
            self.sigma   = 2.8284
            self.tau     = 150.
            self.rhoStar = 0.5
            self.D       = 0.0046 # in sec
            self.beta    = 0.5
            self.b       = 5.
        elif plasticityCase == 'DPDprime':
            print('DPDprime')
            self.tauCa = 0.02 # in sec
            self.Cpre = 1.
            self.Cpost = 2.
            self.thetaD  = 1
            self.thetaP  = 2.5
            self.gammaD  = 50.
            self.gammaP  = 600.
            self.sigma   = 2.8284
            self.tau     = 150.
            self.rhoStar = 0.5
            self.D       = 0.0022 # in sec
            self.beta    = 0.5
            self.b       = 5.
        elif plasticityCase == 'P':
            print('P')
            self.tauCa = 0.02 # in sec
            self.Cpre = 2.
            self.Cpost = 2.
            self.thetaD  = 1.
            self.thetaP  = 1.3
            self.gammaD  = 160.
            self.gammaP  = 257.447
            self.sigma   = 2.8284
            self.tau     = 150.
            self.rhoStar = 0.5
            self.D       = 0.0 # in sec
            self.beta    = 0.5
            self.b       = 5.
        elif plasticityCase == 'D':
            print('D')
            self.tauCa = 0.02 # in sec
            self.Cpre = 0.6
            self.Cpost = 0.6
            self.thetaD  = 1.
            self.thetaP  = 1.3
            self.gammaD  = 500.
            self.gammaP  = 550.
            self.sigma   = 5.6568
            self.tau     = 150.
            self.rhoStar = 0.5
            self.D       = 0. # in sec
            self.beta    = 0.5
            self.b       = 5.
        elif plasticityCase == 'Dprime':
            print('Dprime')
            self.tauCa = 0.02 # in sec
            self.Cpre = 1.
            self.Cpost = 2.
            self.thetaD  = 1.
            self.thetaP  = 3.5
            self.gammaD  = 60.
            self.gammaP  = 600.
            self.sigma   = 2.8284
            self.tau     = 150.
            self.rhoStar = 0.5
            self.D       = 0. # in sec
            self.beta    = 0.5
            self.b       = 5.
        elif plasticityCase == 'hippocampal slices':
            print('hippocampal slices')
            self.tauCa = 0.0488373 # in sec
            self.Cpre = 1.
            self.Cpost = 0.275865
            self.thetaD  = 1.
            self.thetaP  = 1.3
            self.gammaD  = 313.0965
            self.gammaP  = 1645.59
            self.sigma   = 9.1844
            self.tau     = 688.355
            self.rhoStar = 0.5
            self.D       = 0.0188008 # in sec
            self.beta    = 0.7
            self.b       = 5.28145
        elif plasticityCase == 'hippocampal cultures':
            print('hippocampal cultures')
            self.tauCa = 0.0119536 # in sec
            self.Cpre = 0.58156
            self.Cpost = 1.76444
            self.thetaD  = 1.
            self.thetaP  = 1.3
            self.gammaD  = 61.141
            self.gammaP  = 113.6545
            self.sigma   = 2.5654
            self.tau     = 33.7596
            self.rhoStar = 0.5
            self.D       = 0.01 # in sec
            self.beta    = 0.5
            self.b       = 36.0263
        elif plasticityCase == 'cortical slices':
            print('cortical slices')
            self.tauCa = 0.0226936 # in sec
            self.Cpre = 0.5617539
            self.Cpost = 1.23964
            self.thetaD  = 1.
            self.thetaP  = 1.3
            self.gammaD  = 331.909
            self.gammaP  = 725.085
            self.sigma   = 3.3501
            self.tau     = 346.3615
            self.rhoStar = 0.5
            self.D       = 0.0046098 # in sec
            self.beta    = 0.5
            self.b       = 5.40988
        
        elif fromFile:
            #sol = 0
            exec('sol = parameter_fit_solutions.%s' % plasticityCase,globals())
            #solll = sol[0]
            #pdb.set_trace()
            if len(sol[0]) == 7 :
                #print sol
                self.tauCa = sol[0][0]
                self.Cpre  = sol[0][1]
                self.Cpost = sol[0][2]
                self.thetaD = 1.
                self.thetaP = tP #1.2885
                self.gammaD = sol[0][3]
                self.gammaP = sol[0][4]
                self.sigma  = 1.
                self.tau    = sol[0][5]
                self.rhoStar= 0.5
                self.D      = sol[0][6]
                self.beta   = 0.5
                self.b      = 2.
            elif len(sol[0]) == 8:
                #print sol
                self.tauCa = sol[0][0]
                self.Cpre  = sol[0][1]
                self.Cpost = sol[0][2]
                self.thetaD = 1.
                self.thetaP = sol[0][3]
                self.gammaD = sol[0][4]
                self.gammaP = sol[0][5]
                self.sigma  = 1.
                self.tau    = sol[0][6]
                self.rhoStar= 0.5
                self.D      = sol[0][7]
                self.beta   = 0.5
                self.b      = 2.
            self.mse= sol[1]
        else:
            print('Choose from one of the available parameter sets!')
            sys.exit(1)



    

        

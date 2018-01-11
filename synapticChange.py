import numpy as np
import math
import sys
import pdb

import parameter_fit_solutions
            

##########################################################
# Parameter sets for the calcium and synaptic weight dynamics
class synapticChange():
    ''' 
        class to calculate the change in synaptic strenght 
    '''
    ###############################################################################
    def __init__(self, plasticityCase,fromFile=False,nonlinear=1.,dataSet='jesper'):
        
        # chose parameters from predefined set or from file
        self.choseParameterSet(plasticityCase,fromFile)
        
        # read in experimental data
        dataDir = 'experimental_data/'
        #dataDir = 'experimental_data/'

        self.Npairs = 5

        if dataSet == 'jesper':
            self.tauRec = 0.148919
            self.U = 0.383753
            self.Npresentations = 15
            jesperReg = np.loadtxt(dataDir+'sjoestroem_regular_all.dat')
            jesperStoch = np.loadtxt(dataDir+'sjoestroem_stochastic.dat')
        elif dataSet == 'henry':
            self.tauRec = 0.525
            self.U = 0.46
            self.Npresentations = 10
            jesperReg = np.loadtxt(dataDir+'henry_regular.dat')
            jesperStoch = np.loadtxt(dataDir + 'sjoestroem_stochastic.dat')

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
        
        
    ##########################################################
    # chose parameter set or read file
    def choseParameterSet(self, plasticityCase,fromFile=False):
        if plasticityCase == 'DP':
            print 'DP'
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
            print 'DPD'
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
            print 'DPDprime'
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
            print 'P'
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
            print 'D'
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
            print 'Dprime'
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
            print 'hippocampal slices'
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
            print 'hippocampal cultures'
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
            print 'cortical slices'
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
            exec('sol = parameter_fit_solutions.%s' % plasticityCase)
            
            if len(sol[0]) == 7 :
                #print sol
                self.tauCa = sol[0][0]
                self.Cpre  = sol[0][1]
                self.Cpost = sol[0][2]
                self.thetaD = 1.
                self.thetaP = 1.3
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
            print 'Choose from one of the available parameter sets!'
            sys.exit(1)



    

        

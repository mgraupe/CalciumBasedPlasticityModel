from scipy import *
from numpy import *
import pdb
import sys
import subprocess
#import commands

class timeAboveThreshold():
        ''' 
            class to calculate the fraction of time \alpha the calcium trace spends above threshold 
        '''
        ###############################################################################
        def __init__(self, tauCa, Cpre, Cpost, thetaD, thetaP, nonlinear=1.,Nves=0):
                self.tauCa = tauCa
                self.Cpre = Cpre
                self.Cpost = Cpost
                self.thetaD  = thetaD
                self.thetaP  = thetaP
                # determine eta based on nonlinearity factor and amplitudes
                self.nonlinear = nonlinear
                self.eta = (self.nonlinear*(self.Cpost + self.Cpre) - self.Cpost)/self.Cpre - 1.
                self.Nvesicles = Nves
                
        ###############################################################################
        # regular spike-pairs vs. frequency
        def spikePairFrequency(self,deltaT,frequency):
                #
                
                interval = 1./frequency
                
                timeAbove = zeros(2)
                
                # in case deltaT is larger then one interval
                if ( fabs(deltaT) > 1./frequency ):
                        deltaT = -(fabs(deltaT) - 1./frequency)
                
                # determine amplitude of the discontinous points of the calcium trace
                # post-pre
                if ( exp(1./(frequency*self.tauCa)) == NaN ) :
                        A = 0.
                        B = self.Cpost*exp(-fabs(deltaT)/self.tauCa)
                else :
                        A  = (self.Cpost + self.Cpre*exp(fabs(deltaT)/self.tauCa))/(exp(1./(frequency*self.tauCa)) - 1.)
                        B  = (self.Cpre + self.Cpost*exp((1./frequency -fabs(deltaT))/self.tauCa))/(exp(1./(frequency*self.tauCa)) - 1.)
                C  = A + self.Cpost
                D  = B + self.Cpre

                # pre-post
                if ( exp(1./(frequency*self.tauCa)) == NaN ) :
                        E = 0.
                        F = self.Cpre*exp(-fabs(deltaT)/self.tauCa)
                else :
                        E  = (self.Cpre + self.Cpost*exp(fabs(deltaT)/self.tauCa))/(exp(1./(frequency*self.tauCa)) - 1.)
                        F  = (self.Cpost + self.Cpre*exp((1./frequency -fabs(deltaT))/self.tauCa))/(exp(1./(frequency*self.tauCa)) - 1.)
                        
                G  = E + self.Cpre
                H  = F + self.Cpost
                
                # loop over depression and potentiation threshold
                for i in range(2):
                        if i == 0:
                                Ct = self.thetaD
                        elif i==1:
                                Ct = self.thetaP
                        # post-pre
                        if (deltaT < 0.):
                                if ( A <= Ct and B > Ct ) :
                                        I = self.tauCa*log(D/Ct) + fabs(deltaT)
                                elif ( A > Ct and B > Ct) :
                                        I = 1./frequency
                                elif ( A > Ct and B <= Ct ) :
                                        I = self.tauCa*log(C/Ct) + fabs(1./frequency) - fabs(deltaT)
                                elif ( A <= Ct and B <= Ct and D > Ct and C > Ct ) :
                                        I = self.tauCa*log(C/Ct) + self.tauCa*log(D/Ct)
                                elif ( A <= Ct and B <= Ct and D <= Ct and C > Ct ) :
                                        I = self.tauCa*log(C/Ct)
                                elif( A <= Ct and B <= Ct and D > Ct and C <= Ct ) :
                                        I = self.tauCa*log(D/Ct)
                                elif ( A <= Ct and B <= Ct and D <= Ct and C <= Ct ) :
                                        I = 0.
                                else :
                                        print(A, B, C, D, Ct, frequency, deltaT)
                                        print("post-pre : Problem in spikePairFrequency!")
                                        sys.exit(1)
                        # pre-post
                        else:
                                if ( E <= Ct and F > Ct ) :
                                        I = self.tauCa*log(H/Ct) + fabs(deltaT)
                                elif ( E > Ct and F > Ct) :
                                        I = 1./frequency
                                elif ( E > Ct and F <= Ct ) :
                                        I = self.tauCa*log(G/Ct) + fabs(1./frequency) - fabs(deltaT)
                                elif ( E <= Ct and F <= Ct and G > Ct and H > Ct ) :
                                        I = self.tauCa*log(G/Ct) + self.tauCa*log(H/Ct)
                                elif ( E <= Ct and F <= Ct and H <= Ct and G > Ct ) :
                                        I = self.tauCa*log(G/Ct)
                                elif ( E <= Ct and F <= Ct and H > Ct and G <= Ct ) :
                                        I = self.tauCa*log(H/Ct)
                                elif ( E <= Ct and F <= Ct and G <= Ct and H <= Ct ) :
                                        I = 0.
                                else :
                                        print(E, F, G, H, Ct, frequency, deltaT)
                                        print("pre-post : Problem in spikePairFrequency! ")
                                        sys.exit(1)
                        #
                        timeAbove[i] = I
                #
                alphaD = timeAbove[0]/interval
                alphaP = timeAbove[1]/interval
                return (alphaD,alphaP)
        ###############################################################################
        # regular spike-pairs at a given frequency, implements the nonlinear calcium model
        def spikePairFrequencyNonlinear(self,deltaT,frequency):
                #
                interval = 1./frequency
                
                timeAbove = zeros(2)
                
                # in case deltaT is larger then one interval
                if ( fabs(deltaT) > 1./frequency ):
                        deltaT = -(fabs(deltaT) - 1./frequency)
                
                # determine amplitude of the discontinous points of the calcium trace
                # post-pre
                if ( exp(1./(frequency*self.tauCa)) == NaN ) :
                        A = 0.
                        B = self.Cpost*exp(-fabs(deltaT)/self.tauCa)
                else :
                        A  = (self.Cpost + (self.eta*exp(-1./(frequency*self.tauCa))/(1.-exp(-1./(frequency*self.tauCa))) + 1.)*self.Cpre*exp(fabs(deltaT)/self.tauCa))/(exp(1./(frequency*self.tauCa)) - 1.)
                        B  = self.Cpre*(1. + self.eta/(1.-exp(-1./(frequency*self.tauCa))))/(exp(1./(frequency*self.tauCa)) - 1.) + self.Cpost*exp(-fabs(deltaT)/self.tauCa)/(1. - exp(-1./(frequency*self.tauCa)))
                C  = A + self.Cpost + self.eta*self.Cpre*exp((fabs(deltaT)-1./frequency)/self.tauCa)/(1.-exp(-1./(frequency*self.tauCa)))
                D  = B + self.Cpre

                # pre-post
                if ( exp(1./(frequency*self.tauCa)) == NaN ) :
                        E = 0.
                        F = self.Cpre*exp(-fabs(deltaT)/self.tauCa)
                else :
                        E  = (self.Cpre*(1. + self.eta/(1.-exp(-1./(frequency*self.tauCa)))) + self.Cpost*exp(fabs(deltaT)/self.tauCa))/(exp(1./(frequency*self.tauCa)) - 1.)
                        F  = (self.Cpost + self.eta/(1.-exp(-1./(frequency*self.tauCa)))*self.Cpre*exp(-fabs(deltaT)/self.tauCa))/(exp(1./(frequency*self.tauCa)) - 1.) + self.Cpre*exp(-fabs(deltaT)/self.tauCa)/(1. - exp(-1./(frequency*self.tauCa)))
                        
                G  = E + self.Cpre 
                H  = F + self.Cpost + self.eta*self.Cpre*exp(-fabs(deltaT)/self.tauCa)/(1.-exp(-1./(frequency*self.tauCa)))
                
                # loop over depression and potentiation threshold
                for i in range(2):
                        if i == 0:
                                Ct = self.thetaD
                        elif i==1:
                                Ct = self.thetaP
                        # post-pre
                        if (deltaT < 0.):
                                if ( A <= Ct and B > Ct ) :
                                        I = self.tauCa*log(D/Ct) + fabs(deltaT)
                                elif ( A > Ct and B > Ct) :
                                        I = 1./frequency
                                elif ( A > Ct and B <= Ct ) :
                                        I = self.tauCa*log(C/Ct) + fabs(1./frequency) - fabs(deltaT)
                                elif ( A <= Ct and B <= Ct and D > Ct and C > Ct ) :
                                        I = self.tauCa*log(C/Ct) + self.tauCa*log(D/Ct)
                                elif ( A <= Ct and B <= Ct and D <= Ct and C > Ct ) :
                                        I = self.tauCa*log(C/Ct)
                                elif( A <= Ct and B <= Ct and D > Ct and C <= Ct ) :
                                        I = self.tauCa*log(D/Ct)
                                elif ( A <= Ct and B <= Ct and D <= Ct and C <= Ct ) :
                                        I = 0.
                                else :
                                        print(A, B, C, D, Ct, frequency, deltaT)
                                        print ("post-pre : Problem in spikePairFrequency!")
                                        sys.exit(1)
                        # pre-post
                        else:
                                if ( E <= Ct and F > Ct ) :
                                        I = self.tauCa*log(H/Ct) + fabs(deltaT)
                                elif ( E > Ct and F > Ct) :
                                        I = 1./frequency
                                elif ( E > Ct and F <= Ct ) :
                                        I = self.tauCa*log(G/Ct) + fabs(1./frequency) - fabs(deltaT)
                                elif ( E <= Ct and F <= Ct and G > Ct and H > Ct ) :
                                        I = self.tauCa*log(G/Ct) + self.tauCa*log(H/Ct)
                                elif ( E <= Ct and F <= Ct and H <= Ct and G > Ct ) :
                                        I = self.tauCa*log(G/Ct)
                                elif ( E <= Ct and F <= Ct and H > Ct and G <= Ct ) :
                                        I = self.tauCa*log(H/Ct)
                                elif ( E <= Ct and F <= Ct and G <= Ct and H <= Ct ) :
                                        I = 0.
                                else :
                                        print(E, F, G, H, Ct, frequency, deltaT)
                                        print("pre-post : Problem in spikePairFrequency! ")
                                        sys.exit(1)
                        #
                        timeAbove[i] = I
                #
                alphaD = timeAbove[0]/interval
                alphaP = timeAbove[1]/interval
                return (alphaD,alphaP)
        
        ###############################################################################
        # one presynaptic spike and a two postsynaptic spikes, i.e., a postsynaptic burst
        def preSpikePostPair(self,deltaT,frequency,deltaBurst):
                #
                interval = 1./frequency
                
                timeAbove = zeros(2)
                #
                if ( fabs(deltaT) > 1./frequency ):
                        deltaT = -(fabs(deltaT) - 1./frequency)
                        
                #########################################
                # post-post-pre
                M = self.Cpost*exp(-(fabs(deltaBurst))/self.tauCa)
                N = self.Cpost*exp(-(fabs(deltaBurst) + fabs(deltaT))/self.tauCa) + self.Cpost*exp(-fabs(deltaT)/self.tauCa)

                if ( exp(1./(frequency*self.tauCa)) == NaN ):
                        A = 0.
                else :
                        A  = (self.Cpost + self.Cpost*exp(fabs(deltaBurst)/self.tauCa) + self.Cpre*exp((fabs(deltaBurst) + fabs(deltaT))/self.tauCa))/(exp(1./(frequency*self.tauCa)) -1.)
                        M = (self.Cpost + self.Cpre*exp(fabs(deltaT)/self.tauCa) + self.Cpost*exp((1./frequency - fabs(deltaBurst))/self.tauCa))/(exp(1./(frequency*self.tauCa)) -1.)
                        N = (self.Cpre + self.Cpost*exp((1./frequency - fabs(deltaT) - fabs(deltaBurst))/self.tauCa) + self.Cpost*exp((1./frequency - fabs(deltaT))/self.tauCa))/(exp(1./(frequency*self.tauCa)) -1.)

                O  = M + self.Cpost
                P  = N + self.Cpre
                
                #########################################
                # post-pre-post
                Q = self.Cpost*exp(-(fabs(deltaBurst) - fabs(deltaT))/self.tauCa)
                R = self.Cpost*exp(- fabs(deltaBurst)/self.tauCa) + self.Cpre*exp(-fabs(deltaT)/self.tauCa)
                if ( exp(1./(frequency*self.tauCa)) == NaN ) :
                        A = 0.
                else :
                        A  = (self.Cpost + self.Cpre*exp((fabs(deltaBurst)-fabs(deltaT))/self.tauCa) + self.Cpost*exp(fabs(deltaBurst)/self.tauCa))/(exp(1./(frequency*self.tauCa)) -1.)
                        Q = (self.Cpre + self.Cpost*exp(fabs(deltaT)/self.tauCa) + self.Cpost*exp((1./frequency - (fabs(deltaBurst) - fabs(deltaT)))/self.tauCa))/(exp(1./(frequency*self.tauCa)) -1.)
                        R = (self.Cpost + self.Cpost*exp((1./frequency - fabs(deltaBurst))/self.tauCa) + self.Cpre*exp((1./frequency - fabs(deltaT))/self.tauCa))/(exp(1./(frequency*self.tauCa)) -1.)

                S  = Q + self.Cpre
                T  = R + self.Cpost
                
                
                #########################################
                # pre-post-post
                U = self.Cpre*exp(-(fabs(deltaT)-fabs(deltaBurst))/self.tauCa)
                V = self.Cpre*exp(- fabs(deltaT)/self.tauCa) + self.Cpost*exp(-fabs(deltaBurst)/self.tauCa)
                
                if ( exp(1./(frequency*self.tauCa)) == NaN ):
                        A = 0.
                else :
                        A  = (self.Cpre + self.Cpost*exp((fabs(deltaT)-fabs(deltaBurst))/self.tauCa) + self.Cpost*exp(fabs(deltaT)/self.tauCa))/(exp(1./(frequency*self.tauCa)) -1.)
                        U = (self.Cpost + self.Cpost*exp(fabs(deltaBurst)/self.tauCa) + self.Cpre*exp((1./frequency - (fabs(deltaT) - fabs(deltaBurst)))/self.tauCa))/(exp(1./(frequency*self.tauCa)) -1.)
                        V = (self.Cpost + self.Cpre*exp((1./frequency - fabs(deltaT))/self.tauCa) + self.Cpost*exp((1./frequency - fabs(deltaBurst))/self.tauCa))/(exp(1./(frequency*self.tauCa)) -1.)

                W  = U + self.Cpost
                X  = V + self.Cpost
	

                for i in range(2):
                        if i == 0:
                                Ct = self.thetaD
                        elif i==1:
                                Ct = self.thetaP
                        # 2 x post and 1 x pre	
                        # post-post-pre
                        if ( deltaT < 0.) :
                                # print "post-post-pre", deltaT, deltaBurst
                                if ( M > Ct and N > Ct and A > Ct) :
                                        Int = 1./frequency
                                elif ( M > Ct and N > Ct and A <= Ct) :
                                        Int = self.tauCa*log(P/Ct) + fabs(deltaT) + fabs(deltaBurst)
                                elif ( M > Ct and N <= Ct and P > Ct) :
                                        Int = self.tauCa*log(O/Ct) + fabs(deltaBurst) + self.tauCa*log(P/Ct)
                                elif ( M > Ct and P <= Ct) :
                                        Int = self.tauCa*log(O/Ct) + fabs(deltaBurst)
                                elif ( (A+self.Cpost) > Ct and M <= Ct and N > Ct ) :
                                        Int = self.tauCa*log((A+self.Cpost)/Ct) + self.tauCa*log(P/Ct) + fabs(deltaT)
                                elif ( (A+self.Cpost) > Ct and M <= Ct and N <= Ct and P > Ct ) :
                                        Int = self.tauCa*log((A+self.Cpost)/Ct) + self.tauCa*log(O/Ct) +  self.tauCa*log(P/Ct)
                                elif ( (A+self.Cpost) > Ct and M <= Ct and P <= Ct ) :
                                        Int = self.tauCa*log((A+self.Cpost)/Ct) + self.tauCa*log(O/Ct)
                                elif ( (A+self.Cpost) > Ct and O <= Ct and P <= Ct) :
                                        Int = self.tauCa*log((A+self.Cpost)/Ct)
                                elif ( (A+self.Cpost) > Ct and O <= Ct and P > Ct) :
                                        Int = self.tauCa*log((A+self.Cpost)/Ct) + self.tauCa*log(P/Ct)
                                elif ( (A+self.Cpost) <= Ct and N > Ct ) :
                                        Int = self.tauCa*log(P/Ct) + fabs(deltaT)
                                elif ( (A+self.Cpost) <= Ct and O > Ct and N <= Ct and P > Ct) :
                                        Int = self.tauCa*log(O/Ct) + self.tauCa*log(P/Ct)
                                elif ( (A+self.Cpost) <= Ct and O > Ct and P <= Ct) :
                                        Int = self.tauCa*log(O/Ct)
                                elif ( (A+self.Cpost) <= Ct and O <= Ct and P > Ct) :
                                        Int = self.tauCa*log(P/Ct)
                                elif ( (A+self.Cpost) <= Ct and O <= Ct and P <= Ct) :
                                        Int = 0.
                                else :
                                        print("post-post-pre : Problem !")
                                        print(deltaT, A, (A+self.Cpost), M, N, O, P, Int)
                                        sys.exit(1)
                        # post-pre-post 
                        elif ( deltaT >= 0.  and deltaT <= deltaBurst )  :
                                # print "post-pre-post", deltaT, deltaBurst 
                                if ( Q > Ct and R > Ct and A > Ct) :
                                        Int = 1./frequency
                                elif ( Q > Ct and R > Ct and A <= Ct) :
                                        Int = self.tauCa*log(T/Ct) + fabs(deltaBurst)
                                elif ( Q > Ct and R <= Ct and T > Ct ) :
                                        Int = self.tauCa*log(S/Ct) + fabs(deltaBurst) - fabs(deltaT) +  self.tauCa*log(T/Ct)	
                                elif ( Q > Ct and T <= Ct and A <= Ct ) :
                                        Int = self.tauCa*log(S/Ct) + fabs(deltaBurst) - fabs(deltaT)	
                                elif ( (A+self.Cpost) > Ct and Q <= Ct and R > Ct ) :
                                        Int = self.tauCa*log((A+self.Cpost)/Ct) + self.tauCa*log(T/Ct) + fabs(deltaT)
                                elif ( (A+self.Cpost) > Ct and Q <= Ct and S > Ct and R <= Ct and T > Ct ) :
                                        Int = self.tauCa*log((A+self.Cpost)/Ct) + self.tauCa*log(S/Ct) + self.tauCa*log(T/Ct)
                                elif ( (A+self.Cpost) > Ct and Q <= Ct and S > Ct and T <= Ct ) :
                                        Int = self.tauCa*log((A+self.Cpost)/Ct) + self.tauCa*log(S/Ct)
                                elif ( (A+self.Cpost) > Ct and S <= Ct and T > Ct ) :
                                        Int = self.tauCa*log((A+self.Cpost)/Ct) + self.tauCa*log(T/Ct)
                                elif ( (A+self.Cpost) > Ct and S <= Ct and T <= Ct ) :
                                        Int = self.tauCa*log((A+self.Cpost)/Ct)
                                elif ( (A+self.Cpost) <= Ct and S > Ct and T <= Ct ) :
                                        Int = self.tauCa*log(S/Ct)
                                elif ( (A+self.Cpost) <= Ct and S <= Ct and T > Ct) :
                                        Int = self.tauCa*log(T/Ct)
                                elif ( (A+self.Cpost) <= Ct and S > Ct and R <= Ct and T > Ct) :
                                        Int = self.tauCa*log(S/Ct) + self.tauCa*log(T/Ct)
                                elif ( (A+self.Cpost) <= Ct  and R > Ct) :
                                        Int = self.tauCa*log(T/Ct) + fabs(deltaT)
                                elif ( (A+self.Cpost) <= Ct and S <= Ct and T <= Ct) :
                                        Int = 0.
                                else :
                                        print("post-pre-post : Problem !")
                                        print(deltaT, A, Q, R, (A+self.Cpost), S, T, Int)
                                        sys.exit(1)
                        # pre-post-post 
                        elif ( deltaT >= 0.  and deltaT > deltaBurst)  :
                                # print "pre-post-post", deltaT, deltaBurst
                                if ( U > Ct and V > Ct and A > Ct ) :
                                        Int = 1./frequency
                                elif ( U > Ct and V > Ct and A <= Ct) :
                                        Int = self.tauCa*log(X/Ct) + fabs(deltaT)
                                elif ( U > Ct and V <= Ct and X > Ct ) :
                                        Int = self.tauCa*log(W/Ct) + fabs(deltaT) - fabs(deltaBurst) + self.tauCa*log(X/Ct)
                                elif ( U > Ct and X <= Ct ) :
                                        Int = self.tauCa*log(W/Ct) + fabs(deltaT) - fabs(deltaBurst)
                                elif ( (A+self.Cpre) > Ct and U <= Ct and  V > Ct) :
                                        Int = self.tauCa*log((A+self.Cpre)/Ct) + self.tauCa*log(X/Ct) + fabs(deltaBurst)	
                                elif ( (A+self.Cpre) > Ct and U <= Ct and  V <= Ct and W > Ct and X > Ct ) :
                                        Int = self.tauCa*log((A+self.Cpre)/Ct) + self.tauCa*log(W/Ct) + self.tauCa*log(X/Ct)
                                elif ( (A+self.Cpre) > Ct and U <= Ct and W > Ct and X <= Ct ) :
                                        Int = self.tauCa*log((A+self.Cpre)/Ct) + self.tauCa*log(W/Ct)
                                elif ( (A+self.Cpre) > Ct and W <= Ct and X > Ct ) :
                                        Int = self.tauCa*log((A+self.Cpre)/Ct) + self.tauCa*log(X/Ct)
                                elif ( (A+self.Cpre) > Ct and W <= Ct and X <= Ct ) :
                                        Int = self.tauCa*log((A+self.Cpre)/Ct)
                                elif ( (A+self.Cpre) <= Ct and V > Ct ) :
                                        Int = self.tauCa*log(X/Ct) + fabs(deltaBurst)
                                elif ( (A+self.Cpre) <= Ct  and V <= Ct and W > Ct and X > Ct ) :
                                        Int = self.tauCa*log(W/Ct) + self.tauCa*log(X/Ct)
                                elif ( (A+self.Cpre) <= Ct  and W <= Ct and X > Ct ) :
                                        Int = self.tauCa*log(X/Ct)
                                elif ( (A+self.Cpre) <= Ct  and W > Ct and X <= Ct ) :
                                        Int = self.tauCa*log(W/Ct)
                                elif ( (A+self.Cpre) <= Ct  and W <= Ct and X <= Ct ) :
                                        Int = 0.
                                else :
                                        print("pre-post-post : Problem !")
                                        print(deltaT, A, (A+self.Cpre), U, V, W, X, Int)
                                        sys.exit(1)
                        else :
                                print("Problem in preSpikePostPair routine!")
                                sys.exit(1)
                        #
                        timeAbove[i] = Int
                #
                alphaD = timeAbove[0]/interval
                alphaP = timeAbove[1]/interval
                return (alphaD,alphaP)
        ###############################################################################
        # irregular spike-pairs, the numerical integration is run in an external C++ code for performance improvment
        def irregularSpikePairs(self,deltaT,preRate,postRate,ppp,deltaCa):
                
                # linear calcium dynamics, numerical integration possible
                if self.nonlinear == 1.:
                    #print 'time above threshold : integration for LINEAR calcium dynamics'
                    # the first argument calcium ampliutde has to be smaller than the second
                    if self.Cpre>self.Cpost:
                            arguments = str(deltaT) + ' ' + str(self.tauCa) + ' ' + str(self.Cpost) + ' ' + str(self.Cpre) + ' ' + str(self.thetaD) + ' ' + str(self.thetaP) + ' ' + str(preRate) + ' ' + str(postRate) + ' ' + str(ppp) + ' ' + str(deltaCa)
                    else:
                            arguments = str(deltaT) + ' ' + str(self.tauCa) + ' ' + str(self.Cpre) + ' ' + str(self.Cpost) + ' ' + str(self.thetaD) + ' ' + str(self.thetaP) + ' ' + str(preRate) + ' ' + str(postRate) + ' ' + str(ppp) + ' ' + str(deltaCa)
        
                    #print arguments
                    # depcreciated in python 3
                    # (out,err) = commands.getstatusoutput('./timeAboveThreshold/poissonPairs_timeAboveThreshold ' + arguments)
                    tt = subprocess.check_output('./timeAboveThreshold/poissonPairs_timeAboveThreshold ' + arguments)
                    pdb.set_trace()
                    alphaD = float(err.split()[0])
                    alphaP = float(err.split()[1])
                    
                # nonlinear calcium dynamics
                else:
                    #print 'time above threshold : NONLINEAR calcium dynamics'
                    
                    # construction of the spike train
                    #tStart = 0.1 # start time at 100 ms
                    
                    random.seed(7)
                    
                    tPre = []
                    tPostCorr = []
                    tPostInd = []
                    
                    tPre.append(0)
                    tPostInd.append(0)
                    
                    for i in range(10000):
                        tPre.append(tPre[-1] + random.exponential(1./preRate))
                        if rand()<ppp:
                            tPostCorr.append(tPre[-1]+deltaT)
                        if (postRate-ppp*preRate) > 0.:
                            tPostInd.append(tPostInd[-1] + random.exponential(1./(postRate-ppp*preRate)))
                    
                    
                    tPost = tPostCorr + tPostInd[1:]
                    
                    tPostSorted = sorted(tPost, key=lambda tPost: tPost)
                    
                    tAll = zeros((len(tPre[1:]) + len(tPostSorted),3))
                    
                    tAll[:,0] = hstack((tPre[1:],tPostSorted))
                    tAll[:,1] = hstack((zeros(len(tPre[1:])),ones(len(tPostSorted))))
                    tAll[:,2] = hstack((repeat(self.Cpre,Npres),repeat(self.Cpost,Npres)))

                    tList = tAll.tolist()
                    tListSorted = sorted(tList, key=lambda tList: tList[0])
                    
                    #tListSorted.append([Npres/freq,2])

                    (tD, tP) = self.eventBasedIntegration(tListSorted)

                    alphaD = tD/tListSorted[-1][0]
                    alphaP = tP/tListSorted[-1][0]
                    #print alphaD, alphaP
                ####################################################################
                
                return (alphaD,alphaP)

        ###############################################################################
        # irregular spike-pairs and deterministic short-term plasticity, event-based integration
        def irregularSpikePairsSTPDeterministic(self, deltaT, preRate, postRate, ppp, tauRec, U):


                # construction of the spike train
                # tStart = 0.1 # start time at 100 ms

                random.seed(7)

                tPre = []
                tPostCorr = []
                tPostInd = []

                tPre.append(0)
                tPostInd.append(0)

                for i in range(50000):
                        tPre.append(tPre[-1] + random.exponential(1. / preRate))
                        if rand() < ppp:
                                tPostCorr.append(tPre[-1] + deltaT)
                        if (postRate - ppp * preRate) > 0.:
                                tPostInd.append(
                                        tPostInd[-1] + random.exponential(1. / (postRate - ppp * preRate)))

                tPost = tPostCorr + tPostInd[1:]

                tPostSorted = sorted(tPost, key=lambda tPost: tPost)

                cpre = zeros(len(tPre[1:]))
                # deterministic STD model implementation
                if U != 0.:
                        cpre[0] = 1.
                        for i in range(1, len(tPre[1:])):
                                cpre[i] = 1. - (1. - (cpre[i - 1] - U * cpre[i - 1])) * exp(-(tPre[i+1]-tPre[i])/tauRec)
                        cpre *= U * self.Cpre
                else:
                        cpre[:] = self.Cpre


                tAll = zeros((len(tPre[1:]) + len(tPostSorted), 3))

                tAll[:, 0] = hstack((tPre[1:], tPostSorted))
                tAll[:, 1] = hstack((zeros(len(tPre[1:])), ones(len(tPostSorted))))
                tAll[:, 2] = hstack((cpre,repeat(self.Cpost,len(tPostSorted))))
                tList = tAll.tolist()
                tListSorted = sorted(tList, key=lambda tList: tList[0])

                # tListSorted.append([Npres/freq,2])

                ###########################################################
                # event-based integration
                (tD, tP) = self.eventBasedIntegration(tListSorted)

                alphaD = tD / tListSorted[-1][0]
                alphaP = tP / tListSorted[-1][0]
                # print alphaD, alphaP

                return (alphaD, alphaP)
        ###############################################################################
        # irregular spike-pairs and deterministic short-term plasticity, event-based integration of calcium and synaptic strength
        # tat.irregularSpikePairsSTPDeterminisitcFullSim(dT-synChange.D,preRate,postRate,p,synChange.tauRec,synChange.U,T_total,rho0,synChange.gammaD,synChange.gammaP)
        def irregularSpikePairsSTPDeterminisitcFullSim(self, deltaT, preRate, postRate, ppp, tauRec, U, T_total, rho0, tau, gammaD, gammaP):

                Nrepetitions = 10000 # 10000

                # construction of the spike train
                # tStart = 0.1 # start time at 100 ms
                #random.seed(7)
                wFinal = []
                for n in range(Nrepetitions):
                        #print(n)
                        tPre = []
                        tPostCorr = []
                        tPostInd = []

                        tPre.append(0)
                        tPostInd.append(0)

                        for i in range(100000):
                                tPre.append(tPre[-1] + random.exponential(1. / preRate))
                                if (postRate - ppp * preRate) > 0.:
                                        tPostInd.append(tPostInd[-1] + random.exponential(1. / (postRate - ppp * preRate)))
                                if tPre[-1]>=T_total: # abort of end of simulation is reached
                                        tPre = tPre[:-1]
                                        break
                                if rand() < ppp:
                                        tPostCorr.append(tPre[-1] + deltaT)

                        tPost = tPostCorr + tPostInd[1:]
                        tPre = tPre[1:]

                        tPostSorted = sorted(tPost, key=lambda tPost: tPost)
                        tPostSorted = [item for item in tPostSorted if (item>0. and item<T_total)] # remove negative time points, which can occur due to negative delta t values

                        tAll = zeros((len(tPre) + len(tPostSorted), 2))

                        tAll[:, 0] = hstack((tPre, tPostSorted))
                        tAll[:, 1] = hstack((zeros(len(tPre)), ones(len(tPostSorted))))
                        tList = tAll.tolist()
                        tListSorted = sorted(tList, key=lambda tList: tList[0])
                        tListSorted.append([T_total, 2])
                        dynamics = self.fullEventBasedSimulation(tListSorted,rho0,tau,gammaD,gammaP,tauRec, U)
                        wFinal.append(dynamics[-1][5])
                        #pdb.set_trace()
                #
                #dyn = asarray(dyn)
                #plt.plot(dyn[:,0])
                #plt.show()
                #pdb.set_trace()

                return mean(wFinal)
        ###############################################################################
        #synCh = tat.regularSpikePairsSTPDeterminisitcFullSim(dT - synChange.D, preRate, postRate, p, synChange.tauRec, synChange.U, T_total, rho0, synChange.tau, synChange.gammaD,synChange.gammaP)
        def regularSpikePairsSTPDeterminisitcFullSim(self,deltaT, preRate, postRate, ppp, tauRec, U, T_total, rho0, tau, gammaD, gammaP):
                # set timing of spikes
                tStart = 0.15  # start time at 100 ms
                # Npres = Npres * 12

                Npairs = int(T_total*preRate)

                #tStartPacket = tStart + arange(Npres) * presentationInterval
                #tStartPacket = repeat(tStartPacket, Npairs)
                tPre = tStart + arange(Npairs) / preRate
                tPost = tPre + deltaT

                tAll = zeros((2 * Npairs, 2))

                # print cpre
                tAll[:, 0] = hstack((tPre, tPost))
                tAll[:, 1] = hstack((zeros(Npairs), ones(Npairs )))
                # tAll[:, 2] = hstack((cpre,repeat(self.Cpost,Npres)))
                tList = tAll.tolist()
                # sort list to pre- and post spike according to increasing spike time
                tListSorted = sorted(tList, key=lambda tList: tList[0])
                # ad an additional time point (t > t_pre,t_post) to evaluate the dynamics after the last spike
                tListSorted.append([T_total, 2])

                dynamics = self.fullEventBasedSimulation(tListSorted,rho0,tau,gammaD,gammaP,tauRec, U)
                #pdb.set_trace()
                return dynamics[-1][5]

        ###############################################################################
        # run the full event based simulation tracking calcium and synaptic weigth
        def fullEventBasedSimulation(self,tListSorted,rho0,tau,gammaD,gammaP,tauRec, U):
                # determine w dynamics parameters based on gammaP, gammaD, and tau
                # for dynamics above thetaD and below thetaP
                w_bar_D = 0.
                tau_prime_D = tau / gammaD
                # for dynamics above thetaP
                w_bar_PD = gammaP / (gammaP + gammaD)
                tau_prime_PD = tau / (gammaP + gammaD)

                # run simulation
                dyn = []
                # CaTotal, CaPre, x ... preResources, CaPost, time, w
                dyn.append([0., 0., 1., 0., 0., rho0])
                pre = 0
                post = 0
                tpreOld = None
                for i in tListSorted:
                        #
                        caTotOld = dyn[-1][0]
                        caPreOld = dyn[-1][1]
                        caPreSTDOld = dyn[-1][2]
                        caPostOld = dyn[-1][3]
                        tOld = dyn[-1][4]
                        w = dyn[-1][5]
                        # caTotTemp  = caTotOld*exp(-(i[0]-tOld)/self.tauCa)
                        caPreTemp = caPreOld * exp(-(i[0] - tOld) / self.tauCa)
                        caPostTemp = caPostOld * exp(-(i[0] - tOld) / self.tauCa)
                        caTotTemp = caPreTemp + caPostTemp
                        # time above potentiation threshold
                        if caTotOld > self.thetaP:
                                tendP = i[0] if (caTotTemp > self.thetaP) else ((self.tauCa) * log(caTotOld / self.thetaP) + tOld)
                                high = True
                        else:
                                high = False
                        # time above depression threshold
                        if (caTotOld > self.thetaD) and (caTotTemp < self.thetaP):
                                tstartD = tOld if (caTotOld < self.thetaP) else ((self.tauCa) * log(caTotOld / self.thetaP) + tOld)
                                tendD = i[0] if (caTotTemp > self.thetaD) else ((self.tauCa) * log(caTotOld / self.thetaD) + tOld)
                                low = True
                        else:
                                low = False
                        # time below both thresholds
                        if (caTotOld < self.thetaD) or (caTotTemp < self.thetaD):
                                tstartDet = tOld if (caTotOld < self.thetaD) else ((self.tauCa) * log(caTotOld / self.thetaD) + tOld)
                                deterministic = True
                        else:
                                deterministic = False
                        # performing update of w
                        if high:
                                w = w_bar_PD + (w - w_bar_PD) * exp(-(tendP - tOld) / tau_prime_PD)
                        if low:
                                w = w_bar_D + (w - w_bar_D) * exp(-(tendD - tstartD) / tau_prime_D)
                        if deterministic:
                                # possible update of w below thresholds in case of double-well or piecewise quadratic potential
                                # w = ...
                                # no update required for flat potential - line-attractor
                                pass
                        # postsynaptic spike
                        if i[1] == 1:
                                caPostTemp += self.Cpost + self.eta * caPreTemp
                                post += 1
                                caPreSTDTemp = caPreSTDOld
                        # presynaptic spike
                        # note that the influence of the presynaptically evoked calcium transient depends on w
                        if i[1] == 0:
                                if U != 0:
                                        caPreSTDTemp = 1. if (tpreOld is None) else (1. - (1. - (caPreSTDOld - U * caPreSTDOld)) * exp(-(i[0] - tpreOld) / tauRec))
                                        caPreTemp += w * U * self.Cpre * caPreSTDTemp
                                else:
                                        caPreTemp += w * self.Cpre
                                pre += 1
                                tpreOld = i[0]
                        try :
                                caPreSTDTemp
                        except NameError:
                                caPreSTDTemp = caPreSTDOld
                        else:
                                pass

                        caTotTemp = caPreTemp + caPostTemp
                        dyn.append([caTotTemp, caPreTemp, caPreSTDTemp, caPostTemp, i[0], w])

                return dyn
        ###############################################################################
        # irregular spike-pairs and deterministic short-term plasticity, event-based integration
        def irregularSpikePairsSTPStochastic(self, deltaT, preRate, postRate, ppp, tauRec, pRelease, Nves):

                NpreSpikes = 5000

                NrepetitionsStoch = 100

                q = self.Cpre/Nves

                # construction of the spike train
                # tStart = 0.1 # start time at 100 ms

                random.seed(7)

                tPre = []
                tPostCorr = []
                tPostInd = []

                tPre.append(0)
                tPostInd.append(0)

                for i in range(NpreSpikes):
                        tPre.append(tPre[-1] + random.exponential(1. / preRate))
                        if rand() < ppp:
                                tPostCorr.append(tPre[-1] + deltaT)
                        if (postRate - ppp * preRate) > 0.:
                                tPostInd.append(
                                        tPostInd[-1] + random.exponential(1. / (postRate - ppp * preRate)))

                tPost = tPostCorr + tPostInd[1:]

                tPostSorted = sorted(tPost, key=lambda tPost: tPost)
                tPre = asarray(tPre[1:])

                #######################################
                # stochastic STD model implementation
                #ampStoch = np.zeros((NrepetitionsStoch, len(tPre)))
                timesAbove = zeros((NrepetitionsStoch, 2))
                #tP = 0.
                for r in range(NrepetitionsStoch):
                        Vesicles = ones((len(tPre), Nves))
                        Release = random.rand(len(tPre), Nves) < pRelease
                        #VesTimes = transpose(np.tile(tPre, (Nves, 1)))
                        for i in range(len(tPre)):
                                nRel = sum(Release[i])
                                releaseSites = argwhere(Release[i] == True)
                                emptyTimes = random.exponential(tauRec, nRel)
                                for n in range(nRel):
                                        if not Vesicles[i, releaseSites[n]] == 0.:
                                                recoveryMask = ((tPre - tPre[i]) < emptyTimes[n]) & ((tPre - tPre[i]) > 0)
                                                Vesicles[recoveryMask, releaseSites[n]] = 0
                        #print 'after amp det'
                        cpreStoch = q*(sum(Vesicles*Release,axis=1))

                        ######################################
                        tAll = zeros((len(tPre) + len(tPostSorted), 3))

                        tAll[:, 0] = hstack((tPre, tPostSorted))
                        tAll[:, 1] = hstack((zeros(len(tPre)), ones(len(tPostSorted))))
                        tAll[:, 2] = hstack((cpreStoch,repeat(self.Cpost,len(tPostSorted))))
                        tList = tAll.tolist()
                        tListSorted = sorted(tList, key=lambda tList: tList[0])

                        # tListSorted.append([Npres/freq,2])

                        ###########################################################
                        # event-based integration
                        (timesAbove[r,0], timesAbove[r,1]) = self.eventBasedIntegration(tListSorted)
                        #print r

                alphaD = average(timesAbove[:,0]) / tListSorted[-1][0]
                alphaP = average(timesAbove[:,1]) / tListSorted[-1][0]
                # print alphaD, alphaP

                return (alphaD, alphaP)


        ###############################################################################
        # irregular spike-pairs and short-term plasticity, event-based integration the numerical integration is run in an external C++ code for performance improvment
        def irregularSpikePairsTMM(self, deltaT, preRate, postRate, ppp, U, tauFac, tauDep):


                # construction of the spike train
                # tStart = 0.1 # start time at 100 ms

                random.seed(7)

                tPre = []
                tPostCorr = []
                tPostInd = []

                tPre.append(0)
                tPostInd.append(0)

                for i in range(10000):
                        tPre.append(tPre[-1] + random.exponential(1. / preRate))
                        if rand() < ppp:
                                tPostCorr.append(tPre[-1] + deltaT)
                        if (postRate - ppp * preRate) > 0.:
                                tPostInd.append(
                                        tPostInd[-1] + random.exponential(1. / (postRate - ppp * preRate)))

                tPost = tPostCorr + tPostInd[1:]

                tPostSorted = sorted(tPost, key=lambda tPost: tPost)

                cpre = zeros(len(tPre[1:]))
                x    = zeros(len(tPre[1:]))
                u    = zeros(len(tPre[1:]))
                if U != 0.:
                        u[0] = U
                        x[0] = 1.
                        cpre[0] = u[0] * x[0]
                        for i in range(1, len(tPre[1:])):
                                u[i] = u[i - 1] * exp(-(tPre[i+1]-tPre[i]) / tauFac) + U * (1. - u[i - 1] * exp(-(tPre[i+1]-tPre[i]) / tauFac))
                                x[i] = 1. - (1. - (x[i - 1] - u[i - 1] * x[i - 1])) * exp(-(tPre[i+1]-tPre[i]) / tauDep)
                                cpre[i] = u[i] * x[i]
                        cpre *= self.Cpre
                        #cpre[0] = 1.
                        #for i in range(1, len(tPre[1:])):
                        #        cpre[i] = 1. - (1. - (cpre[i - 1] - U * cpre[i - 1])) * exp(-(tPre[i+1]-tPre[i])/tauRec)
                        #cpre *= U * self.Cpre
                else:
                        cpre[:] = self.Cpre


                tAll = zeros((len(tPre[1:]) + len(tPostSorted), 3))

                tAll[:, 0] = hstack((tPre[1:], tPostSorted))
                tAll[:, 1] = hstack((zeros(len(tPre[1:])), ones(len(tPostSorted))))
                tAll[:, 2] = hstack((cpre,repeat(self.Cpost,len(tPostSorted))))
                tList = tAll.tolist()
                tListSorted = sorted(tList, key=lambda tList: tList[0])

                # tListSorted.append([Npres/freq,2])

                ###########################################################
                # event-based integration
                (tD, tP) = self.eventBasedIntegration(tListSorted)

                alphaD = tD / tListSorted[-1][0]
                alphaP = tP / tListSorted[-1][0]
                # print alphaD, alphaP

                return (alphaD, alphaP)

        ###############################################################################
        # stochastic Sjoestroem 2001 protocol
        def spikePairStochasticFrequency(self,DeltaTStart,DeltaTEnd,freq,Npres):
                tStart = 0.1 # start time at 100 ms
                
                Npres = Npres*12
                timeAbove = zeros((1,2))
                
                tD = 0.
                tP = 0.
                random.seed(7)
                tPre = arange(Npres)/freq + tStart + (DeltaTStart + rand(Npres)*(DeltaTEnd-DeltaTStart))
                tPost = tPre +  (DeltaTStart + rand(Npres)*(DeltaTEnd-DeltaTStart))
                
                tAll = zeros((2*Npres,3))
                
                tAll[:,0] = hstack((tPre,tPost))
                tAll[:,1] = hstack((zeros(Npres),ones(Npres)))
                tAll[:,2] = hstack((repeat(self.Cpre, Npres), repeat(self.Cpost, Npres)))
                tList = tAll.tolist()
                
                tListSorted = sorted(tList, key=lambda tList: tList[0])
                
                tListSorted.append([Npres/freq,2])

                ###########################################################
                # event-based integration
                (tD, tP) = self.eventBasedIntegration(tListSorted)

                alphaD = tD/(float(Npres))
                alphaP = tP/(float(Npres))
                
                return (tD/float(Npres),tP/float(Npres))
            
                #return (alphaD,alphaP)
                
        ###############################################################################
        def spikePairFrequencySTPDeterministic(self, deltaT, freq, Npres, tauRec, U ):
                tStart = 0.1  # start time at 100 ms

                #Npres = Npres * 12
                #timeAbove = zeros((1, 2))

                #random.seed(7)
                tPre  = arange(Npres) / freq + tStart
                tPost = tPre + deltaT

                tAll = zeros((2 * Npres, 3))

                cpre = zeros(Npres)
                if U != 0:
                    cpre[0] = 1.
                    for i in range(1, Npres):
                        cpre[i] = 1. - (1. - (cpre[i-1] - U*cpre[i-1]))*exp(-1./(freq * tauRec))
                    cpre *= U*self.Cpre
                else:
                    cpre[:] = self.Cpre

                #print cpre
                tAll[:, 0] = hstack((tPre, tPost))
                tAll[:, 1] = hstack((zeros(Npres), ones(Npres)))
                tAll[:, 2] = hstack((cpre,repeat(self.Cpost,Npres)))
                tList = tAll.tolist()

                tListSorted = sorted(tList, key=lambda tList: tList[0])

                tListSorted.append([Npres / freq + tStart, 2,0])

                ###########################################################
                # event-based integration
                (tD, tP) = self.eventBasedIntegration(tListSorted)

                return (tD , tP )

        ###############################################################################
        def spikePairFrequencySTPStochastic(self, deltaT, freq, Npres, tauRec, pRelease, Nves):

                tStart = 0.1  # start time at 100 ms
                NrepetitionsStoch = 100
                q = self.Cpre / Nves

                tPre  = arange(Npres) / freq + tStart
                tPost = tPre + deltaT

                #######################################
                # stochastic STD model implementation
                # ampStoch = np.zeros((NrepetitionsStoch, len(tPre)))
                timesAbove = zeros((NrepetitionsStoch, 2))
                # tP = 0.
                for r in range(NrepetitionsStoch):
                        Vesicles = ones((len(tPre), Nves))
                        Release = random.rand(len(tPre), Nves) < pRelease
                        # VesTimes = transpose(np.tile(tPre, (Nves, 1)))
                        for i in range(len(tPre)):
                                nRel = sum(Release[i])
                                releaseSites = argwhere(Release[i] == True)
                                emptyTimes = random.exponential(tauRec, nRel)
                                for n in range(nRel):
                                        if not Vesicles[i, releaseSites[n]] == 0.:
                                                recoveryMask = ((tPre - tPre[i]) < emptyTimes[n]) & (
                                                                (tPre - tPre[i]) > 0)
                                                Vesicles[recoveryMask, releaseSites[n]] = 0
                        # print 'after amp det'
                        cpreStoch = q * (sum(Vesicles * Release, axis=1))

                        ######################################
                        tAll = zeros((len(tPre) + len(tPost), 3))

                        tAll[:, 0] = hstack((tPre, tPost))
                        tAll[:, 1] = hstack((zeros(len(tPre)), ones(len(tPost))))
                        tAll[:, 2] = hstack((cpreStoch, repeat(self.Cpost, len(tPost))))
                        tList = tAll.tolist()
                        tListSorted = sorted(tList, key=lambda tList: tList[0])

                        tListSorted.append([Npres / freq + tStart, 2, 0])

                        # tListSorted.append([Npres/freq,2])

                        ###########################################################
                        # event-based integration
                        (timesAbove[r, 0], timesAbove[r, 1]) = self.eventBasedIntegration(tListSorted)
                        # print r
                return (average(timesAbove[:,0]) , average(timesAbove[:,1]) )


        ###############################################################################
        # (timeD,timeP) = tat.spikePairFrequencyNonlinear(DeltaTStart,DeltaTEnd,D,frequency)
        def spikePairFrequencyTMM(self, deltaT, freq, Npres, U, tauFac, tauDep ):
                tStart = 0.1  # start time at 100 ms

                #Npres = Npres * 12
                #timeAbove = zeros((1, 2))

                tD = 0.
                tP = 0.
                #random.seed(7)
                tPre  = arange(Npres) / freq + tStart
                tPost = tPre + deltaT

                tAll = zeros((2 * Npres, 3))

                cpre = zeros(Npres)
                u   = zeros(Npres)
                x   = zeros(Npres)
                if U != 0:
                        u[0] = U
                        x[0] = 1.
                        cpre[0] = u[0] * x[0]
                        for i in range(1, len(tPre[1:])):
                                u[i] = u[i - 1] * exp(-1./(freq*tauFac)) + U * (1. - u[i - 1] * exp(-1./(freq*tauFac)))
                                x[i] = 1. - (1. - (x[i - 1] - u[i - 1] * x[i - 1])) * exp(
                                        -1./(freq*tauDep))
                                cpre[i] = u[i] * x[i]
                        cpre *= self.Cpre
                        #cpre[0] = 1.
                        #for i in range(1, Npres):
                        #    cpre[i] = 1. - (1. - (cpre[i-1] - U*cpre[i-1]))*exp(-1./(freq * tauRec))
                        #cpre *= U*self.Cpre
                else:
                        cpre[:] = self.Cpre

                #print cpre
                tAll[:, 0] = hstack((tPre, tPost))
                tAll[:, 1] = hstack((zeros(Npres), ones(Npres)))
                tAll[:, 2] = hstack((cpre,repeat(self.Cpost,Npres)))
                tList = tAll.tolist()

                tListSorted = sorted(tList, key=lambda tList: tList[0])

                tListSorted.append([Npres / freq + tStart, 2,0])

                ###########################################################
                # event-based integration
                (tD, tP) = self.eventBasedIntegration(tListSorted)

                return (tD , tP )

        ###############################################################################
        # (timeD,timeP) = tat.spikePairFrequencyNonlinear(DeltaTStart,DeltaTEnd,D,frequency)
        def spikePairStochasticFrequencySTPDeterministic(self, DeltaTStart, DeltaTEnd, freq, Npres, tauRec, U):
            tStart = 0.1  # start time at 100 ms

            # Npres = Npres * 12
            # timeAbove = zeros((1, 2))
            Naverages = 100.
            tD = 0.
            tP = 0.
            for i in arange(Naverages):
                # random.seed(7)
                tPre = arange(Npres) / freq + tStart + (DeltaTStart + rand(Npres) * (DeltaTEnd - DeltaTStart))
                tPost = tPre + (DeltaTStart + rand(Npres) * (DeltaTEnd - DeltaTStart))

                # tPre  = arange(Npres) / freq + tStart
                # tPost = tPre + deltaT

                tAll = zeros((2 * Npres, 3))

                cpre = zeros(Npres)
                cpre[0] = self.Cpre * U
                for i in range(1, Npres):
                    cpre[i] = cpre[i - 1] * (1. - U * exp(-(tPre[i] - tPre[i - 1]) / (tauRec)))
                tAll[:, 0] = hstack((tPre, tPost))
                tAll[:, 1] = hstack((zeros(Npres), ones(Npres)))
                tAll[:, 2] = hstack((cpre, repeat(self.Cpost, Npres)))
                tList = tAll.tolist()

                tListSorted = sorted(tList, key=lambda tList: tList[0])

                tListSorted.append([Npres / freq + tStart, 2])

                ###########################################################
                # event-based integration
                (tDTemp, tPTemp) = self.eventBasedIntegration(tListSorted)

                tD += tDTemp
                tP += tPTemp

            return (tD / Naverages, tP / Naverages)



        ###############################################################################
        # (timeD,timeP) = tat.spikePairFrequencyNonlinear(DeltaTStart,DeltaTEnd,D,frequency)
        def eventBasedIntegration(self, tListSorted):

                # event-based integration
                ca = []
                # CaTotal, CaPre, CaPost, time
                ca.append([0.,0.,0.,0.])
                pre = 0
                post = 0
                tD = 0.
                tP = 0.
                for i in tListSorted:
                        #
                        caTotOld    = ca[-1][0]
                        caPreOld    = ca[-1][1]
                        caPostOld   = ca[-1][2]
                        tOld        = ca[-1][3]
                        #caTotTemp  = caTotOld*exp(-(i[0]-tOld)/self.tauCa)
                        caPreTemp  = caPreOld*exp(-(i[0]-tOld)/self.tauCa)
                        caPostTemp = caPostOld*exp(-(i[0]-tOld)/self.tauCa)
                        caTotTemp  = caPreTemp + caPostTemp
                        if caTotOld > self.thetaD:
                                if caTotTemp > self.thetaD:
                                        tD += i[0]-tOld
                                else:
                                        tD += (self.tauCa)*log(caTotOld/self.thetaD)
                        if caTotOld > self.thetaP:
                                if caTotTemp > self.thetaP:
                                        tP += i[0]-tOld
                                else:
                                        tP += (self.tauCa)*log(caTotOld/self.thetaP)
                        # postsynaptic spike
                        if i[1] == 1:
                                caPostTemp += i[2] + self.eta*caPreTemp
                                post+=1
                        # presynaptic spike
                        if i[1] == 0:
                                caPreTemp += i[2] #self.Cpre
                                pre+=1
                        caTotTemp = caPreTemp + caPostTemp
                        ca.append([caTotTemp,caPreTemp,caPostTemp,i[0]])
                        #
                        #pdb.set_trace()

                return (tD, tP)

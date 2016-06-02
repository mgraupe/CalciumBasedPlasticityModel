from scipy import *
from numpy import *
import pdb
import sys
import commands

class timeAboveThreshold():
        ''' 
            class to calculate the fraction of time \alpha the calcium trace spends above threshold 
        '''
        ###############################################################################
        def __init__(self, tauCa, Cpre, Cpost, thetaD, thetaP, nonlinear=1.):
                self.tauCa = tauCa
                self.Cpre = Cpre
                self.Cpost = Cpost
                self.thetaD  = thetaD
                self.thetaP  = thetaP
                # determine eta based on nonlinearity factor and amplitudes
                self.eta = (nonlinear*(self.Cpost + self.Cpre) - self.Cpost)/self.Cpre - 1.
                
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
                                        print A, B, C, D, Ct, frequency, deltaT 
                                        print "post-pre : Problem in spikePairFrequency!"
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
                                        print E, F, G, H, Ct, frequency, deltaT
                                        print "pre-post : Problem in spikePairFrequency! " 
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
                                        print A, B, C, D, Ct, frequency, deltaT 
                                        print "post-pre : Problem in spikePairFrequency!"
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
                                        print E, F, G, H, Ct, frequency, deltaT
                                        print "pre-post : Problem in spikePairFrequency! " 
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
                                        print "post-post-pre : Problem !"
                                        print deltaT, A, (A+self.Cpost), M, N, O, P, Int
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
                                        print "post-pre-post : Problem !"
                                        print deltaT, A, Q, R, (A+self.Cpost), S, T, Int 
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
                                        print "pre-post-post : Problem !"
                                        print deltaT, A, (A+self.Cpre), U, V, W, X, Int
                                        sys.exit(1)
                        else :
                                print "Problem in preSpikePostPair routine!"
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
                
                # the first argument calcium ampliutde has to be smaller than the second
                if self.Cpre>self.Cpost:
                        arguments = str(deltaT) + ' ' + str(self.tauCa) + ' ' + str(self.Cpost) + ' ' + str(self.Cpre) + ' ' + str(self.thetaD) + ' ' + str(self.thetaP) + ' ' + str(preRate) + ' ' + str(postRate) + ' ' + str(ppp) + ' ' + str(deltaCa)
                else:
                        arguments = str(deltaT) + ' ' + str(self.tauCa) + ' ' + str(self.Cpre) + ' ' + str(self.Cpost) + ' ' + str(self.thetaD) + ' ' + str(self.thetaP) + ' ' + str(preRate) + ' ' + str(postRate) + ' ' + str(ppp) + ' ' + str(deltaCa)
    
                print arguments
                (out,err) = commands.getstatusoutput('./timeAboveThreshold/poissonPairs_timeAboveThreshold ' + arguments) 
                alphaD = float(err.split()[0])
                alphaP = float(err.split()[1])
                
                return (alphaD,alphaP)
            
                



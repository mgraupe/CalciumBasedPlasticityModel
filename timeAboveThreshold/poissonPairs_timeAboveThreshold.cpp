/*************************************************************************
file:         PoissonPairs_timeAboveThreshold.cpp
author:       Michael Graupner 
mail:         michael.graupner@parisdescartes.fr
version:      1.0
Time-stamp:   <15/07/20 graupner>
---------------------------------------------------------------------

DESCRIPTION: The time the calcium trace spends above depression and potentiation thresholds is calculated. The stimulation protocol consits of spike-pairs with fixed \Delta t but irregular occurrence. The presynaptic neurons fires with rate \nu_pre and the postsynaptic neurons with \nu_post. The probability that a presynaptic spike is followed by a postsynaptic spike with delay \Delta t is given by 'p'. Both, pre- and postsynaptic neurons are firing following Poisson statistics. 

Model as defined in Graupner&Brunel, 2012 PNAS. 

    
ARGUMENTS:
1: time lag of pre-post spike-pair
2: time constant of the intracellular calcium
3: presynaptically induced calcium amplitude 
4: postsynaptically induced calcium amplitude
5: depression threshold
6: potentiation threshold
7: firing rate of the presynaptic neuron
8: firing rate of the postsynaptic neuron
9: proability that pres-spike is followed by post-spike with delay \Delta t

RETURNS:
timeAboveDepressionThreshold
timeAbovePotentiaionThreshold


**************************************************************************/

#include <cmath>
#include <map>
#include <string>
#include <sstream>
#include <fstream>    
#include <iostream>  
#include <algorithm>
#include <cstring>

using namespace std;

// ------------------- global variables---------------------------------------

const double Ca_end = 8.;
//const double delta_Ca = 0.0001; //0.0001; // 0.0001; 


//--------------------analytical part of integral------------------------------

double analytical_integral (double b, long b_steps, long i, double delta_Ca, double rate) {
	//
	double integral = 0.;
	
	for ( long k = 0 ; k <= (i - b_steps) ; k++ ) {
		if ( rate >= 1. ) {
			integral += delta_Ca*(1./(pow((b + delta_Ca*(double(k))),rate)) - 1./(pow(b,rate)))*(pow((delta_Ca*(double(k))),(rate - 1.)));
		}
		else {
			if ( !k == 0 ) {
				integral += delta_Ca*(1./(pow((b + delta_Ca*(double(k))),rate)) - 1./(pow(b,rate)))/(pow((delta_Ca*(double(k))),(1. - rate)));
			}
		}		
	}
	integral += pow((delta_Ca*(double(i)) - b),rate)/(rate*pow(b,rate));
	
	return integral;
}

void time_above_treshold(double * time_t, double rate_pre, double rate_post, double ppp, double delta_tt, double tauCa, double theta_d, double theta_p, double C_pre, double C_post, double delta_Ca){

  
	
	double Cpre, Cpost, rate_1, rate_2, delta_t, p, tau_ca;
	double c1, c2, c3, c4, b1, b2, b3, b4;
	
    
    // auxillary variables
    double Ctresh[2], rate;
    double threshold, time_above, time_above_a;
	long steps, threshold_steps; 
	double* rho, * integral1,* integral2,* integral3,* integral4; 
    long b1_steps, b2_steps, b3_steps, b4_steps;
	double normalization, time_below;
 	
    //if (write_output) {
    //ofstream amp_f("amplitude_distribution.dat");
    //ofstream thresh_f("time_above_threshold.dat");
    //}

    p          = ppp; 
	rate_1     = rate_pre; 
	rate_2     = rate_post - p*rate_pre; 
	tau_ca     = tauCa; 
	Cpre       = C_pre;
	Cpost      = C_post;
    Ctresh[0]  = theta_d; 
	Ctresh[1]  = theta_p; 
	delta_t    = delta_tt; 
    
    //cout << p << " " << rate_1 << " " << rate_2 << " " << tau_ca << " " << Cpre << " " << Cpost << " " << delta_t << endl;
	
	 
	steps = (long)(Ca_end/delta_Ca  + 0.5);
	rho = new double[steps+1];
	integral1 = new double[steps+1];
	integral2 = new double[steps+1];
	integral3 = new double[steps+1];
	integral4 = new double[steps+1];
	//convolution = new double[steps+1];
	for (long i= 0; i<(steps+1) ; i++) {
		rho[i] = 0.;
		integral1[i] = 0.;
		integral2[i] = 0.;
		integral3[i] = 0.;
		integral4[i] = 0.;
		//convolution[i] = 0.;
	}
	

	//delta_Ca = delta_Ca/stretch;
	if (delta_t >= 0.) {
		b1 = Cpre*exp(-delta_t/tau_ca);
		b2 = Cpre;
		b3 = Cpost;
		b4 = Cpre*exp(-delta_t/tau_ca) + Cpost;
		
		c1 = p*rate_1*tau_ca;
		c2 = -rate_1*tau_ca;
		c3 = -rate_2*tau_ca;
		c4 = -p*rate_1*tau_ca;
	}
	else if (delta_t<0. && (Cpre + Cpost*exp(delta_t/tau_ca)) >= Cpost && Cpost*exp(delta_t/tau_ca) < Cpre) {
		
		b1 = Cpost*exp(delta_t/tau_ca);
		b2 = Cpre;
		b3 = Cpost;
		b4 = Cpost*exp(delta_t/tau_ca) + Cpre;
		
		c1 = p*rate_1*tau_ca;
		c2 = -(1.-p)*rate_1*tau_ca;
		c3 = -(p*rate_1+rate_2)*tau_ca;
		c4 = -p*rate_1*tau_ca;
	}
	else if (delta_t<0. && (Cpre + Cpost*exp(delta_t/tau_ca)) >= Cpost && Cpost*exp(delta_t/tau_ca) >= Cpre) {
		
		b1 = Cpre;
		b2 = Cpost*exp(delta_t/tau_ca);
		b3 = Cpost;
		b4 = Cpost*exp(delta_t/tau_ca) + Cpre;
		
		c1 = -(1.-p)*rate_1*tau_ca;
		c2 = p*rate_1*tau_ca;
		c3 = -(p*rate_1+rate_2)*tau_ca;
		c4 = -p*rate_1*tau_ca;
	}
	else if (delta_t<0. && (Cpre + Cpost*exp(delta_t/tau_ca)) < Cpost && Cpost*exp(delta_t/tau_ca) < Cpre) {
		
		b1 = Cpost*exp(delta_t/tau_ca);
		b2 = Cpre;
		b3 = Cpost*exp(delta_t/tau_ca) + Cpre;
		b4 = Cpost;
		
		c1 = p*rate_1*tau_ca;
		c2 = -(1.-p)*rate_1*tau_ca;
		c3 = -p*rate_1*tau_ca;
		c4 = -(p*rate_1+rate_2)*tau_ca;
	}
	else if (delta_t<0. && (Cpre + Cpost*exp(delta_t/tau_ca)) < Cpost && Cpost*exp(delta_t/tau_ca) >= Cpre) {
		
		b1 = Cpre;
		b2 = Cpost*exp(delta_t/tau_ca);
		b3 = Cpost*exp(delta_t/tau_ca) + Cpre;
		b4 = Cpost;
		
		c1 = -(1.-p)*rate_1*tau_ca; 
		c2 = p*rate_1*tau_ca;
		c3 = -p*rate_1*tau_ca;
		c4 = -(p*rate_1+rate_2)*tau_ca;
	}
	else {
		cout << "problem in amplitude calculations" << endl;
		exit(1);
	}
	
	if ( !(b1 <= b2 && b2 <= b3 && b3 <= b4) ) {
		cout << "Wrong amplitude relations" << endl;
		cout << delta_t << "  " << b1 << "  " << b2 << "  " << b3 << "  " << b4 << endl;
		exit(1);
	}
	//c = exp(-rate*EULER)/(exp(gamma(rate))); 
	rate = tau_ca*(rate_1+rate_2);
	b1_steps = long(b1/delta_Ca + 0.5);
	b2_steps = long(b2/delta_Ca + 0.5);
	b3_steps = long(b3/delta_Ca + 0.5);
	b4_steps = long(b4/delta_Ca + 0.5);
	//small_amp  = amplitude;
	//kappa = 1.; // exp(-rate*EULER)/(pow(amplitude,rate)*exp(gamma(rate))); 
	//cout << "calculating : amplitudes = " << b1 << " " << b2 << " " << b3 << " " << b4 << " for pre- and post rates = " << rate_1 << "   " << rate_2 << "  " << rate << endl;
	
	for( long i=0;i<=steps;i++){

		//multiple = i/unit_steps;

		if ( i <= b1_steps ) {
			//rho[i] = exp(-rate*EULER)*pow(((double)i)*delta_Ca,(rate-1.))/gamma(rate);
			//cout << i << endl;
			if (rate >= 1.) {
				rho[i] = pow((((double)i)*delta_Ca),(rate-1.));
			}
			else {
				if ( i == 0 )
					rho[i] = 1./delta_Ca; //exp(-rate*EULER)*pow(delta_Ca,rate)/(exp(gamma(rate))*rate);
				else 
					rho[i] = 1./(pow((((double)i)*delta_Ca),(1.-rate)));
			}
			//cout << "ende" << endl;
		}

		// density function between amplitude and 2*amplitude
		else if (i>b1_steps && i <= b2_steps ) {
			//cout << " ! " << endl;
			//integral = 0.;
			//multiple = i/unit_steps;
			//cout << i << "   " << unit_steps << "   " <<  multiple << " " << integral[i] << endl;
			if (i <= 2*b1_steps ) {
				integral1[i] = analytical_integral(b1, b1_steps, i, delta_Ca,rate);
				rho[i] = pow((((double)i)*delta_Ca),(rate-1.))*(1. + c1*integral1[i]);
			}
			else {
				for ( long k = ((2*b1_steps)+1) ; k <= i ; k++ ) {
					//cout << i << "  " <<  k-unit_steps << endl;
					integral1[i] += delta_Ca*rho[(k - b1_steps)]/(pow((((double)k)*delta_Ca),rate)); 
				}

				integral1[i] += integral1[(2*b1_steps)];
				//cout << integral[i] << "  " << integral[(multiple*unit_steps)] << endl;
				//exit(1);
				rho[i] = pow((((double)i)*delta_Ca),(rate-1.))*(1. + c1*integral1[i]);
				
			}
			
		}
		// density function for amplitude < 2*amplitude
		else if (i>b2_steps && i <= b3_steps ) {	
				//integral = 0.;
			//multiple = i/unit_steps;
			//cout << i << "   " << unit_steps << "   " <<  multiple << "\t" ;
			//exit(1);
			if (i <= 2*b1_steps ) {
				integral1[i] = analytical_integral(b1, b1_steps, i, delta_Ca,rate);
				
			}
			else {
				for ( long k = ((2*b1_steps)+1) ; k <= i ; k++ ) {
					//cout << i << "  " << k - b1_steps << "  " << k - b2_steps << endl;
					integral1[i] += delta_Ca*rho[(k - b1_steps)]/(pow((((double)k)*delta_Ca),rate));
				}
				//cout << integral1[i] << "  " << integral2[i] << "  " ;
				integral1[i] += integral1[(2*b1_steps)];
				//cout << integral1[i] << "  " << integral2[i] << endl;
			}
			if (i <= (b1_steps + b2_steps) ) {
				integral2[i] = analytical_integral(b2, b2_steps, i, delta_Ca,rate);
			}
			else {
				for ( long k = ((b1_steps+b2_steps)+1) ; k <= i ; k++ ) {
					integral2[i] += delta_Ca*rho[(k - b2_steps)]/(pow((((double)k)*delta_Ca),rate));
				}
				integral2[i] += integral2[(b1_steps + b2_steps)];
			}
			
			rho[i] = pow((((double)i)*delta_Ca),(rate-1.))*(1. + c1*integral1[i] + c2*integral2[i]);
			
		}
		else if (i>b3_steps && i <= b4_steps ) {	
				//integral = 0.;
			//multiple = i/unit_steps;
			//cout << i << "   " << unit_steps << "   " <<  multiple << "\t" ;
			//exit(1);
			// first integral 
			if (i <= 2*b1_steps ) {
				integral1[i] = analytical_integral(b1, b1_steps, i, delta_Ca,rate);
			}
			else {
				for ( long k = ((2*b1_steps)+1) ; k <= i ; k++ ) {
					integral1[i] += delta_Ca*rho[(k - b1_steps)]/(pow((((double)k)*delta_Ca),rate));
				}
				integral1[i] += integral1[(2*b1_steps)];
			}
			// second integral 
			if (i <= (b1_steps + b2_steps) ) {
				integral2[i] = analytical_integral(b2, b2_steps, i, delta_Ca,rate);
			}
			else {
				for ( long k = ((b1_steps+b2_steps)+1) ; k <= i ; k++ ) {
					integral2[i] += delta_Ca*rho[(k - b2_steps)]/(pow((((double)k)*delta_Ca),rate));
				}
				integral2[i] += integral2[(b1_steps + b2_steps)];
			}
			// third integral 
			if (i <= (b1_steps + b3_steps) ) {
				integral3[i] = analytical_integral(b3, b3_steps, i, delta_Ca,rate);
			}
			else {
				for ( long k = ((b1_steps+b3_steps)+1) ; k <= i ; k++ ) {
					integral3[i] += delta_Ca*rho[(k - b3_steps)]/(pow((((double)k)*delta_Ca),rate));
				}
				integral3[i] += integral3[(b1_steps + b3_steps)];
			}

			rho[i] = pow((((double)i)*delta_Ca),(rate-1.))*(1. + c1*integral1[i] + c2*integral2[i] + c3*integral3[i]);
			
		}
		else if (i>b4_steps) {	
				//integral = 0.;
			//multiple = i/unit_steps;
			//cout << i << "   " << unit_steps << "   " <<  multiple << "\t" ;
			//exit(1);
			if (i <= 2*b1_steps ) {
				integral1[i] = analytical_integral(b1, b1_steps, i, delta_Ca,rate);
			}
			else {
				for ( long k = ((2*b1_steps)+1) ; k <= i ; k++ ) {
					integral1[i] += delta_Ca*rho[(k - b1_steps)]/(pow((((double)k)*delta_Ca),rate));
				}
				integral1[i] += integral1[(2*b1_steps)];
			}
			// second integral 
			if (i <= (b1_steps + b2_steps) ) {
				integral2[i] = analytical_integral(b2, b2_steps, i, delta_Ca,rate);
			}
			else {
				for ( long k = ((b1_steps+b2_steps)+1) ; k <= i ; k++ ) {
					integral2[i] += delta_Ca*rho[(k - b2_steps)]/(pow((((double)k)*delta_Ca),rate));
				}
				integral2[i] += integral2[(b1_steps + b2_steps)];
			}
			// third integral 
			if (i <= (b1_steps + b3_steps) ) {
				integral3[i] = analytical_integral(b3, b3_steps, i, delta_Ca,rate);
			}
			else {
				for ( long k = ((b1_steps+b3_steps)+1) ; k <= i ; k++ ) {
					integral3[i] += delta_Ca*rho[(k - b3_steps)]/(pow((((double)k)*delta_Ca),rate));
				}
				integral3[i] += integral3[(b1_steps + b3_steps)];
			}
			// fourth integral
			if (i <= (b1_steps + b4_steps) ) {
				integral4[i] = analytical_integral(b4, b4_steps, i, delta_Ca,rate);
			}
			else {
				for ( long k = ((b1_steps+b4_steps)+1) ; k <= i ; k++ ) {
					integral4[i] += delta_Ca*rho[(k - b4_steps)]/(pow((((double)k)*delta_Ca),rate));
				}
				integral4[i] += integral4[(b1_steps + b4_steps)];
			}
			
			rho[i] = pow((((double)i)*delta_Ca),(rate-1.))*(1. + c1*integral1[i] + c2*integral2[i] + c3*integral3[i] + c4*integral4[i]);
			
		}
		// the density function of the convolution can be calculated until I = 1 only, sigh!
		//else if ( i > unit_steps && j == 2 ) {
		//	break;
		//}
		else {
			cout << "Problem in density function integration ! " << endl;
			exit(1);
		}
		//if (write_output) {
        //amp_f << i << "\t" << ((double)i)*delta_Ca << "\t" << rho[i] << "\t" << integral1[i] << "\t" << integral2[i] << "\t" << integral3[i] << "\t" << integral4[i] << endl;
        //}
				
	}

	for (int j = 0 ; j < 2 ; j++ ) { 
		threshold = Ctresh[j];
		//steps = (long)(Ca_end/delta_Ca);
		threshold_steps = (threshold/delta_Ca + 0.5);
		//double start = (0.5/delta_Ca + 0.5);
		//cout << threshold_steps << endl;
		time_above = time_above_a = time_below = 0.;
		normalization = 0.;
		if ( threshold_steps < b1_steps ) {
			time_below = pow(threshold,rate)/(rate);
		}
		else  {
			time_below = pow(b1,rate)/(rate);
		}
		
		//cout << time_above << endl;
		for (long i = (b1_steps+1) ; i <= threshold_steps ; i++ ) {
			time_below += rho[i]*delta_Ca;
		}
		normalization = time_below;
		for (long i = (threshold_steps+1) ; i <= steps ; i++ ) {
 			normalization += rho[i]*delta_Ca;
 		}
		
		time_above = (1. - time_below/normalization);
		
 		// this methods is used for calculating \bar{rho}
		for (long i = threshold_steps ; i <= steps ; i++ ) {
 			time_above_a += delta_Ca*rho[i]/normalization;
 		}
		//cout << keep << "  " << time_above_3 << "  " << (keep - time_above_3) << endl;
		//keep = (1. - time_above_3);
		//cout << threshold << "  " << time_above << "  " << time_above_a << "  " << time_below << "  " << normalization <<  endl;
		time_t[j] = time_above_a;
        //if (write_output) {
        //thresh_f << threshold << "  " << time_above << "  " << time_above_a << "  " << time_below << "  " << normalization << endl;
        //}
	}
    //amp_f.close();
    //thresh_f.close();
    
	//delete[] rho,integral1,integral2,integral3,integral4, convolution;
	
    return;

}

// ------------------- main program ------------------------------------------
 
int main(int argc, char* argv[]){
    // time_above_treshold(double * time_t, double rate_pre, double rate_post, double ppp, double delta_tt, double theta_d, double theta_p, double C_pre, double C_post, bool write_output){
    
    double deltaT     = atof(argv[1]);
    double tauCa      = atof(argv[2]);
    double Cpre       = atof(argv[3]);
    double Cpost      = atof(argv[4]);
    double theta_d    = atof(argv[5]);
    double theta_p    = atof(argv[6]);
    double rate_pre   = atof(argv[7]);
    double rate_post  = atof(argv[8]);
    double p          = atof(argv[9]);
    double delta_Ca   = atof(argv[10]);
    
    // AUX variables
    double timeAboveThres[2];
    double timeAboveDepressionThreshold, timeAbovePotentiationThreshold;
    
    ofstream integral_f("timeAboveThreshold.dat");
    
    time_above_treshold(timeAboveThres,rate_pre,rate_post,p,deltaT,tauCa,theta_d,theta_p,Cpre,Cpost,delta_Ca);

    timeAboveDepressionThreshold   = timeAboveThres[0];
    timeAbovePotentiationThreshold = timeAboveThres[1];
    
    integral_f << rate_pre << "\t" << rate_post << "\t" <<  p  << "\t" << deltaT << "\t" << Cpre << "\t" << Cpost  << "\t" << theta_d << "\t" << theta_p << "\t" << timeAboveDepressionThreshold << "\t" << timeAbovePotentiationThreshold << endl;
    
    cout << timeAboveDepressionThreshold << "\t" << timeAbovePotentiationThreshold << endl;
    
    integral_f.close();
    
    return 0;
}



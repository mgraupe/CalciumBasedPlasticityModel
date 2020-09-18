/*************************************************************************
file:         camkmotifscan.cpp
author:       Michael Graupner 
mail:         michael.graupner@univ-paris5.fr
version:      0.1
last-change:  05.12.2005 mg
---------------------------------------------------------------------
description:  program simulates the activity of CaMKII in response to 
              a scan of stimulation motifs (at the moment only STDP)
**************************************************************************/  

#include <iostream>
#include <fstream>
#include <cmath>     
#include <string>

#include "nrutils.hpp"
#include "parameter.hpp"
#include "motif.hpp"
#include "compose.hpp"

using namespace std;
int seed;
int idum;

//-main----------------------------------------------------------------------------------

int main ()
{
  //free(idum);
  //int test = -8;
  //cout << test << "  " << ran1(test) << endl;
  Parameter par("camkmotifscan.par");
  seed = par.get_param<int>("seed");
  

  using String::compose;

  if ( !seed ) 
   	idum = - (time(NULL));
  else 
    idum = 7;

  //cout << gasdev(idum) << endl;
  //exit(1);
  //time_t Tval = 0;
  //Tval = time(NULL);
  //printf("%d \n",Tval);
  double scan_start, scan_end, delta_scan, scan_param, t_stim_start;
  int scan_steps;
  double delay;
  double cond_up, init_down, syn_strength;

  double rho_unstable = par.get_param<double>("rho_unstable");
  ofstream final_state_f("final_camkII_state.dat");
  ofstream end_stimulation_state_f("end_stimulation_camkII_state.dat");
  ofstream dynamics_f;
  
  scan_start  = par.get_param<double>("start_delta_t_1");
  scan_end    = par.get_param<double>("end_delta_1");
  scan_steps  = par.get_param<int>("scan_steps_1");
  t_stim_start= par.get_param<int>("t_start_1");

  init_down   = par.get_param<double>("init_down");
  cond_up     = par.get_param<double>("cond_up");
  
  delta_scan=double((scan_end-scan_start)/scan_steps);

  delay      = par.get_param<double>("delay");

  const double t_start   = par.get_param<double>("t_start");
  const double h         = par.get_param<double>("dt");
  //const long write_steps = par.get_param<long>("write_steps");
  const double epsilon   = par.get_param<double>("epsilon");
  double dt, dtdid, dtnext;
  //long next_write; 
  
  string file_name;
  const char* file_dyn;

  unsigned int runs = par.get_param<int>("runs");
  double average1, average2, average3, average4, average5; 
  double sigma1, sigma2, sigma3a, sigma3b, sigma4, sigma5;
  double *** end_stimulation_state, *** final_state; 
  end_stimulation_state = new double**[2];
  final_state = new double**[2];
  for (int i=0 ; i<2 ; i++ ) {
	end_stimulation_state[i] = new double*[(scan_steps+1)];
	final_state[i] = new double*[(scan_steps+1)];
	for ( int k=0 ; k< (scan_steps+1) ; k++ ) {
		end_stimulation_state[i][k] = new double[runs];
		final_state[i][k] = new double[runs];
		for ( int l=0 ; l < runs ; l++ ) {
			end_stimulation_state[i][k][l] = 0.;
			final_state[i][k][l] = 0.;
		}
	}
   }
  
   //double end_stimulation_state[2][scan_steps+1][runs], final_state[2][scan_steps+1][runs];

	

  	// loop for up-switching "k=0" and down-shifting "k=1"
  	for (int k=0; k<2; k++) {  
    	scan_param=scan_start;
    
    	// loop for different delta-t between pre- and post-spike
    	for (long j=0; j<=scan_steps; ++j) {
			//final_state << scan_param << "\t" << k;
			//end_stimulation_state << scan_param << "\t" << k;
			// loop over "n" runs under the same conditiions but different random seeds
			for (unsigned int n=0 ; n<runs ; n++) {
				double t   = t_start;
				//next_write = write_steps;
				dt         = h;
				dtdid      = h;
				dtnext     = h;
				pre_post_spikes pre_post(par,(scan_param-delay),k);

				//file_name=compose("rho_dyn_%1.dat",n);
				//file_dyn=file_name.c_str();
				//dynamics_f.open(file_dyn);

				//dynamics_f << t << "\t" << pre_post.yt[0] << endl;
				// integration
				for (long i=1;; ++i) {
					
                    

					pre_post.rkqs(t, dt, epsilon, dtdid, dtnext);
					//dynamics_f << t << "\t" << pre_post.yt[0] << endl;
					// check if pre-stimuli steady-state is reached
					//if (pre_post.stop && pre_post.stage==1) {
					//  t=0.;  
						//next_write = (i + write_steps);
					//  pre_post.protocol=true;
						//pre_post.write_files(t);   
					//}
					// end of stimulation protocol
					if (t>pre_post.t_stim_end) {
						dt         = 100*h;
						dtdid      = 100*h;
						dtnext     = 100*h;
						
					}
					// convergence to final steady-state after stimulation protocol
					if (pre_post.stop && pre_post.stage==2) {
						end_stimulation_state[k][j][n] = pre_post.phossum;
						//cout << scan_param << "\t" << pre_post.ca_max << "\t" << pre_post.buff_max << endl;
						break;
	  				}
					//else if (i>=next_write || pre_post.pre_spike || pre_post.post_spike) {
						//pre_post.write_files(t); 
						//next_write += write_steps;
						//}
	  				dt = dtnext;
				}
				//dynamics_f.close();
				final_state[k][j][n] = pre_post.phossum ;
				//cout << scan_param << "   " << k << "  " << pre_post.phossum << endl;
				//end_state[n] = pre_post.phossum;
				pre_post.free();
				//cout << scan_param << "\t" << pre_post.time_dephos << "\t" << pre_post.time_phos << "\t" << endl;
      		}
			//final_state << "\t" << average1 << endl;
			//end_stimulation_state << "\t" << average2 << "\t" << average3 << endl;
			scan_param+=delta_scan;
			
    	}
    	//final_state << endl << endl << endl;
  	}

	// loops to calculate the average transitions results
	scan_param=scan_start;
	for (long j=0; j<=scan_steps; ++j) { 
		average1 = average2 = average3 = average4 = average5 = 0. ;
		end_stimulation_state_f << scan_param << "\t";
		final_state_f << scan_param << "\t" ;
		for (unsigned int n=0 ; n<runs ; n++) {     
			average1 += end_stimulation_state[0][j][n]/(double(runs));
			average2 += end_stimulation_state[1][j][n]/(double(runs));
			//average4 += final_state[0][j][n]/((double)runs);
			//average5 += final_state[1][j][n]/((double)runs);
	
			if ( final_state[0][j][n] > 0.5 ) {
				average3 += 1./((double)runs);
				average4 += 1./((double)runs);
			}
			//if ( final_state[1][j][n] > 0.5 )
			//    average3 += 1./((double)runs) ;
			//else if (final_state[0][j][n] < 0.5 && final_state[1][j][n] > 0.5)
			//    average3 += 0.;
			if ( final_state[1][j][n] <= 0.5 ) {
				average3 += -1./((double)runs);
				average5 += -1./((double)runs); 
			}
			//cout << scan_param << "  " <<  final_state[0][j][n] << "  " << final_state[1][j][n] << endl;
    	}
	
		sigma1 = sigma2 = sigma3a = sigma3b = sigma4 = 0.;
		for (unsigned int n=0 ; n<runs ; n++) {     
			sigma1 += pow((end_stimulation_state[0][j][n] - average1),2);
			sigma2 += pow((end_stimulation_state[1][j][n] - average2),2);
			//sigma3a += pow((final_state[0][j][n] - average4),2);
			//sigma3b += pow((final_state[1][j][n] - average5),2);
//          if ( final_state[0][j][n] > 0.5 ) {
// 			sigma3a += pow((1. - average3),2);
// 			sigma3b += pow((final_state[0][j][n] - average3),2);
// 			//cout << n << "  up  " << average3 << "  " << final_state[0][j][n] << "  " << sigma3a << endl;
// 		}
// 		else {
// 			sigma3a += pow((average3),2);
// 			sigma3b += pow((final_state[0][j][n] - average3),2);
// 			//cout << n << "  down  " << average3 << "  " << final_state[0][j][n] << "  " << sigma3a << endl;
// 		}
// 			
// 		//if ( final_state[1][j][n] > 0.5 )
//         //    average3 += 1./((double)runs) ;
//         //else if (final_state[0][j][n] < 0.5 && final_state[1][j][n] > 0.5)
//         //    average3 += 0.;
// 		// in case of transtion
//         if ( final_state[1][j][n] < 0.5 ) {
// 			sigma3a += pow((1. - average3),2);
//             sigma3b += pow((final_state[1][j][n] + 1. - average3),2);
// 			//cout << n << "  down  " << average3 << "  " << final_state[1][j][n] << "  " << sigma3a << endl;
// 		}
// 		// in case no transition occurred
// 		else {
// 			sigma3a += pow((average3),2);
//             sigma3b += pow((final_state[1][j][n] - 1. - average3),2);
// 			//cout << n << "  up  " << average3 << "  " << final_state[1][j][n] << "  " << sigma3a << endl;
// 		}
			
		}
		sigma1 = sqrt(sigma1/((double)runs));
		sigma2 = sqrt(sigma2/((double)runs));
		sigma3a = sqrt(average4*(1.-average4)/((double)runs));
		sigma3b = sqrt(fabs(average5)*(1.-fabs(average5))/((double)runs));
		sigma4  = sqrt((average4*(1.-average4) +  fabs(average5)*(1.-fabs(average5))));
		sigma4  = sigma4/(sqrt((double)runs));

		//cout << average4 << "  " << average5 << "  " << fabs(average5) << "   " << sigma4 << endl;
		end_stimulation_state_f << average1 << "\t" << average2 << "\t" << sigma1 << "\t" << sigma2  << endl;

		sigma5 = (cond_up - 1.)*(sqrt((init_down*average4*(1.-average4) + (1.-init_down)*fabs(average5)*(1.-fabs(average5)))/((double)runs)));

		syn_strength = (cond_up + (cond_up - 1.)*(init_down*(fabs(average4) -1.) - (1.-init_down)*fabs(average5)))/((1.-init_down)*cond_up + init_down);
		//  ((init_down*(1.-average4) + (1.-init_down)*fabs(average5)) + (init_down*average4 + (1.-init_down)*(1.-fabs(average5)))*cond_up)/(init_down + (1.-init_down)*cond_up);

		final_state_f << average3 << "\t" << (average4 + average5) << "\t" <<  sigma4  << "\t" << syn_strength << "\t" <<  sigma5 << "\t" << average4  << "\t" << sigma3a << "\t" << average5 << "\t" << sigma3b << "\t" << (sigma3a + sigma3b) <<  endl;

		scan_param+=delta_scan;
	
 	}
  	final_state_f.close();
  	end_stimulation_state_f.close();
  	return (0);
}

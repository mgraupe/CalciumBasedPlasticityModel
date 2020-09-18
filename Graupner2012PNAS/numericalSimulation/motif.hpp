/*********************************************************************
file:         motif.hpp                                             
authors:      Michael Graupner                                        
mail:         michael.graupner@univ-paris5.fr                        
version:      0.1                                                    
last change:  23.02.2005 mg                                           
--------------------------------------------------------------------- 
description:  class which describes the dynamics of CaMKII 
              phosphorylation    
*********************************************************************/

#ifndef motif_H
#define motif_H

#include <fstream>
#include <cmath>
#include <string>

#include "parameter.hpp"
#include "compose.hpp"
#include "nr.hpp"

using namespace NR;
extern int idum;

//-pre_post_spikes-----------------------------------

class pre_post_spikes {
public:
  pre_post_spikes(Parameter&, const double&, const int&);
  ~pre_post_spikes();
  void evaluate_event(const double&, const double&);
  void rkqs(double &, const double, const double, double &, double &);
  void write_files(const double&);
  void free();
  bool stop, protocol, pre_spike, post_spike;
  int stage, success;
  double t_stim_end, phossum, pp1, PP1_rest;
  double *yt;
  double can_int, pka_int, t_integrate;
  double time_phos, time_dephos;  
  double buff_max, ca_max;

private:   
  void functions(double , double [], double []); 
  void rk4(double [], double [], int , double, double, double []);
  void rkck(const double [], const double [], const double, const int, const double, double []);
  void update(double);
  
  double t_start, rate, delta_t, ddt, dt_done;
  double t_compare, t_comp_dist, epsilon_comp;
  double t_pre_spike, t_post_spike, inter_pair, inter_spike, t_duration, t_stim_ap;
  long   n_presentation, post_count, pre_count;
  double C_pre, tau_pre, tau_pre_rise, C_pre_scale, tau_post, C_post;
  double buff, K_buff;
  bool   pre_s, post_s, docking;
  double B0, kon, koff;
  
  double C_pre_mean, C_post_mean, C_pre_quant, C_post_quant, sigma_pre, sigma_post;
  double pre_open_prob, post_open_prob;
  double n_binom_post, n_binom_pre, post_number, pre_number;
  double test, C_pre_drawn;
  bool  draw_pre, draw_post;
  double p_release_0, tau_facilitation_recovery, p_release_old;
  double tau_depression_recovery, dock_threshold, facilitation;
  int seed, n_vesicles, n_release;

  double Ca,ca_rest;
  double t_stim_start, t_initialization, t_simulation;
  double t_write, dt_write;

  double phos, dephos, C, F;
  double tau_rho, rho_up, rho_unstable, alpha, beta, eta, sigma;
  double Ct_phos, Ct_dephos;
  double D, gamma, phossum_old; // phossum_old1, phossum_old2, phossum_old3; 
  int state, scenario;
  int order; 
 
  
  double epsilon;

  int nvar, nvar_total;
  
  double *dydt, *yout, *dy, *yerr, *yscal, *zufall;
  bool *vesicle_docked, *vesicle_docked_old, *vesicle_do;

};

#endif

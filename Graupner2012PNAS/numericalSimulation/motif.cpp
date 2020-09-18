/*********************************************************************
file:         motif.cpp                                             
authors:      Michael Graupner                                        
mail:         michael.graupner@univ-paris5.fr                        
version:      0.1                                                    
last change:  05.12.2005 mg                                           
--------------------------------------------------------------------- 
description:  class which describes the dynamics of the transmembrane
              voltage, the calcium concentration and CaMKII 
              phosphorylation            
*********************************************************************/

#include "motif.hpp"


//--pre_post_spikes------------------------------------------------------------

pre_post_spikes::pre_post_spikes(Parameter& pre_post_spikes_par, const double& scan_delta_t, const int& run) {
 
  using String::compose;
  
  t_start        = pre_post_spikes_par.get_param<double>("t_start_1");
  rate           = pre_post_spikes_par.get_param<double>("rate_1");
  delta_t        = scan_delta_t;
  n_presentation = pre_post_spikes_par.get_param<long>("n_presentation_1");
  epsilon   = pre_post_spikes_par.get_param<double>("epsilon");
  //t_comp_dist    = pre_post_spikes_par.get_param<double>("t_comp_dist");
  dt_done        = 0.;
  //t_initialization = pre_post_spikes_par.get_param<double>("t_initialization");
  t_simulation     = pre_post_spikes_par.get_param<double>("t_simulation");

  ca_rest        = pre_post_spikes_par.get_param<double>("ca_rest");
  C_pre     = pre_post_spikes_par.get_param<double>("C_pre"); 
  tau_pre        = pre_post_spikes_par.get_param<double>("tau_pre");
  //tau_pre_rise   = pre_post_spikes_par.get_param<double>("tau_pre_rise");
  C_post    = pre_post_spikes_par.get_param<double>("C_post"); 
  tau_post       = pre_post_spikes_par.get_param<double>("tau_post"); 
//   buff      = pre_post_spikes_par.get_param<double>("buff");    
//   K_buff    = pre_post_spikes_par.get_param<double>("K_buff"); 

  //C_pre_scale   = 6.5068; //(1./(1./tau_pre_rise - 1./tau_pre))*exp(- log(tau_pre_rise/tau_pre)/(1. - tau_pre/tau_pre_rise)) - exp(-log(tau_pre_rise/tau_pre)/(tau_pre_rise/tau_pre - 1.));
  // calcium noise variables 
  seed          = pre_post_spikes_par.get_param<int>("seed");
  //sigma_ca      = pre_post_spikes_par.get_param<double>("sigma_ca");
  
  // stochastic vesicle release 
//   p_release_0    = pre_post_spikes_par.get_param<double>("p_release_0");
//   tau_facilitation_recovery = pre_post_spikes_par.get_param<double>("tau_facilitation_recovery");
//   facilitation   = pre_post_spikes_par.get_param<double>("facilitation");
//   n_vesicles     = pre_post_spikes_par.get_param<int>("n_vesicles");
//   tau_depression_recovery = pre_post_spikes_par.get_param<double>("tau_depression_recovery");

//   vesicle_docked = new bool[n_vesicles];
//   vesicle_docked_old = new bool[n_vesicles];
//   vesicle_do = new bool[n_vesicles];
//   zufall = new double[n_vesicles];
//   for (int j=0; j<n_vesicles; j++ ) 
//     vesicle_docked[j] = true;
//   n_release = 0;
//   success = 0;
// 
//   sigma_pre      = C_pre_mean*pre_post_spikes_par.get_param<double>("sigma_pre_p");
//   sigma_post      = C_post_mean*pre_post_spikes_par.get_param<double>("sigma_post_p");

  draw_pre = draw_post = true;

  // change seed of the random numbers 
  for (int i=0; i<=seed ; i++) {
    test = NR::gasdev(idum);
  }	
  

  //t_compare    = t_comp_dist;
  
  t_pre_spike  = t_start;
  t_post_spike = t_start + delta_t;
  inter_spike  = delta_t;
  inter_pair   = 1./rate;
  t_duration   = inter_pair*(double(n_presentation-1));
  pre_s        = false;
  post_s       = false;
  pre_spike    = false;
  post_spike   = false;
  pre_count    = 0;
  post_count   = 0; 
  
 
  t_stim_start = 0.; //t_start + delta_t;
  t_stim_end   = t_start + (1./rate)*(double(n_presentation));// t_start + (1./rate)*(double(n_presentation));

 
  
  // CaMKII parameter ==============================================
  tau_rho   = pre_post_spikes_par.get_param<double>("tau_rho");
  rho_up    = pre_post_spikes_par.get_param<double>("rho_up");
  rho_unstable = pre_post_spikes_par.get_param<double>("rho_unstable");
  phos      = pre_post_spikes_par.get_param<double>("phos");
  dephos    = pre_post_spikes_par.get_param<double>("dephos");
  Ct_phos   = pre_post_spikes_par.get_param<double>("Ct_phos");
  Ct_dephos = pre_post_spikes_par.get_param<double>("Ct_dephos");
  sigma     = pre_post_spikes_par.get_param<double>("sigma");

  dt_write  = pre_post_spikes_par.get_param<double>("dt_write");

  state     = run;   
  
// 	B0    = pre_post_spikes_par.get_param<double>("B0");
// 	kon   = pre_post_spikes_par.get_param<double>("kon");
// 	koff  = pre_post_spikes_par.get_param<double>("koff");

  
  Ca = ca_rest;
  
  // general initialization
  // allocating memory 
  nvar_total = 4;
  nvar = nvar_total;

  yt   = new double[nvar];
  dydt = new double[nvar];
  yout = new double[nvar];
  dy   = new double[nvar];
  yerr = new double[nvar];
  yscal= new double[nvar];

  if (state==0)                                // if starting from low-phos. stable fixpoint
    yt[0] = 0.;
  else if (state==1)                           // if starting from highly-phos. stable fixpoint
    yt[0] = rho_up;    
  else {
    cout << "Problem in motif class initialization! " << endl;
    exit(1);
  }
 
  yt[1] = 0.;       
  yt[2] = 0.;
  yt[3] = 0.;
  //yt[4] = 0.;
  //yt[3] = p_release_0;
     
  phossum_old = 0.;     // for convergence check
  stop        = false;
  protocol    = true;
  stage       = 1;
  t_write     = 0.;
  buff_max   = 0.;
  ca_max = 0.;
}

//-----------------------------------------------------------------------------

pre_post_spikes::~pre_post_spikes() {
  
  //delete  [ ] yt;
  //delete  [ ] dydt;
  //delete  [ ] yout;
  //delete  [ ] dy;  
  //delete  [ ] yerr;
  //delete  [ ] yscal;

   
}



//-----------------------------------------------------------------------------

void pre_post_spikes::evaluate_event(const double& t, const double& dt) {
  
  ddt = dt;
  // during stimulation protocol check for occurences of spikes
  if (protocol) { 
//     for (int j=0 ; j<n_vesicles ; j++ ) 
//             vesicle_docked_old[j] = vesicle_docked[j];
// 
//     dock_threshold = dt/tau_depression_recovery;
//     for (int i=0 ; i<n_vesicles ; i++ ) {
//             if ( !vesicle_docked[i] && (dock_threshold > (NR::ran1(idum))) ) {
//                 vesicle_docked[i] = true; 
//                 docking = true;
//             }
//     }
    // presynaptic spike
    if((t+dt)>=t_pre_spike && t_pre_spike<=(t_duration + t_start + 0.5)) {   
        pre_s=true;
//         if (draw_pre) {
//             C_pre_drawn  = C_pre_mean + sigma_pre*(NR::gasdev(idum));
//             for (int i=0 ; i<n_vesicles ; i++ ) 
//                 zufall[i] = (NR::ran1(idum));
//             draw_pre = false;
//         }
//         n_release = 0;
//         for (int i=0 ; i<n_vesicles ; i++ ) {
//             if ( vesicle_docked_old[i] && (yt[3] > zufall[i]) ){
//                 n_release ++;
//                 vesicle_docked[i] = false;
//                 docking = true;                    
//                 //cout << t <<  "  Release ! " << vesicle_docked_old[0] << "->"  << vesicle_docked[0] << "  " << vesicle_docked_old[1] << "->"  << vesicle_docked[1] << "  " << C_pre <<  endl;
//             }
//         }
        //cout << t << "  " << vesicle_docked_old[0] << "->"  << vesicle_docked[0] << "  " << vesicle_docked_old[1] << "->"  << vesicle_docked[1] << "  " << zufall[0] << "  " << zufall[1] <<  endl;
//         if ( (n_release == 0) || (C_pre < 0.))
//             C_pre = 0.;
//         else 
//             C_pre = C_pre_drawn;
        yt[1] += C_pre; //C_nmda;                                  // kick to NMDA calcium
//         p_release_old = yt[3];
//         yt[3] += facilitation*yt[3];
        pre_spike=true;
        //cout << "presynaptic spike at : " << t << endl;
	   //cout << t << "\t" << C_nmda << endl;
    }
    else 
        pre_s=false;

    // postsynaptic spike
    if((t+dt)>=t_post_spike && t_post_spike<=(t_duration + t_start + inter_spike + 0.5) && t_post_spike >= 0.) {  
	 post_s  = true;
// 	 if (draw_post) {
// 		//n_binom_post = NR::bnldev(post_open_prob,post_number,idum);
// 		C_post      = C_post_mean + sigma_post*(NR::gasdev(idum));
//         //C_post_mean*(NR::gamdev(order,idum));// + sigma_post*(NR::gasdev(idum));
// 		if (C_post < 0.)
// 	  		C_post = 0.;
// 		//cout << C_post << endl;
// 		draw_post=false;
//       }
      yt[2] += C_post;
      post_spike=true;
      //cout << "postsynaptic spike at : " << t <<  endl;
    }
    else 
      post_s=false;

  }
  // variable for gaussian noise
  if (sigma > 0. ) {
  	eta = (NR::gasdev(idum))/(sqrt(ddt)); 
  }
  // evaluate differential equations
  functions(t,yt,dydt);
  rkck(yt,dydt,t,nvar,dt,yout); 
  
}

//-----------------------------------------------------------------------------

void pre_post_spikes::update(double t) {
  

  //time_phos = yt[3]/(double(n_presentation));
  //time_dephos = yt[4]/(double(n_presentation)); 

  if(!protocol)
    Ca = ca_rest;
  else if (protocol)
    Ca = ca_rest + yt[1] + yt[2]; 

   //if (t >= t_write ) {
	//cout << t << "  " << protocol << "  " << Ca << "  " << yt[1] << "  " << yt[2] << "  " << ddt << "  " << yt[0]  << "  " << C_pre << endl;
 	//t_write += dt_write;
  //}
  phossum=yt[0];
	
	if (buff_max < yt[3]) {
		buff_max = yt[3];
	}
	if (ca_max < Ca) {
		ca_max = Ca;
	}
   // if presyn. spike: increase the occuring time for next pre.-spike
  if (pre_spike) {
    t_pre_spike+=inter_pair;
    pre_count+=1;
    draw_pre = true;
    //if ( C_pre ) 
    //   success ++;
    //cout << t << "  " << Ca << "  " << yt[3] << "  " << vesicle_docked[0] << "  " << vesicle_docked[1] << "  " << n_release << "  " << yt[3] << endl;
  }
  // if postsyn. spike: increase the occuring time for next post.-spike
  if (post_spike) {
    t_post_spike+=inter_pair;
    post_count+=1;
    draw_post = true;
  }
  
  // checking for convergence either in the initialzing phase or after the stimulation protocol: "protocol=false"
 /* if (!protocol) {
    if(t > t_initialization && stage == 0) {
      stop = true;
      stage += 1;
      //motif_camkII_f << endl << endl << endl;
      
	//cout << "End of initialization: system has reached pre-stimulus steady-state for delta t = ";
	//cout << delta_t;
	//if (state)
	  //cout << " (down)." << endl;
	//else
	  //cout << " (up)." << endl;
	t_compare = t_stim_end;
	nvar = nvar_total;
   } */
   if ( t > t_simulation ) {
	stop = true;
    stage += 1;
	//cout << "End of simulation: system has reached final steady-state for delta t = ";
	//cout << delta_t;
	//if (state)
	  //cout << " (down)." << endl;
	//else
	  //cout << " (up)." << endl;
    }
    else 
		stop = false;
  
  
}
//-----------------------------------------------------------------------------

void pre_post_spikes::free() {
  
  delete  [ ] yt;
  delete  [ ] dydt;
  delete  [ ] yout;
  delete  [ ] dy;  
  delete  [ ] yerr;
  delete  [ ] yscal;
}

//-----------------------------------------------------------------------------

void pre_post_spikes::functions(double t, double yt[], double dydt[]) {
  
  
  if(!protocol)
    Ca = ca_rest;
  else if (protocol)
    Ca = ca_rest + yt[1] + yt[2]; 

  if ( Ca >= Ct_phos ) 
    alpha = phos;
  else 
    alpha = 0.;

  if ( Ca >= Ct_dephos )
    beta  = dephos; 
  else 
    beta  = 0.;
  
    // current
    dydt[0] = (alpha*(1.- yt[0]) - beta*yt[0] + sigma*sqrt(beta/dephos)*sqrt(tau_rho)*eta - epsilon*yt[0]*(1.-yt[0])*(0.5-yt[0]))/tau_rho;
     //alphad 
    //dydt[0] = (alpha*(1.- yt[0]) - beta*yt[0] + sigma*(beta/dephos)*sqrt(tau_rho)*eta - epsilon*yt[0]*(1.-yt[0])*(0.5-yt[0]))/tau_rho;
    //alphad + alphap
    //dydt[0] = (alpha*(1.- yt[0]) - beta*yt[0] + sigma*(beta/dephos+alpha/phos)*sqrt(tau_rho)*eta - epsilon*yt[0]*(1.-yt[0])*(0.5-yt[0]))/tau_rho;

   //dydt[0] = (phos*(1.- yt[0])/10. - dephos*yt[0]/10. + sigma*sqrt(tau_rho)*eta - epsilon*yt[0]*(1.-yt[0])*(0.5-yt[0]))/tau_rho;

   
	//cout << epsilon << endl;
   // during stimulation protocol evaluate calcium dynamics  =================================
  dydt[1] = -yt[1]/tau_pre; 
  dydt[2] = -yt[2]/tau_post;
  dydt[3] = (ca_rest + yt[1]+yt[2]);
  //dydt[3] = alpha/phos;
  //dydt[4] = beta/dephos;
    //dydt[3] = (p_release_0 - yt[3])/tau_facilitation_recovery;
  
}

//---------------------------------------------------------------------------------------------------------------------------

void pre_post_spikes::rkqs(double &t, const double htry, const double eps, double &hdid, double &hnext) {
  
  const double SAFETY=0.9, PGROW=-0.2, PSHRNK=-0.25, ERRCON=1.89e-4, TINY=1.0e-5;
  int i;
  double errmax,h,htemp,tnew;

  pre_spike = pre_s  = false;
  post_spike= post_s = false;
  //docking = false;
  if(!protocol)
    nvar = 1;
  h=htry;
//   for (;;) {
//     if (pre_spike) {
//       n_release = 0;
//       yt[1] -= C_pre; 
//       //yt[3] = p_release_old;
//     }
//     if (post_spike) {
// 	  yt[2] -= C_post;
//     }
//     pre_spike = pre_s  = false;
//     post_spike= post_s = false;
//     //docking = false;
    evaluate_event(t,h);
    // scaling used to monitor accuracy 
//     for (i=0; i<nvar; i++)
//       yscal[i] = fabs(yt[i]) + fabs(dydt[i]*h)+TINY;
//     errmax=0.0;
//     for (i=0;i<nvar;i++) {
//       errmax=MAX(errmax,fabs(yerr[i]/yscal[i]));
//       //errmax=MAX(errmax,fabs(yerr[i]));
//     }
//     errmax /= eps;
//     if (errmax <= 1.0) break;
//     htemp=SAFETY*h*pow(errmax,PSHRNK);
//     h=(h >= 0.0 ? MAX(htemp,0.1*h) : MIN(htemp,0.1*h));
    tnew=t+h;
//     if (tnew == t) nrerror("stepsize underflow in rkqs");
//   }
//   if (errmax > ERRCON) hnext=SAFETY*h*pow(errmax,PGROW);
//   else hnext=5.0*h;
  t += (hdid=h);
  dt_done = h;
  for (i=0;i<nvar;i++) yt[i]=yout[i];
  pre_s = post_s = false;
  update(t);
}

//-----------------------------------------------------------------------------

void pre_post_spikes::rkck(const double y[], const double dydx[], const double x, const int nva,
		    const double h, double yout[]) 
{
  static const double a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875,
    b21=0.2, b31=3.0/40.0, b32=9.0/40.0, b41=0.3, b42 = -0.9,
    b43=1.2, b51 = -11.0/54.0, b52=2.5, b53 = -70.0/27.0,
    b54=35.0/27.0, b61=1631.0/55296.0, b62=175.0/512.0,
    b63=575.0/13824.0, b64=44275.0/110592.0, b65=253.0/4096.0,
    c1=37.0/378.0, c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0,
    dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0, dc5 = -277.00/14336.0, dc6=c6-0.25;
  int i;
  
  //int n=y.size();
  double ak2[nva],ak3[nva],ak4[nva],ak5[nva],ak6[nva],ytemp[nva];
  
  for (i=0;i<nva;i++)
    ytemp[i]=y[i]+b21*h*dydx[i];
  functions(x+a2*h,ytemp,ak2);
  for (i=0;i<nva;i++)
    ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
  functions(x+a3*h,ytemp,ak3);
  for (i=0;i<nva;i++)
    ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
  functions(x+a4*h,ytemp,ak4);
  for (i=0;i<nva;i++)
    ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  functions(x+a5*h,ytemp,ak5);
  for (i=0;i<nva;i++)
    ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
  functions(x+a6*h,ytemp,ak6);
  for (i=0;i<nva;i++)
    yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
  for (i=0;i<nva;i++) {
    yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
  }
}



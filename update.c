
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "lattice.h"
#include "macros.h"
#include "ranlxd.h"
#include "local_meas.h"
#include "glob_meas.h"
#include "action.h"
#include "update.h"


double dS_dphi( site *s, int i ){
  int mu,j;
  double dSdp;
  site *sp,*sm;

  /*  dS/dphi_i = dS/d(d_mu phi_j) d(d_mu phi_j)/dphi_i  */
  dSdp = 0.0;
#ifdef DUAL_LATTICE
  for(mu=0;mu<4;mu++){
    for(j=0;j<16;j++){
      sm = (site *)((*s).dneib[j]);
      dSdp += dsign[mu][j] * 0.125 * dS_ddmuphi(sm,i,mu);
    }
  }
#else
  for(mu=0;mu<4;mu++){
    sp = (site *)((*s).neib[mu]);
    sm = (site *)((*s).neib[OPP_DIR(mu)]);
    dSdp += 0.5 * ( dS_ddmuphi(sm,i,mu) - dS_ddmuphi(sp,i,mu) );
  }
#endif

  return(dSdp);
}


double dS_dpsi( site *s ){
  int mu,j;
  double dSdp;
  site *sp,*sm;

  /*  dS/dpsi = dS/d(d_mu psi) d(d_mu psi)/dpsi  */
  dSdp = 0.0;
#ifdef DUAL_LATTICE
  for(mu=0;mu<4;mu++){
    for(j=0;j<16;j++){
      sm = (site *)((*s).dneib[j]);
      dSdp += dsign[mu][j] * 0.125 * dS_ddmupsi(sm,mu);
    }
  }
#else
  for(mu=0;mu<4;mu++){
    sp = (site *)((*s).neib[mu]);
    sm = (site *)((*s).neib[OPP_DIR(mu)]);
    dSdp += 0.5 * ( dS_ddmupsi(sm,mu) - dS_ddmupsi(sp,mu) );
  }
#endif

  return(dSdp);
}


int Met_update( ){
  register int n,is,i;
  register site *s;
  double dSdp,dp,dS,xx;
  int met;

  met = 0;

  for( n=0 ; n<params.Vol4 ; n++ ){
    /*  pick a random site  */
    //    is = (int)( params.Vol4 * drand48() );
    ranlxd( &xx , 1 );
    is = (int)( params.Vol4 * xx );
    s = &(lattice[is]);
    /*  update psi field  */
    dSdp = dS_dpsi(s);
    //    dp = -dSdp * params.dt * 2.0*(1.0-drand48());
    ranlxd( &xx , 1 );
    dp = -dSdp * params.dt * 2.0*(1.0-xx);
    dS = dSdp * dp;
    if( dS > 0.0 ){
      //      xx = drand48();
      ranlxd( &xx , 1 );
      if( xx < exp(-dS) ){
	(*s).psi += dp;
	met += 1;
      }
    }
    else{
      (*s).psi += dp;
      met += 1;
    }
    /*  update phi_i fields  */
    if( params.upd_phi ) 
      for(i=0;i<3;i++){
	dSdp = dS_dphi(s,i);
	//	dp = -dSdp * params.dt * 2.0*(1.0-drand48());
	ranlxd( &xx , 1 );
	dp = -dSdp * params.dt * 2.0*(1.0-xx);
	dS = dSdp * dp;
	if( dS > 0.0 ){
	  //	  xx = drand48();
	  ranlxd( &xx , 1 );
	  if( xx < exp(-dS) ){
	    (*s).phi[i] += dp;
	    met += 1;
	  }
	}
	else{
	  (*s).phi[i] += dp;
	  met += 1;
	}
      }
  }

  /*  calculate new field derivs  */
  if( params.bc != PERIODIC ) init_bound();
  field_derivs();
  Saction = action();

  return(met);
}


/*  gaussian random conjugate momenta (for MD updates)  */
void random_conj_mom( ){
  register int is,i;
  register site *s;
  double xx[4];

  if( params.upd_phi ){
    /*  #pragma omp parallel for private(is,s,xx)  */
    FORALLSITES(is,s) 
      //      for(i=0;i<4;i++) xx[i] = drand48();

      /*  #pragma omp critical  */
      /*  #pragma omp ordered  */
      ranlxd( xx , 4 );

      (*s).psiCM = sqrt(-2.0*log(xx[0]))*cos(2.0*M_PI*xx[1]);
      (*s).phiCM[0] = sqrt(-2.0*log(xx[0]))*sin(2.0*M_PI*xx[1]);
      (*s).phiCM[1] = sqrt(-2.0*log(xx[2]))*cos(2.0*M_PI*xx[3]);
      (*s).phiCM[2] = sqrt(-2.0*log(xx[2]))*sin(2.0*M_PI*xx[3]);
    ENDALLSITES 
  }
  else{
    /*  #pragma omp parallel for private(is,s,xx)  */
    FORALLSITES(is,s) 
      //      for(i=0;i<2;i++) xx[i] = drand48();

      /*  #pragma omp critical  */
      /*  #pragma omp ordered  */
      ranlxd( xx , 2 );

      (*s).psiCM = sqrt(-2.0*log(xx[0]))*cos(2.0*M_PI*xx[1]);
      (*s).phiCM[0] = (*s).phiCM[1] = (*s).phiCM[2] = 0.0;
    ENDALLSITES 
  }

  return;
}


int HMC_update( ){
  register int is,i;
  register site *s;
  double tau,S_old,H_old,xx;
  int met;

  met = 0;
  S_old = Saction;
#pragma omp parallel for private(is,s,i)
  FORALLSITES(is,s) 
    (*s).psi_old = (*s).psi;
    for(i=0;i<3;i++) (*s).phi_old[i] = (*s).phi[i];
  ENDALLSITES 

  /*  Molecular Dynamics (MD) Updates:  */

  /*  initially random conjugate momenta  */
  random_conj_mom();
  H_old = hamiltonian( S_old );

  /*  initial half-step for conjugate momenta  */
#pragma omp parallel for private(is,s,i)
  FORALLSITES(is,s) 
    (*s).psiCM -= 0.5 * params.dt * dS_dpsi(s);
    if( params.upd_phi ) 
      for(i=0;i<3;i++) (*s).phiCM[i] -= 0.5 * params.dt * dS_dphi(s,i);
  ENDALLSITES 

  /*  MD trajectory (leap frog)  */
  for( tau=0.0 ; tau<(1.0-0.5*params.dt) ; tau+=params.dt ){

    /*  update fields  */
#pragma omp parallel for private(is,s,i)
    FORALLSITES(is,s) 
      (*s).psi += params.dt * (*s).psiCM;
      if( params.upd_phi ) 
	for(i=0;i<3;i++) (*s).phi[i] += params.dt * (*s).phiCM[i];
    ENDALLSITES 

    /*  calculate new field derivs  */
    if( params.bc != PERIODIC ) init_bound();
    field_derivs();

    /*  update conjugate momenta  */
    if( tau < (1.0-1.5*params.dt) ){
#pragma omp parallel for private(is,s,i)
      FORALLSITES(is,s) 
	(*s).psiCM -= params.dt * dS_dpsi(s);
	if( params.upd_phi ) 
	  for(i=0;i<3;i++) (*s).phiCM[i] -= params.dt * dS_dphi(s,i);
      ENDALLSITES 
    }

  }

  /*  final half-step for conjugate momenta  */
#pragma omp parallel for private(is,s,i)
  FORALLSITES(is,s) 
    (*s).psiCM -= 0.5 * params.dt * dS_dpsi(s);
    if( params.upd_phi ) 
      for(i=0;i<3;i++) (*s).phiCM[i] -= 0.5 * params.dt * dS_dphi(s,i);
  ENDALLSITES 

  Saction = action();
  Hamilto = hamiltonian( Saction );

  /*  Metropolis accept/reject:  */
  if( params.Met_step ){
    if( Hamilto > H_old ){
      //      xx = drand48();
      ranlxd( &xx , 1 );
      if( xx > exp(H_old-Hamilto) ){
#pragma omp parallel for private(is,s,i)
	FORALLSITES(is,s) 
	  (*s).psi = (*s).psi_old;
	  for(i=0;i<3;i++) (*s).phi[i] = (*s).phi_old[i];
	ENDALLSITES 
	if( params.bc != PERIODIC ) init_bound();
	field_derivs();
	Saction = S_old;
	Hamilto = H_old;
      }
      else met = 1;
    }
    else met = 1;
  }

  return(met);
}



/*  "Ideal-gas" action: dS = S_1 d^4x b^4/3  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lattice.h"
#include "macros.h"
#include "action.h"

#define THIRD 0.33333333333333
#define TWOTHIRD 0.66666666666667
#define FOURTHIRD 1.3333333333333
#define FOURNINTH 0.44444444444444


void action_init( ){
  params.b0 = pow( fabs(params.S1) , -0.75 );
  params.y0 = 0.0;
  params.bc = PERIODIC;
  return;
}


/*  dL(x) = F(b,X,y ; x) [= p - T*dp/dT ; from EoS (thermo/hydro matching)]  */
double dLagrangian( site *s ){
  double dL;
  dL = params.S1 * pow( (*s).b , FOURTHIRD );
  return(dL);
}


/*  dF/db , dF/dX , dF/dy  */
void dF_dbXy( double *dF, site *s ){
  dF[0] = params.S1 * FOURTHIRD * pow( (*s).b , THIRD );
  dF[1] = 0.0;
  dF[2] = 0.0;
  return;
}


/*  d^2F/db^2 , d^2F/dbdX , d^2F/dbdy ,d^2F/dX^2 , d^2F/dXdy , d^2F/dy^2  */
void d2F_dbXy2( double *d2F, site *s ){
  d2F[0] = params.S1 * FOURNINTH * pow( (*s).b , -TWOTHIRD );
  d2F[1] = 0.0;
  d2F[2] = 0.0;
  d2F[3] = 0.0;
  d2F[4] = 0.0;
  d2F[5] = 0.0;
  return;
}


void EoS_match( double *n_rho_p_T , site *s ){
  double F,dF[3];
  F = dLagrangian( s );
  dF_dbXy( dF , s );
  n_rho_p_T[0] = 0.0;
  //  n_rho_p_T[0] = -dF[2] + 2.0 * dF[1] * (*s).y;
  n_rho_p_T[1] = -(*s).y * dF[2] + F + 2.0 * (*s).y * (*s).y * dF[1];
  n_rho_p_T[2] = -F + (*s).b * dF[0];
  n_rho_p_T[3] = dF[0];
  return;
}


void Tmn_match( Ttensor *Tmn , site *s ){
  int mu,nu,munu;
  double n_rho_p_T[4],dF[3],F;
  EoS_match( n_rho_p_T , s );
  //  F = dLagrangian( s );
  //  dF_dbXy( dF , s );
  /*  T_mu,nu = (p + rho)*u_mu*u_nu - p*g_mu,nu  */
  for( mu=0,munu=0 ; mu<4 ; mu++ ) 
    for( nu=mu ; nu<4 ; nu++,munu++ ){
      (*Tmn).dir12[munu] = ( n_rho_p_T[1] + n_rho_p_T[2] ) * 
	(*s).uu[mu] * (*s).uu[nu];
      if ( nu == mu ) (*Tmn).dir12[munu] -= n_rho_p_T[2];
      /*
      (*Tmn).dir12[munu] = ( (*s).y * dF[2] - (*s).b * dF[0] ) * 
	(*s).uu[mu] * (*s).uu[nu];
      (*Tmn).dir12[munu] -= 2.0 * dF[1] * (*s).dpsi[mu] * (*s).dpsi[nu];
      if( nu == mu ) (*Tmn).dir12[munu] += F - (*s).b * dF[0];
      */
    }
  return;
}


double dS_ddmuphi( site *s, int i, int mu ){
  int j;
  double dSdd;

  /*  dS/d(d_mu phi_i) ~ 4/3 b^4/3 (B^-1)_ij (d_mu phi_j)  */
  dSdd = 0.0;
  for(j=0;j<3;j++) dSdd += (*s).Binv[i][j] * (*s).dphi[j][mu];
  dSdd *= params.S1 * FOURTHIRD * pow( (*s).b , FOURTHIRD );

  return(dSdd);
}


double dS_ddmupsi( site *s, int mu ){
  double dSdd;

  /*  dS/d(d_mu psi) = 0  */
  dSdd = 0.0;

  return(dSdd);
}


void init_bound( ){
  return;
}



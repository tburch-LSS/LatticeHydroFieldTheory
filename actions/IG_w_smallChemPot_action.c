/*  "Ideal-gas" action with small chemical potential: 
    dS = S_1 d^4x ( b^4/3 + S_2 b y + S_3 b^2/3 y^2)  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lattice.h"
#include "macros.h"
#include "action.h"

#define THIRD 0.33333333333333
#define TWOTHIRD 0.66666666666667
#define FOURTHIRD 1.3333333333333
#define TWONINTH 0.22222222222222
#define FOURNINTH 0.44444444444444


void action_init( ){
  params.b0 = pow( fabs(params.S1) , -0.75 );
  params.y0 = 0.0;
  params.bc = PERIODIC;
  return;
}


/*  dL(x) = F(b,X,y ; x) [= p - T*dp/dT ; from EoS (thermo/hydro matching)]  */
double dLagrangian( site *s ){
  double b23,dL;
  b23 = pow( (*s).b , TWOTHIRD );
  dL = params.S1 * ( b23 * b23 
		     + params.S2 * (*s).b * (*s).y 
		     + params.S3 * b23 * (*s).y * (*s).y );
  return(dL);
}


/*  dF/db , dF/dX , dF/dy  */
void dF_dbXy( double *dF, site *s ){
  double b13;
  b13 = pow( (*s).b , THIRD );
  dF[0] = params.S1 * ( FOURTHIRD * b13 + params.S2 * (*s).y + 
			params.S3 * TWOTHIRD * (*s).y * (*s).y / b13 );
  dF[1] = 0.0;
  dF[2] = params.S1 * ( params.S2 * (*s).b + 
			2.0 * params.S3 * b13 * b13 * (*s).y );
  return;
}


/*  d^2F/db^2 , d^2F/dbdX , d^2F/dbdy ,d^2F/dX^2 , d^2F/dXdy , d^2F/dy^2  */
void d2F_dbXy2( double *d2F, site *s ){
  double b13,b23;
  b13 = pow( (*s).b , THIRD ); b23 = b13 * b13;
  d2F[0] = params.S1 * ( FOURNINTH / b23 - params.S3 * TWONINTH / ( b23 * b23 ) );
  d2F[1] = 0.0;
  d2F[2] = params.S1 * ( params.S2 + 
			 params.S3 * FOURTHIRD * (*s).y / b13 );
  d2F[3] = 0.0;
  d2F[4] = 0.0;
  d2F[5] = params.S1 * ( 2.0 * params.S3 * b13 * b13 );
  return;
}


void EoS_match( double *n_rho_p_T , site *s ){
  double F,dF[3];
  F = dLagrangian( s );
  dF_dbXy( dF , s );
  //  n_rho_p_T[0] = 0.0;
  n_rho_p_T[0] = -dF[2] + 2.0 * dF[1] * (*s).y;
  n_rho_p_T[1] = -(*s).y * dF[2] + F + 2.0 * (*s).y * (*s).y * dF[1];
  n_rho_p_T[2] = -F + (*s).b * dF[0];
  n_rho_p_T[3] = dF[0];
  return;
}


void Tmn_match( Ttensor *Tmn , site *s ){
  int mu,nu,munu;
  double n_rho_p_T[4],dF[3],F;
  //  EoS_match( n_rho_p_T , s );
  F = dLagrangian( s );
  dF_dbXy( dF , s );
  for( mu=0,munu=0 ; mu<4 ; mu++ ) 
    for( nu=mu ; nu<4 ; nu++,munu++ ){
      /*
      (*Tmn).dir12[munu] = ( n_rho_p_T[1] + n_rho_p_T[2] ) * 
	(*s).uu[mu] * (*s).uu[nu];
      if ( nu == mu ) (*Tmn).dir12[munu] -= n_rho_p_T[2];
      */
      (*Tmn).dir12[munu] = ( (*s).y * dF[2] - (*s).b * dF[0] ) * 
	(*s).uu[mu] * (*s).uu[nu];
      (*Tmn).dir12[munu] -= 2.0 * dF[1] * (*s).dpsi[mu] * (*s).dpsi[nu];
      if( nu == mu ) (*Tmn).dir12[munu] += F - (*s).b * dF[0];
    }
  return;
}


double dS_ddmuphi( site *s, int i, int mu ){
  int j;
  double b23,dSdd;

  /*  dS/d(d_mu phi_i) ~ 4/3 b^4/3 (B^-1)_ij (d_mu phi_j) + ...  */
  b23 = pow( (*s).b , TWOTHIRD );
  dSdd = params.S1 * ( FOURTHIRD * b23 * b23 
		       //		       + params.S2 * (*s).b * (*s).y 
		       - params.S3 * FOURTHIRD * b23 * (*s).y * (*s).y ) * (*s).AImu[i][mu];

  return(dSdd);
}


double dS_ddmupsi( site *s, int mu ){
  double dSdd;

  /*  dS/d(d_mu psi) =   */
  dSdd = params.S1 * ( params.S2 * (*s).b + 
		       params.S3 * 2.0 * pow( (*s).b , TWOTHIRD ) * (*s).y ) * (*s).uu[mu];

  return(dSdd);
}


void init_bound( ){
  return;
}



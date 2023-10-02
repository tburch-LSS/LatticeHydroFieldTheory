/*  "alternative" relativistic, degenerate fermion action: 
    p = C_1 mu^4 + A_2 mu^2 T^2 
    p - T dp/dT = C_1 mu^4 - A_2 mu^2 T^2 
    p - T dp/dT = C_1 mu^4 - C_2 (s/mu)^2 
    dS = d^4x [S_1 X^2 - S_2 (b/y)^2]  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lattice.h"
#include "macros.h"
#include "action.h"


void action_init( ){
  if( params.S2 > 0.0 ) params.b0 = 1.0 / sqrt( params.S1 * params.S2 );
  else params.b0 = 1.0;
  //  params.y0 = pow( params.S1 , -0.25 );
  params.y0 = 1.0;
  params.bc = PERIODIC;
  //  params.bc = TIME_NeuBC;
  //  params.bc = TIME_MixBC;
  //  params.bc = ALL_NeuBC;
  return;
}


/*  dL(x) = F(b,X,y ; x) [= p - T dp/dT ; from EoS (thermo/hydro matching)]  */
double dLagrangian( site *s ){
  double dL;
  dL = params.S1 * (*s).X * (*s).X - params.S2 * (*s).b * (*s).b / ( (*s).y * (*s).y );
  return(dL);
}


/*  dF/db , dF/dX , dF/dy  */
void dF_dbXy( double *dF, site *s ){
  double tmp;
  tmp = (*s).y * (*s).y;
  dF[0] = -params.S2 * 2.0 * (*s).b / tmp;
  dF[1] = params.S1 * 2.0 * (*s).X;
  dF[2] = params.S2 * 2.0 * (*s).b * (*s).b / ( tmp * (*s).y );
  return;
}


/*  d^2F/db^2 , d^2F/dbdX , d^2F/dbdy ,d^2F/dX^2 , d^2F/dXdy , d^2F/dy^2  */
void d2F_dbXy2( double *d2F, site *s ){
  double tmp;
  tmp = (*s).y * (*s).y;
  d2F[0] = -params.S2 * 2.0 / tmp;
  d2F[1] = 0.0;
  d2F[2] = params.S2 * 2.0 * (*s).b / ( tmp * (*s).y );
  d2F[3] = params.S1 * 2.0;
  d2F[4] = 0.0;
  d2F[5] = -params.S2 * 6.0 * (*s).b * (*s).b / ( tmp * tmp );
  return;
}


void EoS_match( double *n_rho_p_T , site *s ){
  double F,dF[3];
  F = dLagrangian( s );
  dF_dbXy( dF , s );
  n_rho_p_T[0] = 2.0 * dF[1] * (*s).y;
  n_rho_p_T[1] = 3.0 * F;
  n_rho_p_T[2] = F;
  n_rho_p_T[3] = 0.0;
  /*
  n_rho_p_T[0] = dF[2] - 2.0 * dF[1] * (*s).y;
  n_rho_p_T[1] = (*s).y * dF[2] - F - 2.0 * (*s).y * (*s).y * dF[1];
  n_rho_p_T[2] = F - (*s).b * dF[0];
  n_rho_p_T[3] = -dF[0];
  */
  return;
}


void Tmn_match( Ttensor *Tmn , site *s ){
  int mu,nu,munu;
  double n_rho_p_T[4],dF[3],F,vv[4];
  //  EoS_match( n_rho_p_T , s );
  F = dLagrangian( s );
  dF_dbXy( dF , s );
  /*  T_mu,nu = (p + rho)*u_mu*u_nu - p*g_mu,nu  */
  for( mu=0 ; mu<4 ; mu++ ) vv[mu] = (*s).dpsi[mu] / sqrt( (*s).X );
  for( mu=0,munu=0 ; mu<4 ; mu++ ) 
    for( nu=mu ; nu<4 ; nu++,munu++ ){
      (*Tmn).dir12[munu] = ( 4.0 * F ) * vv[mu] * vv[nu];
      if ( nu == mu ) (*Tmn).dir12[munu] -= F;
      /*
      (*Tmn).dir12[munu] = ( n_rho_p_T[1] + n_rho_p_T[2] ) * 
        (*s).uu[mu] * (*s).uu[nu];
      if ( nu == mu ) (*Tmn).dir12[munu] -= n_rho_p_T[2];
      */
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

  /*  dS/d(d_mu phi_i) ~ -2 (b/y)^2 (B^-1)_ij (d_mu phi_j)  */
  dSdd = 0.0;
  for(j=0;j<3;j++) dSdd += (*s).Binv[i][j] * (*s).dphi[j][mu];
  dSdd *= -2.0 * params.S2 * (*s).b * (*s).b / ( (*s).y * (*s).y );

  return(dSdd);
}


double dS_ddmupsi( site *s, int mu ){
  double dSdd;

  /*  dS/d(d_mu psi) ~ 4 (d_mu psi) X + 2C u_mu b^2 y^-3  */
  dSdd = 4.0 * (*s).dpsi[mu] * params.S1 * (*s).X + 
    2.0 * (*s).uu[mu] * params.S2 * (*s).b * (*s).b / ( (*s).y * (*s).y * (*s).y );

  return(dSdd);
}


void init_bound( ){
  register int is,t,i,j,k,js,x,l;
  register site *s,*sn;
  register double cp;
  int dr[4];
  /*
  cp = 0.0;
  FORALLSITES(is,s){
    cp += (*s).y;
  }
  cp /= (double) Nsite;
  if( cp <= 0.0 ) cp = 1.0;
  */
  if( params.bc == TIME_MixBC ){
    for(i=0;i<3;i++) dr[i] = 0;
    dr[3] = -1;
    k = 0;
    FORALLSPACE(is,s,0){
      sn = (site *) (*s).neib[TUP];
      latt_bound[k].psi = (*sn).psi - 2.0 * params.chem_pot;
      //      latt_bound[k].psi = (*sn).psi - 2.0 * cp;
      js = new_index_from_dr(is,dr);
      sn = &(lattice[js]);
      for(i=0;i<3;i++) 
	latt_bound[k].phi[i] = (*sn).phi[i];
      for(i=0;i<4;i++) latt_bound[k].dpsi[i] = 0.0;
      k++;
    }
    dr[3] = 1;
    t = params.Nt - 1;
    FORALLSPACE(is,s,t){
      sn = (site *) (*s).neib[TDOWN];
      latt_bound[k].psi = (*sn).psi + 2.0 * params.chem_pot;
      //      latt_bound[k].psi = (*sn).psi + 2.0 * cp;
      js = new_index_from_dr(is,dr);
      sn = &(lattice[js]);
      for(i=0;i<3;i++) 
	latt_bound[k].phi[i] = (*sn).phi[i];
      for(i=0;i<4;i++) latt_bound[k].dpsi[i] = 0.0;
      k++;
    }
  }
  else if( params.bc == TIME_NeuBC || params.bc == ALL_NeuBC ){
    k = 0;
    FORALLSPACE(is,s,0){
      sn = (site *) (*s).neib[TUP];
      latt_bound[k].psi = (*sn).psi - 2.0 * params.chem_pot;
      //      latt_bound[k].psi = (*sn).psi - 2.0 * cp;
      for(i=0;i<3;i++){
	latt_bound[k].phi[i] = 2.0 * (*s).phi[i] - (*sn).phi[i];
	for(j=0;j<3;j++) latt_bound[k].Binv[i][j] = 0.0;
      }
      for(i=0;i<4;i++) latt_bound[k].dpsi[i] = 0.0;
      latt_bound[k].X = 1.0;
      k++;
    }
    t = params.Nt - 1;
    FORALLSPACE(is,s,t){
      sn = (site *) (*s).neib[TDOWN];
      latt_bound[k].psi = (*sn).psi + 2.0 * params.chem_pot;
      //      latt_bound[k].psi = (*sn).psi + 2.0 * cp;
      for(i=0;i<3;i++){
	latt_bound[k].phi[i] = 2.0 * (*s).phi[i] - (*sn).phi[i];
	for(j=0;j<3;j++) latt_bound[k].Binv[i][j] = 0.0;
      }
      for(i=0;i<4;i++) latt_bound[k].dpsi[i] = 0.0;
      latt_bound[k].X = 1.0;
      k++;
    }
  }
  if( params.bc == ALL_NeuBC ){
    for(i=0;i<4;i++) dr[i] = 0;
    dr[0] = -1;
    FORALLXBOUND(is,s,0){
      sn = (site *) (*s).neib[XUP];
      for(i=0;i<3;i++){
	latt_bound[k].phi[i] = 2.0 * (*s).phi[i] - (*sn).phi[i];
	for(j=0;j<3;j++) latt_bound[k].Binv[i][j] = 0.0;
      }
      js = new_index_from_dr(is,dr);
      sn = &(lattice[js]);
      latt_bound[k].psi = (*sn).psi;
      k++;
    }
    dr[0] = 1;
    x = params.Nt - 1;
    FORALLXBOUND(is,s,x){
      sn = (site *) (*s).neib[XDOWN];
      for(i=0;i<3;i++){
	latt_bound[k].phi[i] = 2.0 * (*s).phi[i] - (*sn).phi[i];
	for(j=0;j<3;j++) latt_bound[k].Binv[i][j] = 0.0;
      }
      js = new_index_from_dr(is,dr);
      sn = &(lattice[js]);
      latt_bound[k].psi = (*sn).psi;
      k++;
    }
    for(i=0;i<4;i++) dr[i] = 0;
    dr[1] = -1;
    FORALLYBOUND(is,s,l,0){
      sn = (site *) (*s).neib[YUP];
      for(i=0;i<3;i++){
	latt_bound[k].phi[i] = 2.0 * (*s).phi[i] - (*sn).phi[i];
	for(j=0;j<3;j++) latt_bound[k].Binv[i][j] = 0.0;
      }
      js = new_index_from_dr(is,dr);
      sn = &(lattice[js]);
      latt_bound[k].psi = (*sn).psi;
      k++;
    }
    dr[1] = 1;
    x = params.Nt - 1;
    FORALLYBOUND(is,s,l,x){
      sn = (site *) (*s).neib[YDOWN];
      for(i=0;i<3;i++){
	latt_bound[k].phi[i] = 2.0 * (*s).phi[i] - (*sn).phi[i];
	for(j=0;j<3;j++) latt_bound[k].Binv[i][j] = 0.0;
      }
      js = new_index_from_dr(is,dr);
      sn = &(lattice[js]);
      latt_bound[k].psi = (*sn).psi;
      k++;
    }
    for(i=0;i<4;i++) dr[i] = 0;
    dr[2] = -1;
    FORALLZBOUND(is,s,l,0){
      sn = (site *) (*s).neib[ZUP];
      for(i=0;i<3;i++){
	latt_bound[k].phi[i] = 2.0 * (*s).phi[i] - (*sn).phi[i];
	for(j=0;j<3;j++) latt_bound[k].Binv[i][j] = 0.0;
      }
      js = new_index_from_dr(is,dr);
      sn = &(lattice[js]);
      latt_bound[k].psi = (*sn).psi;
      k++;
    }
    dr[2] = 1;
    x = params.Nt - 1;
    FORALLZBOUND(is,s,l,x){
      sn = (site *) (*s).neib[ZDOWN];
      for(i=0;i<3;i++){
	latt_bound[k].phi[i] = 2.0 * (*s).phi[i] - (*sn).phi[i];
	for(j=0;j<3;j++) latt_bound[k].Binv[i][j] = 0.0;
      }
      js = new_index_from_dr(is,dr);
      sn = &(lattice[js]);
      latt_bound[k].psi = (*sn).psi;
      k++;
    }
  }
  return;
}




#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "lattice.h"
#include "macros.h"
#include "action.h"
#include "local_meas.h"
#include "glob_meas.h"


/*  S = Int dt Lagr(b,X,y) = Int d^4x dLagr(b,X,y;x)  */
double action( ){
  register int is;
  register site *s;
  double S;

  S = 0.0;
#pragma omp parallel for private(is,s) reduction(+:S)
  FORALLSITES(is,s) 
    S += dLagrangian( s );
  ENDALLSITES 

  return(S);
}


/*  H = Sum_x,i pi_i(x)^2 + S  */
double hamiltonian( double S ){
  register int is;
  register site *s;
  double H;

  H = 0.0;
#pragma omp parallel for private(is,s) reduction(+:H)
  FORALLSITES(is,s) 
    H += 0.5 * ( (*s).psiCM*(*s).psiCM + (*s).phiCM[0]*(*s).phiCM[0] + 
		 (*s).phiCM[1]*(*s).phiCM[1] + (*s).phiCM[2]*(*s).phiCM[2] );
  ENDALLSITES 

  H += S;
  return(H);
}


/*  < dS/db > , < b dS/db > , < dS/dX > , < X dS/dX > , & 
    < dS/dy > , < y dS/dy > , etc.  */
void dS_dbdXdy( double *output ){
  register int is,n;
  register site *s;
  double dF[6],n_rho_p_T[4];
  double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8;

  tmp1 = tmp2 = tmp3 = tmp4 = tmp5 = tmp6 = tmp7 = tmp8 = 0.0;

#pragma omp parallel for private(is,s,dF) reduction(+:tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8)
  FORALLSITES(is,s) 
    /*  1st derivs: d/db , d/dX, d/dy  */
    dF_dbXy( dF , s );
    tmp1 += dF[0];
    tmp2 += (*s).b * dF[0];
    tmp3 += dF[1];
    tmp4 += (*s).X * dF[1];
    tmp5 += dF[2];
    tmp6 += (*s).y * dF[2];
    /*  oddballs: y d/dX , y^2 d/dX  */
    tmp7 += (*s).y * dF[1];
    tmp8 += (*s).y * (*s).y * dF[1];
  ENDALLSITES 

  output[0] = tmp1;
  output[1] = tmp2;
  output[2] = tmp3;
  output[3] = tmp4;
  output[4] = tmp5;
  output[5] = tmp6;
  output[6] = tmp7;
  output[7] = tmp8;

  tmp1 = tmp2 = tmp3 = tmp4 = tmp5 = tmp6 = tmp7 = tmp8 = 0.0;

#pragma omp parallel for private(is,s,dF) reduction(+:tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8)
  FORALLSITES(is,s) 
    /*  2nd derivs: d^2/db^2 , d^2/dbdX , d^2/dbdy , d^2/dX^2 , d^2/dXdy , 
	d^2/dy^2  */
    d2F_dbXy2( dF , s );
    tmp1 += dF[0];
    tmp2 += (*s).b * (*s).b * dF[0];
    tmp3 += dF[1];
    tmp4 += (*s).b * (*s).X * dF[1];
    tmp5 += dF[2];
    tmp6 += (*s).b * (*s).y * dF[2];
    tmp7 += dF[3];
    tmp8 += (*s).X * (*s).X * dF[3];
  ENDALLSITES 

  output[8] = tmp1;
  output[9] = tmp2;
  output[10] = tmp3;
  output[11] = tmp4;
  output[12] = tmp5;
  output[13] = tmp6;
  output[14] = tmp7;
  output[15] = tmp8;

  tmp1 = tmp2 = tmp3 = tmp4 = tmp5 = tmp6 = tmp7 = tmp8 = 0.0;

#pragma omp parallel for private(is,s,dF) reduction(+:tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7)
  FORALLSITES(is,s) 
    /*  2nd derivs: d^2/db^2 , d^2/dbdX , d^2/dbdy , d^2/dX^2 , d^2/dXdy , 
	d^2/dy^2  */
    d2F_dbXy2( dF , s );
    tmp1 += dF[4];
    tmp2 += (*s).X * (*s).y * dF[4];
    tmp3 += dF[5];
    tmp4 += (*s).y * (*s).y * dF[5];
    /*  oddballs: b y^2 d^2/dbdX , y^4 d^2/dX^2 , y^3 d^2/dXdy  */
    tmp5 += (*s).b * (*s).y * (*s).y * dF[1];
    tmp6 += (*s).y * (*s).y * (*s).y * (*s).y * dF[3];
    tmp7 += (*s).y * (*s).y * (*s).y * dF[4];
  ENDALLSITES 

  output[16] = tmp1;
  output[17] = tmp2;
  output[18] = tmp3;
  output[19] = tmp4;
  output[20] = tmp5;
  output[21] = tmp6;
  output[22] = tmp7;

  tmp1 = tmp2 = tmp3 = tmp4 = tmp5 = tmp6 = tmp7 = tmp8 = 0.0;

#pragma omp parallel for private(is,s,n_rho_p_T) reduction(+:tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8)
  FORALLSITES(is,s) 
    /*  total scalar invariants: <b> , <X> , <y> , <xi> , <n> ,<rho> , <p> , <T>  */
    tmp1 += (*s).b;
    tmp2 += (*s).X;
    tmp3 += (*s).y;
    tmp4 += sqrt( (*s).X - (*s).y * (*s).y );
    EoS_match( s , n_rho_p_T );
    tmp5 += n_rho_p_T[0];
    tmp6 += n_rho_p_T[1];
    tmp7 += n_rho_p_T[2];
    tmp8 += n_rho_p_T[3];
  ENDALLSITES 

  output[23] = tmp1;
  output[24] = tmp2;
  output[25] = tmp3;
  output[26] = tmp4;
  output[27] = tmp5;
  output[28] = tmp6;
  output[29] = tmp7;
  output[30] = tmp8;

  for(n=0;n<NUM_LOCAL_OBS;n++) output[n] /= (double) Nsite;

  return;
}


/*  flow tensor: Omega_mu,nu = u_mu u_nu + g_mu,nu 
    = (B^-1)_ij (d_mu phi^i) (d_nu phi^j)  */
void Correlators( Ttensor *Tmn , Ttensor_corr *TmnTmn , 
		  flow *Om , flow_corr *OmOm , 
		  phon_corr *PiPi , scal_corrs *scsc ){
  register int is,mu,nu,munu,munu2,i,j,k,js,t;
  register site *s,*s2;
  field_offset nvfo,svfo,fophi;
  int dr[4];
  double dTT,dO,dO2,dOO,dF[6],o1,o2,F,F2;
  Ttensor Tmn_src;
  flow Om_src;
  double Sclr_src[8],n_rho_p_T[8];
  double corr1,corr2,corr3,corr4,corr5,corr6,corr7,corr8;
  Ttensor *Tmn_x;
  flow *Om_x;

  Tmn_x = (Ttensor *) malloc( Nsite * sizeof(Ttensor) );
  Om_x = (flow *) malloc( Nsite * sizeof(flow) );

  /*  T_mu,nu(x)  &  < T_mu,nu >  ;  Omega_mu,nu(x)  &  < Omega_mu,nu >  ; 
      xi_mu(x)  &  < xi_mu > :  */

  for( munu=0 ; munu<10 ; munu++ ) (*Tmn).dir12[munu] = 0.0;
  for( munu=0 ; munu<10 ; munu++ ) (*Om).ndir12[munu] = 0.0;
  for( mu=0 ; mu<4 ; mu++ ) { (*Om).u1[mu] = 0.0; (*Om).J1[mu] = 0.0; 
    (*Om).uL1[mu] = 0.0; (*Om).JL1[mu] = 0.0; 
    (*Om).sdir[mu] = 0.0; (*Om).pi[mu] = 0.0; }
  for( munu=0 ; munu<6 ; munu++ ) (*Om).nvort[munu] = 0.0;
  for( munu=0 ; munu<6 ; munu++ ) (*Om).svort[munu] = 0.0;

  nvfo = F_OFFSET( uu );
  svfo = F_OFFSET( dpsi );

#pragma omp parallel for private(is,s,mu,nu,munu,i,j,dO,fophi)
  FORALLSITES(is,s) 
    Tmn_match( s , &(Tmn_x[is]) );

    eucl_to_lrnz_vel( s , Om_x[is].uL1 );
    //    lrnz_flow_match( s , &(Om_x[is]) );

    vorticity( nvfo , s , Om_x[is].nvort );
    //    vorticity( svfo , s , Om_x[is].svort );

    for( mu=0,munu=0 ; mu<4 ; mu++ ){
      Om_x[is].u1[mu] = (*s).uu[mu];
      Om_x[is].J1[mu] = (*s).b * (*s).uu[mu];
      Om_x[is].JL1[mu] = (*s).b * Om_x[is].uL1[mu];
      Om_x[is].sdir[mu] = 0.0;
      if( mu < 3 ){
	(*Om).pi[mu] += (*s).phi[mu];
	//	Om_x[is].pi[mu] = (*s).phi[mu];
	fophi = F_OFFSET(phi) + mu*sizeof(double);
	Om_x[is].pi[mu] = avgBC( fophi , s );
      }
      else{
	(*Om).pi[mu] += (*s).psi;
	//	Om_x[is].pi[mu] = (*s).psi;
	Om_x[is].pi[mu] = avgBC( F_OFFSET(psi) , s );
      }
      for( nu=mu ; nu<4 ; nu++,munu++ ){
	Om_x[is].ndir12[munu] = 0.0;
	//	for(i=0;i<3;i++) 
	//	  for(j=0;j<3;j++){
	//	    dO = (*s).Binv[i][j] * (*s).dphi[i][mu] * (*s).dphi[j][nu];
	//	    Om_x[is].ndir12[munu] += dO;
	//	  }
	if( params.upd_phi ){
	  for(i=0;i<3;i++){
	    dO = (*s).AImu[i][nu] * (*s).dphi[i][mu];
	    Om_x[is].ndir12[munu] += dO;
	  }
	}
	else {
	  Om_x[is].ndir12[munu] = -(*s).uu[mu] * (*s).uu[nu];
	  if( nu == mu ) Om_x[is].ndir12[munu] += 1.0;
	}
      }
    }
    for( mu=0,munu=0 ; mu<4 ; mu++ ) 
      for( nu=mu ; nu<4 ; nu++,munu++ ){
	dO = Om_x[is].ndir12[munu] * (*s).dpsi[nu];
	Om_x[is].sdir[mu] += dO;
	if( nu > mu ){
	  dO = Om_x[is].ndir12[munu] * (*s).dpsi[mu];
	  Om_x[is].sdir[nu] += dO;
	}
      }
  ENDALLSITES 

    for( mu=0,munu=0,munu2=0 ; mu<4 ; mu++ ){
      dO = dO2 = 0.0;
#pragma omp parallel for private(is,s) reduction(+:dO,dO2)
      FORALLSITES(is,s) 
	dO += Om_x[is].u1[mu];
	dO2 += Om_x[is].J1[mu];
      ENDALLSITES 
      (*Om).u1[mu] = dO;
      (*Om).J1[mu] = dO2;
      dO = dO2 = 0.0;
#pragma omp parallel for private(is,s) reduction(+:dO,dO2)
      FORALLSITES(is,s) 
	dO += Om_x[is].uL1[mu];
	dO2 += Om_x[is].JL1[mu];
      ENDALLSITES 
      (*Om).uL1[mu] = dO;
      (*Om).JL1[mu] = dO2;
      for( nu=mu ; nu<4 ; nu++,munu++ ){
	dTT = dOO = dO = dO2 = 0.0;
#pragma omp parallel for private(is,s) reduction(+:dTT,dOO,dO,dO2)
	FORALLSITES(is,s) 
	  dTT += Tmn_x[is].dir12[munu];
	  dOO += Om_x[is].ndir12[munu];
	  if( nu == mu ) dO += Om_x[is].sdir[mu];
	  else{
	    dO += Om_x[is].nvort[munu2];
	    dO2 += Om_x[is].svort[munu2];
	  }
	ENDALLSITES 
	(*Tmn).dir12[munu] = dTT;
	(*Om).ndir12[munu] = dOO;
	if( nu == mu ) (*Om).sdir[mu] = dO;
	else{
	  (*Om).nvort[munu2] = dO;
	  (*Om).svort[munu2] = dO2;
	  munu2++;
	}
      }
    }

  for( munu=0 ; munu<10 ; munu++ ){
    (*Tmn).dir12[munu] /= (double) Nsite;
    (*Om).ndir12[munu] /= (double) Nsite;
  }
  for( mu=0 ; mu<4 ; mu++ ) { (*Om).u1[mu] /= (double) Nsite; (*Om).J1[mu] /= (double) Nsite; 
    (*Om).uL1[mu] /= (double) Nsite; (*Om).JL1[mu] /= (double) Nsite; 
    (*Om).sdir[mu] /= (double) Nsite; (*Om).pi[mu] /= (double) Nsite; }
  for( munu=0 ; munu<6 ; munu++ ){
    (*Om).nvort[munu] /= (double) Nsite;
    (*Om).svort[munu] /= (double) Nsite;
  }


  /*  < T_mu,nu(x) T_mu',nu'(x') >  &  < Omega_mu,nu(x) Omega_mu',nu'(x') >  & 
      < xi_mu(x) xi_nu(x') >  &  < b(x) b(x') >  , etc. :  */

  /*  initialization  */
  //  j = params.Nx/2 + 1 + params.Nt;
  //  j = (params.Nx/2 + 1) * (params.Nt/2 + 1);
  j = params.Nt/2 + 1 + params.Nt;
  for(i=0;i<j;i++){
    for( munu=0 ; munu<10 ; munu++ ){
      for( munu2=0 ; munu2<10 ; munu2++ ){
	TmnTmn[i].dir1234[munu][munu2] = 0.0;
	OmOm[i].ndir1234[munu][munu2] = 0.0;
      }
    }
    for( mu=0 ; mu<4 ; mu++ )
      for( nu=0 ; nu<4 ; nu++ ){
	OmOm[i].uu12[mu][nu] = 0.0;
	OmOm[i].JJ12[mu][nu] = 0.0;
	OmOm[i].uuL12[mu][nu] = 0.0;
	OmOm[i].JJL12[mu][nu] = 0.0;
	OmOm[i].sdir12[mu][nu] = 0.0;
      }
    for( munu=0 ; munu<6 ; munu++ ) 
      for( munu2=0 ; munu2<6 ; munu2++ ){
	OmOm[i].nvort12[munu][munu2] = 0.0;
	OmOm[i].svort12[munu][munu2] = 0.0;
      }
  }

  j = 3 * (params.Nx/2 + 1) * (params.Nt/2 + 1) + params.Nt;
  for(i=0;i<j;i++){
    for( mu=0 ; mu<4 ; mu++ )
      for( nu=0 ; nu<4 ; nu++ ) PiPi[i].dir12[mu][nu] = 0.0;
    scsc[i].bb = 0.0; scsc[i].XX = 0.0;
    scsc[i].yy = 0.0; scsc[i].xixi = 0.0;
    scsc[i].nn = 0.0; scsc[i].rhorho = 0.0;
    scsc[i].pp = 0.0; scsc[i].TT = 0.0;
  }

  /*  correlators as a function of x-x', t-t'  */

  for(k=0,j=0;k<3;k++) for(t=0;t<=params.Nt/2;t++) for(i=0;i<=params.Nx/2;i++,j++){

    /*  separation along lattice k-axis, shifted by t  */
    switch(k){
    case 0:
      dr[0] = i; dr[1] = 0; dr[2] = 0; dr[3] = t;
    case 1:
      dr[0] = 0; dr[1] = i; dr[2] = 0; dr[3] = t;
    case 2:
      dr[0] = 0; dr[1] = 0; dr[2] = i; dr[3] = t;
    }

    corr1 = corr2 = corr3 = corr4 = corr5 = corr6 = corr7 = corr8 = 0.0;

    /*  scalar correlators  */
#pragma omp parallel for private(is,s,js,s2,o1,o2,n_rho_p_T) reduction(+:corr1,corr2,corr3,corr4,corr5,corr6,corr7,corr8)
    FORALLSITES(is,s) 
      js = new_index_from_dr(is,dr);
      s2 = &(lattice[js]);
      corr1 += (*s).b * (*s2).b;
      corr2 += (*s).X * (*s2).X;
      corr3 += (*s).y * (*s2).y;
      o1 = sqrt( (*s).X - (*s).y * (*s).y );
      o2 = sqrt( (*s2).X - (*s2).y * (*s2).y );
      corr4 += o1 * o2;
      EoS_match( s , n_rho_p_T );
      EoS_match( s2 , &(n_rho_p_T[4]) );
      corr5 += n_rho_p_T[0] * n_rho_p_T[4];
      corr6 += n_rho_p_T[1] * n_rho_p_T[5];
      corr7 += n_rho_p_T[2] * n_rho_p_T[6];
      corr8 += n_rho_p_T[3] * n_rho_p_T[7];
    ENDALLSITES 

    scsc[j].bb = corr1;
    scsc[j].XX = corr2;
    scsc[j].yy = corr3;
    scsc[j].xixi = corr4;
    scsc[j].nn = corr5;
    scsc[j].rhorho = corr6;
    scsc[j].pp = corr7;
    scsc[j].TT = corr8;

    /*  vector correlators  */
    for( mu=0 ; mu<4 ; mu++ )
      for( nu=0 ; nu<4 ; nu++ ){

	dO = 0.0;

#pragma omp parallel for private(is,s,js,s2,o1) reduction(+:dO)
	FORALLSITES(is,s) 
	  js = new_index_from_dr(is,dr);
	  s2 = &(lattice[js]);
	  /*
	  if( mu < 3 ) o1 = (*s).phi[mu];
	  else o1 = (*s).psi;
	  if( nu < 3 ) o1 *= (*s2).phi[nu];
	  else o1 *= (*s2).psi;
	  */
	  o1 = Om_x[is].pi[mu] * Om_x[js].pi[nu];
	  dO += o1;
	ENDALLSITES 

	PiPi[j].dir12[mu][nu] = dO;
      }

    /*  separation along t-axis only for the following (note: [j] --> [t])  */
    if( k==0 && i==0 ){

    for( mu=0 ; mu<4 ; mu++ )
      for( nu=0 ; nu<4 ; nu++ ){

	dO = 0.0;

#pragma omp parallel for private(is,s,js,s2,o1) reduction(+:dO)
	FORALLSITES(is,s) 
	  js = new_index_from_dr(is,dr);
	  s2 = &(lattice[js]);
	  dO += Om_x[is].sdir[mu] * Om_x[js].sdir[nu];
	ENDALLSITES 

	OmOm[t].sdir12[mu][nu] = dO;
      }

    for( mu=0 ; mu<4 ; mu++ )
      for( nu=0 ; nu<4 ; nu++ ){

	dO = dOO = 0.0;

#pragma omp parallel for private(is,s,js,s2) reduction(+:dO,dOO)
	FORALLSITES(is,s) 
	  js = new_index_from_dr(is,dr);
	  s2 = &(lattice[js]);
	  dO += Om_x[is].u1[mu] * Om_x[js].u1[nu];
	  dOO += Om_x[is].J1[mu] * Om_x[js].J1[nu];
	ENDALLSITES 

	OmOm[t].uu12[mu][nu] = dO;
	OmOm[t].JJ12[mu][nu] = dOO;

	dO = dOO = 0.0;

#pragma omp parallel for private(is,s,js,s2) reduction(+:dO,dOO)
	FORALLSITES(is,s) 
	  js = new_index_from_dr(is,dr);
	  s2 = &(lattice[js]);
	  dO += Om_x[is].uL1[mu] * Om_x[js].uL1[nu];
	  dOO += Om_x[is].JL1[mu] * Om_x[js].JL1[nu];
	ENDALLSITES 

	OmOm[t].uuL12[mu][nu] = dO;
	OmOm[t].JJL12[mu][nu] = dOO;
      }

    /*  tensor correlators  */
    for( munu=0 ; munu<10 ; munu++ ) 
      for( munu2=0 ; munu2<10 ; munu2++ ){

	dOO = dTT = 0.0;

#pragma omp parallel for private(is,s,js,s2) reduction(+:dTT,dOO)
	FORALLSITES(is,s) 
	  js = new_index_from_dr(is,dr);
	  s2 = &(lattice[js]);
	  dTT += Tmn_x[is].dir12[munu] * Tmn_x[js].dir12[munu2];
	  dOO += Om_x[is].ndir12[munu] * Om_x[js].ndir12[munu2];
	ENDALLSITES 

	TmnTmn[t].dir1234[munu][munu2] = dTT;
	OmOm[t].ndir1234[munu][munu2] = dOO;
      }

    for( munu=0 ; munu<6 ; munu++ ) 
      for( munu2=0 ; munu2<6 ; munu2++ ){

	dOO = dO = 0.0;

#pragma omp parallel for private(is,s,js,s2) reduction(+:dOO,dO)
	FORALLSITES(is,s) 
	  js = new_index_from_dr(is,dr);
	  s2 = &(lattice[js]);
	  dOO += Om_x[is].nvort[munu] * Om_x[js].nvort[munu2];
	  dO += Om_x[is].svort[munu] * Om_x[js].svort[munu2];
	ENDALLSITES 

	OmOm[t].nvort12[munu][munu2] = dOO;
	OmOm[t].svort12[munu][munu2] = dO;
      }

    }  // end if( k==0 && i==0 )

  }

  /*  correlators as a function of t-t':  */
  /*
  for(i=0;i<params.Nt;i++){
    j = i + params.Nx/2 + 1;

    //  separation along lattice t-axis 
    dr[0] = 0; dr[1] = 0; dr[2] = 0; dr[3] = i;

    corr1 = corr2 = corr3 = corr4 = corr5 = corr6 = corr7 = corr8 = 0.0;

    //  scalar correlators 
#pragma omp parallel for private(is,s,js,s2,o1,o2,n_rho_p_T) reduction(+:corr1,corr2,corr3,corr4,corr5,corr6,corr7,corr8)
    FORALLSITES(is,s) 
      js = new_index_from_dr(is,dr);
      s2 = &(lattice[js]);
      corr1 += (*s).b * (*s2).b;
      corr2 += (*s).X * (*s2).X;
      corr3 += (*s).y * (*s2).y;
      o1 = sqrt( (*s).X - (*s).y * (*s).y );
      o2 = sqrt( (*s2).X - (*s2).y * (*s2).y );
      corr4 += o1 * o2;
      EoS_match( s , n_rho_p_T );
      EoS_match( s2 , &(n_rho_p_T[4]) );
      corr5 += n_rho_p_T[0] * n_rho_p_T[4];
      corr6 += n_rho_p_T[1] * n_rho_p_T[5];
      corr7 += n_rho_p_T[2] * n_rho_p_T[6];
      corr8 += n_rho_p_T[3] * n_rho_p_T[7];
    ENDALLSITES 

    OmOm[j].bb = corr1;
    OmOm[j].XX = corr2;
    OmOm[j].yy = corr3;
    OmOm[j].xixi = corr4;
    OmOm[j].nn = corr5;
    OmOm[j].rhorho = corr6;
    OmOm[j].pp = corr7;
    OmOm[j].TT = corr8;

    //  vector correlators 
    for( mu=0 ; mu<4 ; mu++ )
      for( nu=0 ; nu<4 ; nu++ ){

	dO = dOO = 0.0;

#pragma omp parallel for private(is,s,js,s2,o1) reduction(+:dO,dOO)
	FORALLSITES(is,s) 
	  js = new_index_from_dr(is,dr);
	  s2 = &(lattice[js]);
	  dO += Om_x[is].sdir[mu] * Om_x[js].sdir[nu];
	  if( mu < 3 ) o1 = (*s).phi[mu];
	  else o1 = (*s).psi;
	  if( nu < 3 ) o1 *= (*s2).phi[nu];
	  else o1 *= (*s2).psi;
	  dOO += o1;
	ENDALLSITES 

	OmOm[j].sdir12[mu][nu] = dO;
	OmOm[j].pipi.dir12[mu][nu] = dOO;
      }

    //  tensor correlators 
    for( munu=0 ; munu<10 ; munu++ ) 
      for( munu2=0 ; munu2<10 ; munu2++ ){

	dOO = dTT = 0.0;

#pragma omp parallel for private(is,s,js,s2) reduction(+:dTT,dOO)
	FORALLSITES(is,s) 
	  js = new_index_from_dr(is,dr);
	  s2 = &(lattice[js]);
	  dTT += Tmn_x[is].dir12[munu] * Tmn_x[js].dir12[munu2];
	  dOO += Om_x[is].ndir12[munu] * Om_x[js].ndir12[munu2];
	ENDALLSITES 

	TmnTmn[j].dir1234[munu][munu2] = dTT;
	OmOm[j].ndir1234[munu][munu2] = dOO;
      }

    for( munu=0 ; munu<6 ; munu++ ) 
      for( munu2=0 ; munu2<6 ; munu2++ ){

	dOO = dO = 0.0;

#pragma omp parallel for private(is,s,js,s2) reduction(+:dOO,dO)
	FORALLSITES(is,s) 
	  js = new_index_from_dr(is,dr);
	  s2 = &(lattice[js]);
	  dOO += Om_x[is].nvort[munu] * Om_x[js].nvort[munu2];
	  dO += Om_x[is].svort[munu] * Om_x[js].svort[munu2];
	ENDALLSITES 

	OmOm[j].nvort12[munu][munu2] = dOO;
	OmOm[j].svort12[munu][munu2] = dO;
      }

  }
  */


  for( munu=0 ; munu<10 ; munu++ ){
    Tmn_src.dir12[munu] = 0.0;
    Om_src.ndir12[munu] = 0.0;
  }
  for( mu=0 ; mu<4 ; mu++ ){
    Om_src.sdir[mu] = 0.0;
    Om_src.u1[mu] = 0.0;
    Om_src.J1[mu] = 0.0;
    Om_src.uL1[mu] = 0.0;
    Om_src.JL1[mu] = 0.0;
    Om_src.pi[mu] = 0.0;
  }
  for( munu=0 ; munu<6 ; munu++ ){
    Om_src.nvort[munu] = 0.0;
    Om_src.svort[munu] = 0.0;
  }
  for( j=0 ; j<8 ; j++ ) Sclr_src[j] = 0.0;

  //  zero-momentum sources (at t=N_t/2) 
  t = params.Nt / 2;
  FORALLSPACE(is,s,t){

    //  scalar sources 
    Sclr_src[0] += (*s).b;
    Sclr_src[1] += (*s).X;
    Sclr_src[2] += (*s).y;
    Sclr_src[3] += sqrt( (*s).X - (*s).y * (*s).y );
    EoS_match( s , n_rho_p_T );
    Sclr_src[4] += n_rho_p_T[0];
    Sclr_src[5] += n_rho_p_T[1];
    Sclr_src[6] += n_rho_p_T[2];
    Sclr_src[7] += n_rho_p_T[3];

    //  vector sources 
    for( mu=0 ; mu<4 ; mu++ ){
      Om_src.sdir[mu] += Om_x[is].sdir[mu];
      Om_src.u1[mu] += Om_x[is].u1[mu];
      Om_src.J1[mu] += Om_x[is].J1[mu];
      Om_src.uL1[mu] += Om_x[is].uL1[mu];
      Om_src.JL1[mu] += Om_x[is].JL1[mu];
    }
    for( mu=0 ; mu<3 ; mu++ ) Om_src.pi[mu] += (*s).phi[mu];
    Om_src.pi[3] += (*s).psi;

    //  tensor sources 
    for( munu=0 ; munu<10 ; munu++ ){
      Tmn_src.dir12[munu] += Tmn_x[is].dir12[munu];
      Om_src.ndir12[munu] += Om_x[is].ndir12[munu];
    }
    for( munu=0 ; munu<6 ; munu++ ){
      Om_src.nvort[munu] += Om_x[is].nvort[munu];
      Om_src.svort[munu] += Om_x[is].svort[munu];
    }
  }

  //  zero-momentum correlators 
  j = 3 * (params.Nx/2 + 1) * (params.Nt/2 + 1);
  for( t=0 ; t<params.Nt ; t++ , j++ ){
    //    j = t + params.Nx/2 + 1;
    k = t + params.Nt/2 + 1;

    FORALLSPACE(is,s,t){

      //  scalar correlators 
      scsc[j].bb += Sclr_src[0] * (*s).b;
      scsc[j].XX += Sclr_src[1] * (*s).X;
      scsc[j].yy += Sclr_src[2] * (*s).y;
      o2 = sqrt( (*s).X - (*s).y * (*s).y );
      scsc[j].xixi += Sclr_src[3] * o2;
      EoS_match( s , n_rho_p_T );
      scsc[j].nn += Sclr_src[4] * n_rho_p_T[0];
      scsc[j].rhorho += Sclr_src[5] * n_rho_p_T[1];
      scsc[j].pp += Sclr_src[6] * n_rho_p_T[2];
      scsc[j].TT += Sclr_src[7] * n_rho_p_T[3];

      //  vector correlators 
      for( mu=0 ; mu<4 ; mu++ )
	for( nu=0 ; nu<4 ; nu++ ){
	  dOO = Om_src.sdir[mu] * Om_x[is].sdir[nu];
	  OmOm[k].sdir12[mu][nu] += dOO;
	  dOO = Om_src.pi[mu];
	  if( nu < 3 ) dOO *= (*s).phi[nu];
	  else dOO *= (*s).psi;
	  PiPi[j].dir12[mu][nu] += dOO;
	}

      //  tensor correlators 
      for( munu=0 ; munu<10 ; munu++ ) 
	for( munu2=0 ; munu2<10 ; munu2++ ){
	  dTT = Tmn_src.dir12[munu] * Tmn_x[is].dir12[munu2];
	  TmnTmn[k].dir1234[munu][munu2] += dTT;
	  dOO = Om_src.ndir12[munu] * Om_x[is].ndir12[munu2];
	  OmOm[k].ndir1234[munu][munu2] += dOO;
	}
      for( munu=0 ; munu<6 ; munu++ ) 
	for( munu2=0 ; munu2<6 ; munu2++ ){
	  dOO = Om_src.nvort[munu] * Om_x[is].nvort[munu2];
	  OmOm[k].nvort12[munu][munu2] += dOO;
	  dOO = Om_src.svort[munu] * Om_x[is].svort[munu2];
	  OmOm[k].svort12[munu][munu2] += dOO;
	}

    }

  }


  //  j = params.Nx/2 + 1 + params.Nt;
  //  j = (params.Nx/2 + 1) * (params.Nt/2 + 1);

  j = params.Nt/2 + 1 + params.Nt;
  for(i=0;i<j;i++){
    if( i <= params.Nt/2 ) o1 = 1.0 / (double) Nsite;
    else o1 = 1.0 / ( (double) (Nspace*Nspace) );
    for( munu=0 ; munu<10 ; munu++ ){
      for( munu2=0 ; munu2<10 ; munu2++ ){
	TmnTmn[i].dir1234[munu][munu2] *= o1;
	OmOm[i].ndir1234[munu][munu2] *= o1;
      }
    }
    for( munu=0 ; munu<6 ; munu++ ) 
      for( munu2=0 ; munu2<6 ; munu2++ ){
	OmOm[i].nvort12[munu][munu2] *= o1;
	OmOm[i].svort12[munu][munu2] *= o1;
      }
    for( mu=0 ; mu<4 ; mu++ ) 
      for( nu=0 ; nu<4 ; nu++ ){
	OmOm[i].uu12[mu][nu] *= o1;
	OmOm[i].JJ12[mu][nu] *= o1; 
	OmOm[i].uuL12[mu][nu] *= o1;
	OmOm[i].JJL12[mu][nu] *= o1; 
	OmOm[i].sdir12[mu][nu] *= o1;
      }
  }

  j = 3 * (params.Nx/2 + 1) * (params.Nt/2 + 1);
  for(i=0;i<j;i++){
    o1 = 1.0 / (double) Nsite;
    for( mu=0 ; mu<4 ; mu++ ) 
      for( nu=0 ; nu<4 ; nu++ ) 
	PiPi[i].dir12[mu][nu] *= o1;
    scsc[i].bb *= o1; scsc[i].XX *= o1;
    scsc[i].yy *= o1; scsc[i].xixi *= o1;
    scsc[i].nn *= o1; scsc[i].rhorho *= o1;
    scsc[i].pp *= o1; scsc[i].TT *= o1;
  }

  k = 3 * (params.Nx/2 + 1) * (params.Nt/2 + 1) + params.Nt;
  for(i=j;i<k;i++){
    o1 = 1.0 / ( (double) (Nspace*Nspace) );
    for( mu=0 ; mu<4 ; mu++ ) 
      for( nu=0 ; nu<4 ; nu++ ) 
	PiPi[i].dir12[mu][nu] *= o1;
    scsc[i].bb *= o1; scsc[i].XX *= o1;
    scsc[i].yy *= o1; scsc[i].xixi *= o1;
    scsc[i].nn *= o1; scsc[i].rhorho *= o1;
    scsc[i].pp *= o1; scsc[i].TT *= o1;
  }

  free(Om_x);
  free(Tmn_x);
  return;
}



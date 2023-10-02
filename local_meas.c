
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "lattice.h"
#include "macros.h"
#include "ltools.h"
#include "local_meas.h"

#define TINY 1.e-12
#define SMALL 1.e-5


/*  calculate 1st derivs of psi and phi's, and associated quantities  */
void field_derivs( ){
  register int is,i,j,k,l,mu,x,t,js;
  register site *s,*sp,*sm,*sn;
  double BIJ[3][3],detB,tmp[4],dp[3],b2;
  int dr[4];

#ifdef DUAL_LATTICE
#pragma omp parallel for private(is,s,mu,i)
  FORALLSITES(is,s) 
    for(mu=0;mu<4;mu++){
      if( mu < 3 ) (*s).dpsi[mu] = 0.0;
      else (*s).dpsi[mu] = params.chem_pot;
      for(i=0;i<3;i++){
	if( i == mu ) (*s).dphi[i][mu] = 1.0;
	else (*s).dphi[i][mu] = 0.0;
      }
    }
  ENDALLSITES
#endif

#pragma omp parallel for private(is,s,mu,i,j,sn,sm,sp)
  FORALLSITES(is,s) 
    /*  phi^I(x) and psi(x) field derivs  */
#ifdef DUAL_LATTICE
    for(mu=0;mu<4;mu++){
      for(j=0;j<16;j++){
	sn = (site *) (*s).dneibBC[j];
	(*s).dpsi[mu] += dsign[mu][j] * 0.125 * (*sn).psi;
	for(i=0;i<3;i++) (*s).dphi[i][mu] += dsign[mu][j] * 0.125 * (*sn).phi[i];
      }
    }
#else
    for(mu=0;mu<4;mu++){
      sp = (site *)((*s).neib[mu]);
      sm = (site *)((*s).neib[OPP_DIR(mu)]);
      (*s).dpsi[mu] = 0.5 * ( (*sp).psi - (*sm).psi );
      if( mu == 3 ) (*s).dpsi[mu] += params.chem_pot;
      if( params.upd_phi ) 
	for(i=0;i<3;i++){
	  (*s).dphi[i][mu] = 0.5 * ( (*sp).phi[i] - (*sm).phi[i] );
	  if( i == mu ) (*s).dphi[i][mu] += 1.0;
	}
    }
#endif
  ENDALLSITES

#pragma omp parallel for private(is,s,mu,i,j,k,BIJ,tmp,dp,detB,b2)
  FORALLSITES(is,s) 

    /*  scalar invariant: X = d_m psi d^m psi  */
    (*s).X = 0.0;
    for(mu=0;mu<4;mu++) (*s).X += (*s).dpsi[mu] * (*s).dpsi[mu];

    if( params.upd_phi ){
      /*  B^IJ = d_mu phi^I d^mu phi^J  */
      for(i=0;i<3;i++) 
	for(j=i;j<3;j++){
	  BIJ[i][j] = 0.0;
	  for(mu=0;mu<4;mu++) BIJ[i][j] += (*s).dphi[i][mu] * (*s).dphi[j][mu];
	  if(j>i) BIJ[j][i] = BIJ[i][j];
	}

      /*
	the determinant can be very small (locations of very low entropy... 
	centers of vortices?) and we must be more careful about how we calculate 
	(or what we use instead of) the inverse of B_IJ ... 
	perhaps we should switch to something like  A_I,mu = B_IJ^-1 * d_mu phi^J 
	by solving B_IJ A_J,mu = d_mu phi^I: 
      */
      tmp[0] = -BIJ[1][0]/BIJ[0][0];
      BIJ[1][1] += tmp[0]*BIJ[0][1];
      BIJ[1][2] += tmp[0]*BIJ[0][2];
      tmp[1] = -BIJ[2][0]/BIJ[0][0];
      BIJ[2][1] += tmp[1]*BIJ[0][1];
      BIJ[2][2] += tmp[1]*BIJ[0][2];
      tmp[2] = -BIJ[2][1]/BIJ[1][1];
      BIJ[2][2] += tmp[2]*BIJ[1][2];
      for(mu=0;mu<4;mu++){
	dp[0] = (*s).dphi[0][mu];
	dp[1] = (*s).dphi[1][mu] + tmp[0] * dp[0];
	dp[2] = (*s).dphi[2][mu] + tmp[1] * dp[0] + tmp[2] * dp[1];
	(*s).AImu[2][mu] = dp[2] / BIJ[2][2];
	(*s).AImu[1][mu] = ( dp[1] - BIJ[1][2] * (*s).AImu[2][mu] ) / BIJ[1][1];
	(*s).AImu[0][mu] = ( dp[0] - BIJ[0][1] * (*s).AImu[1][mu] - 
			     BIJ[0][2] * (*s).AImu[2][mu] ) / BIJ[0][0];
      }

      /*  fast inversion of 3x3 real symmetric matrix, B^IJ  */
      /*
      tmp[0] = BIJ[1][1]*BIJ[2][2] - BIJ[1][2]*BIJ[1][2];
      tmp[1] = BIJ[0][2]*BIJ[1][2] - BIJ[0][1]*BIJ[2][2];
      tmp[2] = BIJ[0][1]*BIJ[1][2] - BIJ[1][1]*BIJ[0][2];
      */
      /*  determinant  */
      //      detB = BIJ[0][0] * tmp[0] + BIJ[0][1] * tmp[1] + BIJ[0][2] * tmp[2];
      /*
      if( fabs(detB) < TINY ){
	printf("singular matrix? det(B^IJ) = %le\n",detB);
	exit(0);
      }
      */
      /*  inverse  */
      /*
      (*s).Binv[0][0] = tmp[0]/detB;
      (*s).Binv[0][1] = tmp[1]/detB;
      (*s).Binv[0][2] = tmp[2]/detB;
      tmp[0] = BIJ[0][0]*BIJ[2][2] - BIJ[0][2]*BIJ[0][2];
      tmp[1] = BIJ[0][1]*BIJ[0][2] - BIJ[0][0]*BIJ[1][2];
      tmp[2] = BIJ[0][0]*BIJ[1][1] - BIJ[0][1]*BIJ[0][1];
      (*s).Binv[1][1] = tmp[0]/detB;
      (*s).Binv[1][2] = tmp[1]/detB;
      (*s).Binv[2][2] = tmp[2]/detB;
      */
      /*  redundant entries  */
      /*
      (*s).Binv[1][0] = (*s).Binv[0][1];
      (*s).Binv[2][0] = (*s).Binv[0][2];
      (*s).Binv[2][1] = (*s).Binv[1][2];
      */

      /*  the (normal) fluid velocity: u^mu = J^mu / b  ; 
	  J^m = eps^mabc eps_IJK  d_a phi^I  d_b phi^J  d_c phi^K  ; 
	  scalar invariant: b^2 = J_m J^m = det B^IJ  */
      b2 = 0.0;
      for(mu=0;mu<4;mu++){
	(*s).uu[mu] = 0.0;
	i = (mu+1)%4;
	j = (mu+2)%4;
	k = (mu+3)%4;
	(*s).uu[mu] += (*s).dphi[0][i] * (*s).dphi[1][j] * (*s).dphi[2][k];
	(*s).uu[mu] -= (*s).dphi[0][i] * (*s).dphi[1][k] * (*s).dphi[2][j];
	(*s).uu[mu] += (*s).dphi[0][j] * (*s).dphi[1][k] * (*s).dphi[2][i];
	(*s).uu[mu] -= (*s).dphi[0][j] * (*s).dphi[1][i] * (*s).dphi[2][k];
	(*s).uu[mu] += (*s).dphi[0][k] * (*s).dphi[1][i] * (*s).dphi[2][j];
	(*s).uu[mu] -= (*s).dphi[0][k] * (*s).dphi[1][j] * (*s).dphi[2][i];
	b2 += (*s).uu[mu] * (*s).uu[mu];
      }
      /*
      tmp[0] = (detB - b2) / detB;
      if( fabs(tmp[0]) > SMALL ){
	printf("detB= %le  !=  b^2= %le\n",detB,b2);
	exit(0);
      }
      */
      (*s).b = sqrt(b2);
      for(mu=0;mu<4;mu++) (*s).uu[mu] /= (*s).b;
    }
    else{
      for(i=0;i<3;i++) 
	for(mu=0;mu<4;mu++){
	  if( i == mu ) (*s).dphi[i][mu] = 1.0;
	  else (*s).dphi[i][mu] = 0.0;
	}
      //      for(mu=0;mu<3;mu++) (*s).uu[mu] = 0.0;
      //      (*s).uu[3] = 1.0;
      (*s).b = 1.0;
      for(mu=0;mu<4;mu++) (*s).uu[mu] = (*s).dpsi[mu] / sqrt( (*s).X ) ;
    }

    /*  scalar invariant: y = u^m d_m psi  */
    (*s).y = 0.0;
    for(mu=0;mu<4;mu++) (*s).y += (*s).uu[mu] * (*s).dpsi[mu];

  ENDALLSITES


  if( params.bc == TIME_MixBC ){
    for(i=0;i<3;i++) dr[i] = 0;
    dr[3] = -1;
    x = 0;
    FORALLSPACE(is,s,0){
      js = new_index_from_dr(is,dr);
      sm = &(lattice[js]);
      for(i=0;i<3;i++){
	//	for(j=0;j<3;j++) 
	//	  latt_bound[x].Binv[i][j] = (*sm).Binv[i][j];
	for(mu=0;mu<4;mu++) 
	  latt_bound[x].dphi[i][mu] = (*sm).dphi[i][mu];
      }
      latt_bound[x].b = (*sm).b;
      latt_bound[x].X = (*sm).X;
      latt_bound[x].y = (*sm).y;
      x++;
    }
    dr[3] = 1;
    t = params.Nt - 1;
    FORALLSPACE(is,s,t){
      js = new_index_from_dr(is,dr);
      sp = &(lattice[js]);
      for(i=0;i<3;i++){
	//	for(j=0;j<3;j++) 
	//	  latt_bound[x].Binv[i][j] = (*sp).Binv[i][j];
	for(mu=0;mu<4;mu++) 
	  latt_bound[x].dphi[i][mu] = (*sp).dphi[i][mu];
      }
      latt_bound[x].b = (*sp).b;
      latt_bound[x].X = (*sp).X;
      latt_bound[x].y = (*sp).y;
      x++;
    }
  }

  return;
}


double deriv( field_offset fo , site *s , int mu ){
  register int i;
  register site *sp,*sm;
  register double *qp,*qm;
  double dq;

#ifdef DUAL_LATTICE
  dq = 0.0;
  for(i=0;i<16;i++){
    sm = (site *) (*s).dneib[i];
    qm = (double *)F_PT(sm,fo);
    dq -= dsign[mu][i] * 0.125 * qm[0];
  }
#else
  sp = (site *)((*s).neib[mu]);
  sm = (site *)((*s).neib[OPP_DIR(mu)]);
  qp = (double *)F_PT(sp,fo);
  qm = (double *)F_PT(sm,fo);
  dq = 0.5 * ( qp[mu] - qm[mu] );
#endif

  return( dq );
}


double derivBC( field_offset fo , site *s , int mu ){
  register int i;
  register site *sp,*sm;
  register double *qp,*qm;
  double dq;

#ifdef DUAL_LATTICE
  dq = 0.0;
  for(i=0;i<16;i++){
    sp = (site *) (*s).dneibBC[i];
    qp = (double *)F_PT(sp,fo);
    dq += dsign[mu][i] * 0.125 * qp[0];
  }
#else
  sp = (site *)((*s).neib[mu]);
  sm = (site *)((*s).neib[OPP_DIR(mu)]);
  qp = (double *)F_PT(sp,fo);
  qm = (double *)F_PT(sm,fo);
  dq = 0.5 * ( qp[mu] - qm[mu] );
#endif

  return( dq );
}


double avgBC( field_offset fo , site *s ){
  register int i,mu;
  register site *sp,*sm;
  register double *qp,*qm;
  double aq;

  aq = 0.0;
#ifdef DUAL_LATTICE
  for(i=0;i<16;i++){
    sp = (site *) (*s).dneibBC[i];
    qp = (double *)F_PT(sp,fo);
    aq += 0.0625 * qp[0];
  }
#else
  for(mu=0;mu<4;mu++){
    sp = (site *)((*s).neib[mu]);
    sm = (site *)((*s).neib[OPP_DIR(mu)]);
    qp = (double *)F_PT(sp,fo);
    qm = (double *)F_PT(sm,fo);
    aq += 0.125 * ( qp[mu] + qm[mu] );
  }
#endif

  return( aq );
}


/*  divergence of entropy current: div.J(x)  */
double divJ( site *s ){
  register int mu,i;
  register site *sp,*sm;
  double dJ;

  dJ = 0.0;
  for(mu=0;mu<4;mu++){
#ifdef DUAL_LATTICE
    for(i=0;i<16;i++){
      sm = (site *) (*s).dneib[i];
      dJ -= dsign[mu][i] * 0.125 * (*sm).uu[mu] * (*sm).b;
    }
#else
    sp = (site *)((*s).neib[mu]);
    sm = (site *)((*s).neib[OPP_DIR(mu)]);
    dJ += 0.5 * ( (*sp).uu[mu] * (*sp).b - (*sm).uu[mu] * (*sm).b );
#endif
  }

  return(dJ);
}


/*  local vorticity: d_mu u_nu - d_nu u_mu  */
void vorticity( field_offset velFO , site *s , double *vort ){
  register int mu,nu,munu,i;
  register site *sp,*sm;
  register double *vp,*vm,u2,gp,gm,dv;

  for(mu=0,munu=0;mu<3;mu++){
#ifdef DUAL_LATTICE
    /*  Euclidean or Lorentzian version, V_mu^nu:  */
    for(nu=mu+1;nu<4;nu++,munu++){
      vort[munu] = 0.0;
      for(i=0;i<16;i++){
	sm = (site *)((*s).dneib[i]);
	vm = (double *)F_PT(sm,velFO);
	//	gm = gammaL(sm);  // for normal fluid only!!!
	//	u2 = vm[0]*vm[0] + vm[1]*vm[1] + vm[2]*vm[2] + vm[3]*vm[3];
	//	gm = sqrt(u2)/fabs(vm[3]);
	dv = dsign[mu][i] * 0.125 * vm[nu];
	//	if( nu < 3 ) dv = dsign[mu][i] * 0.125 * vm[nu] * gm;
	//	else dv = dsign[mu][i] * 0.125 * gm;
	vort[munu] -= dv;
	dv = dsign[nu][i] * 0.125 * vm[mu];
	//	if( mu < 3 ) dv = dsign[nu][i] * 0.125 * vm[mu] * gm;
	//	else dv = dsign[nu][i] * 0.125 * gm;
	vort[munu] += dv;
      }
    }
#else
    /*  still Euclidean version here:  */
    sp = (site *)((*s).neib[mu]);
    sm = (site *)((*s).neib[OPP_DIR(mu)]);
    vp = (double *)F_PT(sp,velFO);
    vm = (double *)F_PT(sm,velFO);
    for(nu=mu+1;nu<4;nu++,munu++){
      vort[munu] = 0.5 * ( vp[nu] - vm[nu] );
      sp = (site *)((*s).neib[nu]);
      sm = (site *)((*s).neib[OPP_DIR(nu)]);
      vp = (double *)F_PT(sp,velFO);
      vm = (double *)F_PT(sm,velFO);
      vort[munu] -= 0.5 * ( vp[mu] - vm[mu] );
    }
#endif
  }

  return;
}


/*  local fluid boost factor  */
double gammaL( site *s ){
  //  double v2,g;
  //  v2 = 1.0 - (*s).uu[3] * (*s).uu[3];
  //  g = 1.0 / sqrt(1.0 - v2);
  //  return(g);
  return( 1.0 / fabs( (*s).uu[3] ) );
}


/*  returns the local, contravariant, Lorentzian 4-velocity (note: 3 remains the t direction)  */
void eucl_to_lrnz_vel( site *s , double *uu_L ){
  double g;
  g = gammaL(s);
  uu_L[0] = (*s).uu[0] * g;
  uu_L[1] = (*s).uu[1] * g;
  uu_L[2] = (*s).uu[2] * g;
  uu_L[3] = g;
  return;
}


void lrnz_flow_match( site *s , flow *Om ){
  int mu,nu,munu;
  double uu_L[4];
  eucl_to_lrnz_vel( s , uu_L );
  /*  Omega^mu,nu = g^mu,nu - u^mu*u^nu  */
  for( mu=0,munu=0 ; mu<4 ; mu++ ) 
    for( nu=mu ; nu<4 ; nu++,munu++ ){
      (*Om).ndir12[munu] = -uu_L[mu] * uu_L[nu];
      if( nu == mu && nu < 3 ) (*Om).ndir12[munu] -= 1;
      else if( nu == mu && nu == 3 ) (*Om).ndir12[munu] += 1;
    }
  return;
}



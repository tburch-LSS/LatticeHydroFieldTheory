
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "lattice.h"
#include "macros.h"
#include "ranlxd.h"
#include "ltools.h"
#include "io.h"
#include "action.h"
#include "setup.h"


/*  */
void setup( ){
  register int x,y,z,t,is,i,j,k,mu,nu;
  register site *s,*sm,*smm,*smmm,*smmmm;

  /*  set up neighbor pointer arrays 
      (4D torus; t is slowest running, x fastest)  */
  s = &(lattice[0]);
  for( t=0,is=0 ; t<params.Nt ; t++) 
    for( z=0 ; z<params.Nz ; z++ ) 
      for( y=0 ; y<params.Ny ; y++ ) 
	for( x=0 ; x<params.Nx ; x++,is++,s++ ) {
	  //	  (*s).xyzt[0] = x; (*s).xyzt[1] = y;
	  //	  (*s).xyzt[2] = z; (*s).xyzt[3] = t;
	  (*s).index = is;

	  if( x>0 ) (*s).neib[XDOWN] = (char *)(s - 1);
	  else (*s).neib[XDOWN] = (char *)(s + params.Nx - 1);
	  if( x<(params.Nx-1) )	(*s).neib[XUP] = (char *)(s + 1);
	  else (*s).neib[XUP] = (char *)(s - params.Nx + 1);

	  if( y>0 ) (*s).neib[YDOWN] = (char *)(s - params.Nx);
	  else (*s).neib[YDOWN] = (char *)(s + params.Nx*(params.Ny-1));
	  if( y<(params.Ny-1) )	(*s).neib[YUP] = (char *)(s + params.Nx);
	  else (*s).neib[YUP] = (char *)(s - params.Nx*(params.Ny-1));

	  if( z>0 ) (*s).neib[ZDOWN] = (char *)(s - params.Vol2);
	  else (*s).neib[ZDOWN] = (char *)(s + params.Vol2*(params.Nz-1));
	  if( z<(params.Nz-1) )	(*s).neib[ZUP] = (char *)(s + params.Vol2);
	  else (*s).neib[ZUP] = (char *)(s - params.Vol2*(params.Nz-1));

	  if( t>0 ) (*s).neib[TDOWN] = (char *)(s - params.Vol3);
	  else (*s).neib[TDOWN] = (char *)(s + params.Vol4 - params.Vol3);
	  if( t<(params.Nt-1) )	(*s).neib[TUP] = (char *)(s + params.Vol3);
	  else (*s).neib[TUP] = (char *)(s - params.Vol4 + params.Vol3);
	}

#ifdef DUAL_LATTICE
#pragma omp parallel for private(is,s,sm,smm,smmm,smmmm,mu,i,j,k)
  FORALLSITES(is,s) 
    i = 0;
    (*s).dneib[i] = (char *) s;
    for( mu=0 ; mu<4 ; mu++ ){
      sm = (site *)((*s).neib[OPP_DIR(mu)]);
      i++;
      (*s).dneib[i] = (char *) sm;
      for( j=mu+1 ; j<4 ; j++ ){
	smm = (site *)((*sm).neib[OPP_DIR(j)]);
	i++;
	(*s).dneib[i] = (char *) smm;
	for( k=j+1 ; k<4 ; k++ ){
	  smmm = (site *)((*smm).neib[OPP_DIR(k)]);
	  i++;
	  (*s).dneib[i] = (char *) smmm;
	  if( k == ZUP ){
	    smmmm = (site *)((*smmm).neib[TDOWN]);
	    i++;
	    (*s).dneib[i] = (char *) smmmm;
	  }
	}
      }
    }
  ENDALLSITES 

  for( nu=0 ; nu<4 ; nu++ ){
    i = 0;
    dsign[nu][i] = -1;
    for( mu=0 ; mu<4 ; mu++ ){
      i++;
      if( nu == mu ) dsign[nu][i] = 1;
      else dsign[nu][i] = -1;
      for( j=mu+1 ; j<4 ; j++ ){
	i++;
	if( nu == j || nu == mu ) dsign[nu][i] = 1;
	else dsign[nu][i] = -1;
	for( k=j+1 ; k<4 ; k++ ){
	  i++;
	  if( nu == k || nu == j || nu == mu ) dsign[nu][i] = 1;
	  else dsign[nu][i] = -1;
	  if( k == ZUP ){
	    i++;
	    dsign[nu][i] = 1;
	    //	    if( nu == TUP ) dsign[nu][i] = 1;
	    //	    else dsign[nu][i] = -1;
	  }
	}
      }
    }
    //    for(i=0;i<16;i++) printf("dsign[%2d][%2d]= %d\n",nu,i,dsign[nu][i]);
  }

#pragma omp parallel for private(is,s,sm,smm,smmm,smmmm,mu,i,j,k)
  FORALLSITES(is,s) 
    i = 0;
    (*s).dneibBC[i] = (char *) s;
    for( mu=0 ; mu<4 ; mu++ ){
      sm = (site *)((*s).neib[mu]);
      i++;
      (*s).dneibBC[i] = (char *) sm;
      for( j=mu+1 ; j<4 ; j++ ){
	smm = (site *)((*sm).neib[j]);
	i++;
	(*s).dneibBC[i] = (char *) smm;
	for( k=j+1 ; k<4 ; k++ ){
	  smmm = (site *)((*smm).neib[k]);
	  i++;
	  (*s).dneibBC[i] = (char *) smmm;
	  if( k == ZUP ){
	    smmmm = (site *)((*smmm).neib[TUP]);
	    i++;
	    (*s).dneibBC[i] = (char *) smmmm;
	  }
	}
      }
    }
  ENDALLSITES 
#endif

  /*  Neumann or Mixed BC's in time  */
  if( params.bc != PERIODIC ){
    k = 0;
    FORALLSPACE(is,s,0){
      (*s).neib[TDOWN] = (char *) &(latt_bound[k]);
      k++;
    }
    t = params.Nt - 1;
    FORALLSPACE(is,s,t){
      (*s).neib[TUP] = (char *) &(latt_bound[k]);
      k++;
    }
  }
  /*  Neumann BC's everywhere  */
  if( params.bc == ALL_NeuBC ){
    FORALLXBOUND(is,s,0){
      (*s).neib[XDOWN] = (char *) &(latt_bound[k]);
      k++;
    }
    x = params.Nx - 1;
    FORALLXBOUND(is,s,x){
      (*s).neib[XUP] = (char *) &(latt_bound[k]);
      k++;
    }
    FORALLYBOUND(is,s,j,0){
      (*s).neib[YDOWN] = (char *) &(latt_bound[k]);
      k++;
    }
    y = params.Ny - 1;
    FORALLYBOUND(is,s,j,y){
      (*s).neib[YUP] = (char *) &(latt_bound[k]);
      k++;
    }
    FORALLZBOUND(is,s,j,0){
      (*s).neib[ZDOWN] = (char *) &(latt_bound[k]);
      k++;
    }
    z = params.Nz - 1;
    FORALLZBOUND(is,s,j,z){
      (*s).neib[ZUP] = (char *) &(latt_bound[k]);
      k++;
    }
  }

  return;
}


/*  check neighbours  */
void check_setup( ){
  register int is,dir,i,js;
  register site *s,*sp,*sm;
  int dr[4],xyzt[4],xyzt2[4];

  FORALLSITES(is,s) 
    xyzt_from_index(xyzt,is);
    printf("isite= %d : x= %d  y= %d  z= %d  t= %d\n",
	   is,xyzt[0],xyzt[1],xyzt[2],xyzt[3]);
    for(dir=0;dir<4;dir++){
      for(i=0;i<4;i++) dr[i] = 0;
      dr[dir] = 1;
      js = new_index_from_dr(is,dr);
      xyzt_from_index(xyzt2,js);
      printf(" dir= %d : xn= %d  yn= %d  zn= %d  tn= %d\n",
	     dir,xyzt2[0],xyzt2[1],xyzt2[2],xyzt2[3]);
      sp = (site *) (*s).neib[dir];
      if( js != (*sp).index ) 
	fprintf(stderr,"neib index mismatch at %d, dir=%d: %d != %d\n",
		is,dir,js,(*sp).index);
      dr[dir] = -1;
      js = new_index_from_dr(is,dr);
      xyzt_from_index(xyzt2,js);
      printf(" dir= %d : xn= %d  yn= %d  zn= %d  tn= %d\n",
	     OPP_DIR(dir),xyzt2[0],xyzt2[1],xyzt2[2],xyzt2[3]);
      sm = (site *) (*s).neib[OPP_DIR(dir)];
      if( js != (*sm).index ) 
	fprintf(stderr,"neib index mismatch at %d, dir=%d: %d != %d\n",
		is,OPP_DIR(dir),js,(*sm).index);
    }
  ENDALLSITES 

  return;
}


/*  initialize phi_i(x) and psi(x) fields  */
void init_fields( ){
  register int is,i;
  register site *s;
  int xyzt[4];
  double rndm[4],phi0,psi0;

  phi0 = pow( params.b0 , 0.33333333333333 );
  psi0 = params.y0;

  /*  "cold" start  */
  if( params.start == 0 ){
#pragma omp parallel for private(is,s,i)
    FORALLSITES(is,s) 
      //      xyzt_from_index(xyzt,is);
      for(i=0;i<3;i++) (*s).phi[i] = 0.0;
      (*s).psi = 0.0;
    ENDALLSITES 
  }
  /*  file start  */
  else if( params.start == 1 ){
    read_lattice();
  }
  /*  "hot" start  */
  else if( params.start == 2 ){
    /*  #pragma omp parallel for private(is,s,i,rndm)  */
    FORALLSITES(is,s) 
      //      for(i=0;i<4;i++) rndm[i] = drand48();

      /*  #pragma omp critical  */
      /*  #pragma omp ordered  */
      ranlxd( rndm , 4 );

      for(i=0;i<3;i++) (*s).phi[i] = phi0 * rndm[i];
      (*s).psi = psi0 * rndm[3];
    ENDALLSITES 
  }
  else{
    printf("params.start = %d not defined!!! (0=cold,1=file,2=hot)\n",params.start);
    exit(0);
  }

  /*  initialize boundary fields (set according to action)  */
  if( params.bc != PERIODIC ) init_bound();

  /*  calculate initial field derivs  */
  field_derivs();

  return;
}



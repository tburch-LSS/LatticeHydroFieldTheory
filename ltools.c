
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lattice.h"
#include "macros.h"
#include "ltools.h"


/*  lattice array index from x,y,z,t coords  */
int index_from_xyzt( int *xyzt ){
  int index;

  index = xyzt[0] + xyzt[1]*params.Nx + xyzt[2]*params.Vol2 + 
    xyzt[3]*params.Vol3;

  return(index);
}


/*  coords from lattice array index  */
void xyzt_from_index( int *xyzt, int index ){

  xyzt[0] = index % params.Nx;
  xyzt[1] = (index % params.Vol2) / params.Nx;
  xyzt[2] = (index % params.Vol3) / params.Vol2;
  xyzt[3] = index / params.Vol3;

  return;
}


/*  new lattice array index from shift in coords (dr)  */
int new_index_from_dr( int index, int *dr ){
  int xyzt[4],j,ni;

  xyzt_from_index(xyzt,index);
  for(j=0;j<4;j++) xyzt[j] += dr[j];
  if( xyzt[0] < 0 ) xyzt[0] += params.Nx;
  if( xyzt[0] > params.Nx-1 ) xyzt[0] -= params.Nx;
  if( xyzt[1] < 0 ) xyzt[1] += params.Ny;
  if( xyzt[1] > params.Ny-1 ) xyzt[1] -= params.Ny;
  if( xyzt[2] < 0 ) xyzt[2] += params.Nz;
  if( xyzt[2] > params.Nz-1 ) xyzt[2] -= params.Nz;
  if( xyzt[3] < 0 ) xyzt[3] += params.Nt;
  if( xyzt[3] > params.Nt-1 ) xyzt[3] -= params.Nt;
  ni = index_from_xyzt(xyzt);

  return(ni);
}


/*  lattice array index from coords and shift (dr)  */
int index_from_xyzt_dr( int *xyzt, int *dr ){
  int r[4],j,in;

  for(j=0;j<4;j++) r[j] = xyzt[j] + dr[j];
  if( r[0] < 0 ) r[0] += params.Nx;
  if( r[0] > params.Nx-1 ) r[0] -= params.Nx;
  if( r[1] < 0 ) r[1] += params.Ny;
  if( r[1] > params.Ny-1 ) r[1] -= params.Ny;
  if( r[2] < 0 ) r[2] += params.Nz;
  if( r[2] > params.Nz-1 ) r[2] -= params.Nz;
  if( r[3] < 0 ) r[3] += params.Nt;
  if( r[3] > params.Nt-1 ) r[3] -= params.Nt;
  in = index_from_xyzt(r);

  return(in);
}



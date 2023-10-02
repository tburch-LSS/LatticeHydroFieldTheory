
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int setdim();
double combine_avg( double *atmp );
double combine_corr( int LL , int TT , int x , double *atmp , double *ctmp );
void write_corr( int LL , int TT , int n , double *Corr , double **Corr2 );


int setdim(){
  return( 1 );
}


double combine_avg( double *atmp ){
  return( atmp[0] );
}


double combine_corr( int LL , int TT , int x , double *atmp , double *ctmp ){
  double c;

  //  if( x <= LL/2 ) c = 0.0;
  //  else c = ctmp[0];
  //  else c = ctmp[0] - atmp[0]*atmp[0];
  c = ctmp[0];
  //  c = ctmp[0] - atmp[0]*atmp[0];

  return( c );
}


void write_corr( int LL , int TT , int n , double *Corr , double **Corr2 ){
  int min,max,x,y;
  double cov;

  //  min = LL/2 + 1;
  min = 0;
  max = min + TT;

  printf("AVERAGES\n");
  for(x=min;x<max;x++) printf("%d %le\n",x-min,Corr[x]);

  printf("COVARIANCES\n");
  for(x=min;x<max;x++) 
    for(y=min;y<max;y++){
      cov = ( Corr2[x][y] - Corr[x]*Corr[y] ) / n;
      printf("%d %d %le\n",x-min,y-min,cov);
    }

  return;
}



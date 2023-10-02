
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

  //  if( x <= LL/2 ) c = ctmp[0] - atmp[0]*atmp[0];
  if( x <= LL/2 ) c = ctmp[0]/(atmp[0]*atmp[0]) - 1.0;
  //  if( x <= LL/2 ) c = ctmp[0];
  else c = 0.0;

  return( c );
}


void write_corr( int LL , int TT , int n , double *Corr , double **Corr2 ){
  int max,x,y;
  double cov;

  max = LL/2 + 1;

  printf("AVERAGES\n");
  for(x=0;x<max;x++) printf("%d %le\n",x,Corr[x]);

  printf("COVARIANCES\n");
  for(x=0;x<max;x++) 
    for(y=0;y<max;y++){
      cov = ( Corr2[x][y] - Corr[x]*Corr[y] ) / n;
      printf("%d %d %le\n",x,y,cov);
    }

  return;
}



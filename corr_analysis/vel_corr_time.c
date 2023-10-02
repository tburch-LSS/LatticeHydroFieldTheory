
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int setdim();
double combine_avg( double *atmp );
double combine_corr( int LL , int TT , int x , double *atmp , double *ctmp );
void write_corr( int LL , int TT , int n , double *Corr , double **Corr2 );


int setdim(){
  return( 4 );
}


double combine_avg( double *atmp ){
  int i;
  double a;

  a = atmp[3];

  return( a );
}


double combine_corr( int LL , int TT , int x , double *atmp , double *ctmp ){
  int i,j;
  double a,c;

  //  c = ctmp[33] - atmp[3]*atmp[3];
  c = ctmp[33];

  return( c );
}


void write_corr( int LL , int TT , int n , double *Corr , double **Corr2 ){
  int maxt,t,t2;
  double cov;

  maxt = TT/2 + 1;

  printf("AVERAGES\n");
  for(t=0;t<maxt;t++) printf("%d %le\n",t,Corr[t]);

  printf("COVARIANCES\n");
  for(t=0;t<maxt;t++) 
    for(t2=0;t2<maxt;t2++){
      cov = ( Corr2[t][t2] - Corr[t]*Corr[t2] ) / n;
      printf("%d %d %le\n",t,t2,cov);
    }

  return;
}




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

  a = atmp[0];

  return( a );
}


double combine_corr( int LL , int TT , int x , double *atmp , double *ctmp ){
  int i,j;
  double a,c;

  c = ctmp[0] - atmp[0]*atmp[0];
  //  c = ctmp[33];

  return( c );
}


void write_corr( int LL , int TT , int n , double *Corr , double **Corr2 ){
  int maxt,maxx,x,t,x2,t2,i,i2;
  double cov;

  maxt = TT/2 + 1;
  maxx = LL/2 + 1;

  printf("AVERAGES\n");
  for(t=0,i=0;t<maxt;t++) for(x=0;x<maxx;x++,i++) printf("%d %d %le\n",t,x,Corr[i]);

  printf("COVARIANCES\n");
  for(t=0,i=0;t<maxt;t++) 
    for(x=0;x<maxx;x++,i++) 
      for(t2=0,i2=0;t2<maxt;t2++) 
	for(x2=0;x2<maxx;x2++,i2++){
	  cov = ( Corr2[i][i2] - Corr[i]*Corr[i2] ) / n;
	  printf("%d %d %d %d %le\n",t,x,t2,x2,cov);
	}

  return;
}



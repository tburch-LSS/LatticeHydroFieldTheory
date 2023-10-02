
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int setdim();
double combine_avg( double *atmp );
double combine_corr( int LL , int TT , int x , double *atmp , double *ctmp );
void write_corr( int LL , int TT , int n , double *Corr , double **Corr2 );


int setdim(){
  return( 10 );
}


double combine_avg( double *atmp ){
  int i;
  double a;

  //  a = atmp[5];
  a = ( 2.0*atmp[0] - atmp[4] - atmp[7] ) / 3.0;

  return( a );
}


double combine_corr( int LL , int TT , int x , double *atmp , double *ctmp ){
  int i,j;
  double a,c;

  //  c = ctmp[55];
  c = ( 4.0*ctmp[0] - 2.0*ctmp[4] - 2.0*ctmp[7] 
	- 2.0*ctmp[40] + ctmp[44] + ctmp[47] 
	- 2.0*ctmp[70] + ctmp[74] + ctmp[77] ) / 9.0;

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



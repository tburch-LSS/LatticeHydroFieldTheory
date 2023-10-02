
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int setdim();
double combine_avg( double *atmp );
double combine_corr( int LL , int TT , int x , double *atmp , double *ctmp );
void write_corr( int LL , int TT , int n , double *Corr , double **Corr2 );


int setdim(){
  return( 6 );
}


double combine_avg( double *atmp ){
  int i;
  double a;

  a = 0.0;
  //  for(i=0;i<6;i++) a += atmp[i];
  a = atmp[5];
  //  a = atmp[0] + atmp[1] + atmp[3];
  //  a = atmp[2] + atmp[4] + atmp[5];

  return( a );
}


double combine_corr( int LL , int TT , int x , double *atmp , double *ctmp ){
  int i,j;
  double c;

  c = 0.0;
  if( x <= LL/2 ) c = ctmp[55];
  //  if( x <= LL/2 ) c = ctmp[0] - atmp[0]*atmp[0];
  //  if( x <= LL/2 ) c = ctmp[0] + ctmp[11] + ctmp[33];
  //  if( x <= LL/2 ) c = ctmp[0] + ctmp[1] + ctmp[3] 
  //		    + ctmp[10] + ctmp[11] + ctmp[13] 
  //		    + ctmp[30] + ctmp[31] + ctmp[33];
  //  if( x <= LL/2 ) c = ctmp[22] + ctmp[44] + ctmp[55];
  //  if( x <= LL/2 )   for(i=0;i<6;i++) c += ctmp[10*i+i] - atmp[i]*atmp[i];
  //  if( x <= LL/2 ) for(i=0;i<6;i++) c += ctmp[10*i+i];
  //  if( x <= LL/2 ) for(i=0;i<6;i++) for(j=0;j<6;j++) c += ctmp[10*i+j];
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



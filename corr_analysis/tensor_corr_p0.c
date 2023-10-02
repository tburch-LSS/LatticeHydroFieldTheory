
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

  a = 0.0;
  //  for(i=0;i<10;i++) a += atmp[i];
  //  a = atmp[5];
  //  a = atmp[1] + atmp[2] + atmp[5];
  //  a = atmp[3] + atmp[6] + atmp[8];
  //  a = (atmp[4] - atmp[7])/sqrt(2.0);
  a = ( 2.0*atmp[7] - atmp[0] - atmp[4] ) / 3.0;
  //  a = atmp[0] + atmp[4] + atmp[7];

  return( a );
}


double combine_corr( int LL , int TT , int x , double *atmp , double *ctmp ){
  int i;
  double c;

  c = 0.0;
  if( x <= LL/2 ) c = 0.0;
  //  else c = ctmp[11] + ctmp[12] + ctmp[15] + 
  //  	 ctmp[21] + ctmp[22] + ctmp[25] + 
  //  	 ctmp[51] + ctmp[52] + ctmp[55];
  //  else for(i=0;i<10;i++) c += ctmp[10*i+i];
  //  else c = ctmp[55];
  //  else c = ctmp[11] + ctmp[22] + ctmp[55];
  //  else c = ctmp[33] + ctmp[66] + ctmp[88];
  //  else c = 0.5*( ctmp[44] - ctmp[47] 
  //		 - ctmp[74] + ctmp[77] );
  else c = ( 4.0*ctmp[77] - 2.0*ctmp[70] - 2.0*ctmp[74] 
	     - 2.0*ctmp[7] + ctmp[0] + ctmp[4] 
	     - 2.0*ctmp[47] + ctmp[40] + ctmp[44] ) / 9.0;
  //  else c = ctmp[0] + ctmp[4] + ctmp[7] 
  //	 + ctmp[40] + ctmp[44] + ctmp[47] 
  //	 + ctmp[70] + ctmp[74] + ctmp[77];
  //  else c = ctmp[0] - atmp[0]*atmp[0];

  return( c );
}


void write_corr( int LL , int TT , int n , double *Corr , double **Corr2 ){
  int min,max,x,y;
  double cov;

  min = LL/2 + 1;
  max = min + TT;

  printf("AVERAGES:\n");
  for(x=min;x<max;x++) printf("%d %le\n",x-min,Corr[x]);

  printf("COVARIANCE: %d blocks\n",n);
  for(x=min;x<max;x++) 
    for(y=min;y<max;y++){
      cov = ( Corr2[x][y] - Corr[x]*Corr[y] ) / n;
      printf("%d %d %le\n",x-min,y-min,cov);
    }

  return;
}



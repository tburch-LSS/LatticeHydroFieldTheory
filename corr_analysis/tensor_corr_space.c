
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
  //  a = atmp[1];
  //  a = atmp[4] - atmp[7];
  //  a = atmp[1] + atmp[2] + atmp[5];
  //  a = atmp[3] + atmp[6] + atmp[8];
  //  a = atmp[0] + atmp[4] - 2.0*atmp[7];
  //  a = atmp[0] + atmp[4] + atmp[7];
  //  a = atmp[0] + atmp[4] + atmp[7] + atmp[9];
  a = sqrt( atmp[1]*atmp[1] + atmp[2]*atmp[2] + atmp[5]*atmp[5] );

  return( a );
}


double combine_corr( int LL , int TT , int x , double *atmp , double *ctmp ){
  int i,j;
  double a,c;

  a = sqrt( atmp[1]*atmp[1] + atmp[2]*atmp[2] + atmp[5]*atmp[5] );

  c = 0.0;
  //  if( x <= LL/2 ) c = ctmp[11];
  //  if( x <= LL/2 ) c = ctmp[55] - atmp[5]*atmp[5];
  //  if( x <= LL/2 ) c = ctmp[11] + ctmp[22] + ctmp[55];
  if( x <= LL/2 ) c = ( ctmp[11] + ctmp[22] + ctmp[55] ) / a;
  //  if( x <= LL/2 ) c = ctmp[33] + ctmp[66] + ctmp[88];
  //  if( x <= LL/2 ) c = ctmp[44] - ctmp[47] 
  //		    - ctmp[74] + ctmp[77];
  //  if( x <= LL/2 ) c = ctmp[0] + ctmp[4] - 2.0*ctmp[7] 
  //  		    + ctmp[40] + ctmp[44] - 2.0*ctmp[47] 
  // 		    - 2.0*ctmp[70] - 2.0*ctmp[74] + 4.0*ctmp[77];
  //  if( x <= LL/2 ) c = ctmp[0] + ctmp[4] + ctmp[7] 
  //  		    + ctmp[40] + ctmp[44] + ctmp[47] 
  //  		    + ctmp[70] + ctmp[74] + ctmp[77];
  //  if( x <= LL/2 ) c = ctmp[0] + ctmp[4] + ctmp[7]  + ctmp[9] 
  //    		    + ctmp[40] + ctmp[44] + ctmp[47]  + ctmp[49] 
  //    		    + ctmp[70] + ctmp[74] + ctmp[77] + ctmp[79] 
  //    		    + ctmp[90] + ctmp[94] + ctmp[97] + ctmp[99];
  //  if( x <= LL/2 ) for(i=0;i<10;i++) c += ctmp[10*i+i];
  //  if( x <= LL/2 ) for(i=0;i<10;i++) for(j=0;j<10;j++) c += ctmp[10*i+j];
  //  if( x <= LL/2 ) for(i=0;i<10;i++) c += ctmp[10*i+i] - atmp[i]*atmp[i];
  else c = 0.0;
  //  else c = 0.0;

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



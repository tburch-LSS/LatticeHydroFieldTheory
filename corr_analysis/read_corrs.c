
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


/*  external routines for combining averages and correlators into the desired form  */
extern int setdim();
extern double combine_avg( double *tmp );
extern double combine_corr( int LL , int TT , int x , double *atmp , double *ctmp );
extern void write_corr( int LL , int TT , int n , double *Corr , double **Corr2 );


int main( int argc , char **argv ){
  int i,j,x,y,LL,TT,skip,end,max,dim,n;
  //  char filename[256];
  double atmp[10],temp,Avg,Avg2;
  double ctmp[100],*ccorr,*Corr,**Corr2;
  FILE *fp;

  argc = 6;
  sscanf(argv[2],"%d",&LL);
  sscanf(argv[3],"%d",&TT);
  sscanf(argv[4],"%d",&skip);
  sscanf(argv[5],"%d",&end);

  //  max = LL/2 + 1 + TT;
  //  max = (LL/2 + 1) * (TT/2 + 1);
  max = TT/2 + 1;
  //  max = TT;
  dim = setdim();

  ccorr = (double *)malloc(max*sizeof(double));
  Corr = (double *)malloc(max*sizeof(double));
  Corr2 = (double **)malloc(max*sizeof(double *));
  Corr2[0] = (double *)malloc(max*max*sizeof(double));
  for(x=1;x<max;x++) Corr2[x] = Corr2[0] + x*max;

  Avg = 0.0;
  Avg2 = 0.0;
  for(x=0;x<max;x++) Corr[x] = 0.0;
  for(x=0;x<max;x++) 
    for(y=0;y<max;y++) Corr2[x][y] = 0.0;

  fp = fopen(argv[1],"rb");

  n = 0;

  while( fread(&(atmp[0]),8,1,fp)  &&  n < end ){

    n += 1;

    for(i=1;i<dim;i++) fread(&(atmp[i]),8,1,fp);
    if ( n > skip ){
      temp = combine_avg( atmp );
      Avg += temp;
      Avg2 += temp*temp;
    }

    for(x=0;x<=LL/2;x++){
      for(i=0;i<dim;i++) 
	for(j=0;j<dim;j++) fread(&(ctmp[10*i+j]),8,1,fp);
      if ( n > skip ){
	ccorr[x] = combine_corr( LL , TT , x , atmp , ctmp );
	Corr[x] += ccorr[x];
      }
    }

    for(y=x;y<max;y++){
      for(i=0;i<dim;i++) 
	for(j=0;j<dim;j++) fread(&(ctmp[10*i+j]),8,1,fp);
      if ( n > skip ){
	ccorr[y] = combine_corr( LL , TT , y , atmp , ctmp );
	Corr[y] += ccorr[y];
      }
    }

    if ( n > skip ) 
      for(x=0;x<max;x++) 
	for(y=0;y<max;y++) Corr2[x][y] += ccorr[x]*ccorr[y];

  }

  fclose(fp);

  n -= skip;

  Avg /= n;
  Avg2 /= n;
  for(x=0;x<max;x++) Corr[x] /= n;
  for(x=0;x<max;x++) 
    for(y=0;y<max;y++) Corr2[x][y] /= n-1;

  fprintf(stderr,"n = %d ; Avg = %le ; Avg2 = %le\n",n,Avg,Avg2);
  write_corr( LL , TT , n , Corr , Corr2 );

  return(0);
}



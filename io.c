
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "lattice.h"
#include "macros.h"
#include "ltools.h"
#include "io.h"

#define SMALL 1.e-5


void read_lattice( ){
  register int is;
  register site *s;
  double temp;
  char *ifile;
  FILE *fp;

  ifile = strtok(params.infile," ");

  fp = fopen(ifile,"rb");

  FORALLSITES(is,s) 
    fread(&temp,8,1,fp);
    (*s).phi[0] = temp;
    fread(&temp,8,1,fp);
    (*s).phi[1] = temp;
    fread(&temp,8,1,fp);
    (*s).phi[2] = temp;
    //    fread(&temp,8,1,fp);
    //    (*s).psi = temp;
  ENDALLSITES 

  fclose(fp);

  return;
}


void dump_lattice( int traj ){
  register int is;
  register site *s;
  double temp;
  char filename[256],*ofile;
  FILE *fp;

  sprintf(filename,"%s%c%08d",params.outfile,'_',traj);
  ofile = strtok(filename," ");
  fp = fopen(ofile,"wb");

  FORALLSITES(is,s) 
    temp = (*s).phi[0];
    fwrite(&temp,8,1,fp);
    temp = (*s).phi[1];
    fwrite(&temp,8,1,fp);
    temp = (*s).phi[2];
    fwrite(&temp,8,1,fp);
    //    temp = (*s).psi;
    //    fwrite(&temp,8,1,fp);
  ENDALLSITES 

  fclose(fp);

  return;
}


void write_corrs( int traj , double *alo, Ttensor *Tmn , Ttensor_corr *TmnTmn , 
		  flow *Om, flow_corr *OmOm , phon_corr *PiPi , 
		  scal_corrs *scsc ){
  register int i,j,k,mu,nu,munu,munu2;
  double temp,temp2;
  char filename[256],*cfile;
  FILE *fp;

  /*  dump scalar fluctuations to stdout  */
  temp = scsc[0].bb - alo[23]*alo[23];
  //  printf("<b> = %le  (<b^2>-<b>^2)^1/2 = %le\n",alo[23],sqrt(temp));
  printf("<b> = %le %le\n",alo[23],sqrt(temp));
  printf("<b*F_b> = %le\n",alo[1]);
  printf("<F_bb> = %le\n",alo[8]);
  printf("<b^2*F_bb> = %le\n",alo[9]);
  temp = scsc[0].XX - alo[24]*alo[24];
  //  printf("<X> = %le %le\n",alo[24],sqrt(temp));
  temp = scsc[0].yy - alo[25]*alo[25];
  //  printf("<y> = %le %le\n",alo[25],sqrt(temp));
  temp = scsc[0].xixi - alo[26]*alo[26];
  //  printf("<xi> = %le %le\n",alo[26],sqrt(temp));
  temp = scsc[0].nn - alo[27]*alo[27];
  //  printf("<n> = %le %le\n",alo[27],sqrt(temp));
  temp = scsc[0].rhorho - alo[28]*alo[28];
  printf("<rho> = %le %le\n",alo[28],sqrt(temp));
  temp = scsc[0].pp - alo[29]*alo[29];
  printf("<p> = %le %le\n",alo[29],sqrt(temp));
  temp = scsc[0].TT - alo[30]*alo[30];
  printf("<T> = %le %le\n",alo[30],sqrt(temp));

  printf("<u_mu> = ( ");
  for( mu=0 ; mu<4 ; mu++ ) 
    printf("%le ",(*Om).u1[mu]);
  printf(")\n");

  printf("<J_mu> = ( ");
  for( mu=0 ; mu<4 ; mu++ ) 
    printf("%le ",(*Om).J1[mu]);
  printf(")\n");

  printf("<u^mu_L> = ( ");
  for( mu=0 ; mu<4 ; mu++ ) 
    printf("%le ",(*Om).uL1[(mu+3)%4]);
  printf(")\n");

  printf("<J^mu_L> = ( ");
  for( mu=0 ; mu<4 ; mu++ ) 
    printf("%le ",(*Om).JL1[(mu+3)%4]);
  printf(")\n");

  temp = 0.0;
  printf("<Omega_mu,nu>:\n");
  for( mu=0,munu=0 ; mu<4 ; mu++ ) 
    for( nu=mu ; nu<4 ; nu++,munu++ ){
      if( nu == mu ) temp += (*Om).ndir12[munu];
      //      if( nu == mu && nu < 3 ) temp -= (*Om).ndir12[munu];
      //      if( nu == mu && nu == 3 ) temp += (*Om).ndir12[munu];
      temp2 = OmOm[0].ndir1234[munu][munu] - (*Om).ndir12[munu] * (*Om).ndir12[munu];
      printf("%1d %1d %le %le\n",mu,nu,(*Om).ndir12[munu],sqrt(temp2));
    }
  if ( fabs(3.0-temp) > SMALL ) printf("<Tr(Omega)> = %le !!!\n",temp);

  temp = 0.0;
  printf("<T_mu,nu>:\n");
  for( mu=0,munu=0 ; mu<4 ; mu++ ) 
    for( nu=mu ; nu<4 ; nu++,munu++ ){
      if( nu == mu ) temp += (*Tmn).dir12[munu];
      //      if( nu == mu && nu < 3 ) temp -= (*Tmn).dir12[munu];
      //      else if( nu == mu && nu == 3 ) temp += (*Tmn).dir12[munu];
      temp2 = TmnTmn[0].dir1234[munu][munu] - (*Tmn).dir12[munu] * (*Tmn).dir12[munu];
      printf("%1d %1d %le %le\n",mu,nu,(*Tmn).dir12[munu],sqrt(temp2));
    }
  if ( fabs(temp) > SMALL ) printf("<T_mu,mu> = %le !!!\n",temp);

  //  j = params.Nx/2 + 1 + params.Nt;
  //  j = (params.Nx/2 + 1) * (params.Nt/2 + 1);
  j = params.Nt/2 + 1;
  k = j + params.Nt;

  sprintf(filename,"%s%c%c%c%c%c",params.corrfile,'_','T','T','_','t');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  for( munu=0 ; munu<10 ; munu++ ) fwrite(&(Tmn[0].dir12[munu]),8,1,fp);
  for( i=0 ; i<j ; i++ ) 
    for( munu=0 ; munu<10 ; munu++ ) 
      for( munu2=0 ; munu2<10 ; munu2++ ) fwrite(&(TmnTmn[i].dir1234[munu][munu2]),8,1,fp);
  fclose(fp);
  sprintf(filename,"%s%c%c%c%c%c%c",params.corrfile,'_','T','T','_','p','0');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  for( munu=0 ; munu<10 ; munu++ ) fwrite(&(Tmn[0].dir12[munu]),8,1,fp);
  for( i=j ; i<k ; i++ ) 
    for( munu=0 ; munu<10 ; munu++ ) 
      for( munu2=0 ; munu2<10 ; munu2++ ) fwrite(&(TmnTmn[i].dir1234[munu][munu2]),8,1,fp);
  fclose(fp);

  sprintf(filename,"%s%c%c%c%c%c",params.corrfile,'_','O','O','_','t');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  for( munu=0 ; munu<10 ; munu++ ) fwrite(&(Om[0].ndir12[munu]),8,1,fp);
  for( i=0 ; i<j ; i++ ) 
    for( munu=0 ; munu<10 ; munu++ ) 
      for( munu2=0 ; munu2<10 ; munu2++ ) fwrite(&(OmOm[i].ndir1234[munu][munu2]),8,1,fp);
  fclose(fp);
  sprintf(filename,"%s%c%c%c%c%c%c",params.corrfile,'_','O','O','_','p','0');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  for( munu=0 ; munu<10 ; munu++ ) fwrite(&(Om[0].ndir12[munu]),8,1,fp);
  for( i=j ; i<k ; i++ ) 
    for( munu=0 ; munu<10 ; munu++ ) 
      for( munu2=0 ; munu2<10 ; munu2++ ) fwrite(&(OmOm[i].ndir1234[munu][munu2]),8,1,fp);
  fclose(fp);

  sprintf(filename,"%s%c%c%c%c%c",params.corrfile,'_','u','u','_','t');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  for( mu=0 ; mu<4 ; mu++ ) fwrite(&(Om[0].u1[mu]),8,1,fp);
  for( i=0 ; i<j ; i++ ) 
    for( mu=0 ; mu<4 ; mu++ ) 
      for( nu=0 ; nu<4 ; nu++ ) fwrite(&(OmOm[i].uu12[mu][nu]),8,1,fp);
  fclose(fp);
  sprintf(filename,"%s%c%c%c%c%c%c",params.corrfile,'_','u','u','_','p','0');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  for( mu=0 ; mu<4 ; mu++ ) fwrite(&(Om[0].u1[mu]),8,1,fp);
  for( i=j ; i<k ; i++ ) 
    for( mu=0 ; mu<4 ; mu++ ) 
      for( nu=0 ; nu<4 ; nu++ ) fwrite(&(OmOm[i].uu12[mu][nu]),8,1,fp);
  fclose(fp);

  sprintf(filename,"%s%c%c%c%c%c",params.corrfile,'_','J','J','_','t');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  for( mu=0 ; mu<4 ; mu++ ) fwrite(&(Om[0].J1[mu]),8,1,fp);
  for( i=0 ; i<j ; i++ ) 
    for( mu=0 ; mu<4 ; mu++ ) 
      for( nu=0 ; nu<4 ; nu++ ) fwrite(&(OmOm[i].JJ12[mu][nu]),8,1,fp);
  fclose(fp);
  sprintf(filename,"%s%c%c%c%c%c%c",params.corrfile,'_','J','J','_','p','0');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  for( mu=0 ; mu<4 ; mu++ ) fwrite(&(Om[0].J1[mu]),8,1,fp);
  for( i=j ; i<k ; i++ ) 
    for( mu=0 ; mu<4 ; mu++ ) 
      for( nu=0 ; nu<4 ; nu++ ) fwrite(&(OmOm[i].JJ12[mu][nu]),8,1,fp);
  fclose(fp);

  sprintf(filename,"%s%c%c%c%c%c%c",params.corrfile,'_','u','u','L','_','t');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  for( mu=0 ; mu<4 ; mu++ ) fwrite(&(Om[0].uL1[mu]),8,1,fp);
  for( i=0 ; i<j ; i++ ) 
    for( mu=0 ; mu<4 ; mu++ ) 
      for( nu=0 ; nu<4 ; nu++ ) fwrite(&(OmOm[i].uuL12[mu][nu]),8,1,fp);
  fclose(fp);
  sprintf(filename,"%s%c%c%c%c%c%c%c",params.corrfile,'_','u','u','L','_','p','0');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  for( mu=0 ; mu<4 ; mu++ ) fwrite(&(Om[0].uL1[mu]),8,1,fp);
  for( i=j ; i<k ; i++ ) 
    for( mu=0 ; mu<4 ; mu++ ) 
      for( nu=0 ; nu<4 ; nu++ ) fwrite(&(OmOm[i].uuL12[mu][nu]),8,1,fp);
  fclose(fp);

  sprintf(filename,"%s%c%c%c%c%c%c",params.corrfile,'_','J','J','L','_','t');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  for( mu=0 ; mu<4 ; mu++ ) fwrite(&(Om[0].JL1[mu]),8,1,fp);
  for( i=0 ; i<j ; i++ ) 
    for( mu=0 ; mu<4 ; mu++ ) 
      for( nu=0 ; nu<4 ; nu++ ) fwrite(&(OmOm[i].JJL12[mu][nu]),8,1,fp);
  fclose(fp);
  sprintf(filename,"%s%c%c%c%c%c%c%c",params.corrfile,'_','J','J','L','_','p','0');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  for( mu=0 ; mu<4 ; mu++ ) fwrite(&(Om[0].JL1[mu]),8,1,fp);
  for( i=j ; i<k ; i++ ) 
    for( mu=0 ; mu<4 ; mu++ ) 
      for( nu=0 ; nu<4 ; nu++ ) fwrite(&(OmOm[i].JJL12[mu][nu]),8,1,fp);
  fclose(fp);

  sprintf(filename,"%s%c%c%c%c%c",params.corrfile,'_','v','v','_','t');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  for( munu=0 ; munu<6 ; munu++ ) fwrite(&(Om[0].nvort[munu]),8,1,fp);
  for( i=0 ; i<j ; i++ ) 
    for( munu=0 ; munu<6 ; munu++ ) 
      for( munu2=0 ; munu2<6 ; munu2++ ) fwrite(&(OmOm[i].nvort12[munu][munu2]),8,1,fp);
  fclose(fp);
  sprintf(filename,"%s%c%c%c%c%c%c",params.corrfile,'_','v','v','_','p','0');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  for( munu=0 ; munu<6 ; munu++ ) fwrite(&(Om[0].nvort[munu]),8,1,fp);
  for( i=j ; i<k ; i++ ) 
    for( munu=0 ; munu<6 ; munu++ ) 
      for( munu2=0 ; munu2<6 ; munu2++ ) fwrite(&(OmOm[i].nvort12[munu][munu2]),8,1,fp);
  fclose(fp);


  j = (params.Nx/2 + 1) * (params.Nt/2 + 1);
  k = 3*j + params.Nt;

  sprintf(filename,"%s%c%c%c%c%c%c%c%c",params.corrfile,'_','p','i','p','i','_','x','t');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  for( mu=0 ; mu<3 ; mu++ ) fwrite(&(Om[0].pi[mu]),8,1,fp);
  for( i=0 ; i<j ; i++ ) 
    for( mu=0 ; mu<3 ; mu++ ) 
      for( nu=0 ; nu<3 ; nu++ ) fwrite(&(PiPi[i].dir12[mu][nu]),8,1,fp);
  fclose(fp);
  sprintf(filename,"%s%c%c%c%c%c%c%c%c",params.corrfile,'_','p','i','p','i','_','y','t');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  for( mu=0 ; mu<3 ; mu++ ) fwrite(&(Om[0].pi[mu]),8,1,fp);
  for( i=j ; i<2*j ; i++ ) 
    for( mu=0 ; mu<3 ; mu++ ) 
      for( nu=0 ; nu<3 ; nu++ ) fwrite(&(PiPi[i].dir12[mu][nu]),8,1,fp);
  fclose(fp);
  sprintf(filename,"%s%c%c%c%c%c%c%c%c",params.corrfile,'_','p','i','p','i','_','z','t');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  for( mu=0 ; mu<3 ; mu++ ) fwrite(&(Om[0].pi[mu]),8,1,fp);
  for( i=2*j ; i<3*j ; i++ ) 
    for( mu=0 ; mu<3 ; mu++ ) 
      for( nu=0 ; nu<3 ; nu++ ) fwrite(&(PiPi[i].dir12[mu][nu]),8,1,fp);
  fclose(fp);
  sprintf(filename,"%s%c%c%c%c%c%c%c%c",params.corrfile,'_','p','i','p','i','_','p','0');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  for( mu=0 ; mu<3 ; mu++ ) fwrite(&(Om[0].pi[mu]),8,1,fp);
  for( i=3*j ; i<k ; i++ ) 
    for( mu=0 ; mu<3 ; mu++ ) 
      for( nu=0 ; nu<3 ; nu++ ) fwrite(&(PiPi[i].dir12[mu][nu]),8,1,fp);
  fclose(fp);

  /*
  sprintf(filename,"%s%c%c%c%c%c",params.corrfile,'_','p','i','p','i');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  fwrite(&(Om[0].pi[3]),8,1,fp);
  for( i=0 ; i<j ; i++ ) 
    fwrite(&(PiPi[i].dir12[3][3]),8,1,fp);
  fclose(fp);
  */

  sprintf(filename,"%s%c%c%c%c%c%c",params.corrfile,'_','b','b','_','x','t');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  fwrite(&(alo[23]),8,1,fp);
  for( i=0 ; i<j ; i++ ) fwrite(&(scsc[i].bb),8,1,fp);
  fclose(fp);
  sprintf(filename,"%s%c%c%c%c%c%c",params.corrfile,'_','b','b','_','y','t');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  fwrite(&(alo[23]),8,1,fp);
  for( i=j ; i<2*j ; i++ ) fwrite(&(scsc[i].bb),8,1,fp);
  fclose(fp);
  sprintf(filename,"%s%c%c%c%c%c%c",params.corrfile,'_','b','b','_','z','t');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  fwrite(&(alo[23]),8,1,fp);
  for( i=2*j ; i<3*j ; i++ ) fwrite(&(scsc[i].bb),8,1,fp);
  fclose(fp);
  sprintf(filename,"%s%c%c%c%c%c%c",params.corrfile,'_','b','b','_','p','0');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  fwrite(&(alo[23]),8,1,fp);
  for( i=3*j ; i<k ; i++ ) fwrite(&(scsc[i].bb),8,1,fp);
  fclose(fp);

  sprintf(filename,"%s%c%c%c%c%c%c",params.corrfile,'_','e','e','_','x','t');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  fwrite(&(alo[28]),8,1,fp);
  for( i=0 ; i<j ; i++ ) fwrite(&(scsc[i].rhorho),8,1,fp);
  fclose(fp);
  sprintf(filename,"%s%c%c%c%c%c%c",params.corrfile,'_','e','e','_','y','t');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  fwrite(&(alo[28]),8,1,fp);
  for( i=j ; i<2*j ; i++ ) fwrite(&(scsc[i].rhorho),8,1,fp);
  fclose(fp);
  sprintf(filename,"%s%c%c%c%c%c%c",params.corrfile,'_','e','e','_','z','t');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  fwrite(&(alo[28]),8,1,fp);
  for( i=2*j ; i<3*j ; i++ ) fwrite(&(scsc[i].rhorho),8,1,fp);
  fclose(fp);
  sprintf(filename,"%s%c%c%c%c%c%c",params.corrfile,'_','e','e','_','p','0');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  fwrite(&(alo[28]),8,1,fp);
  for( i=3*j ; i<k ; i++ ) fwrite(&(scsc[i].rhorho),8,1,fp);
  fclose(fp);

  /*
  sprintf(filename,"%s%c%c%c",params.corrfile,'_','t','t');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  fwrite(&(alo[30]),8,1,fp);
  for( i=0 ; i<j ; i++ ) fwrite(&(scsc[i].TT),8,1,fp);
  fclose(fp);
  */

  /*
  sprintf(filename,"%s%c%c%c",params.corrfile,'_','y','y');
  cfile = strtok(filename," ");
  fp = fopen(cfile,"ab");
  fwrite(&(alo[25]),8,1,fp);
  for( i=0 ; i<j ; i++ ) fwrite(&(scsc[i].yy),8,1,fp);
  fclose(fp);
  */

  return;
}



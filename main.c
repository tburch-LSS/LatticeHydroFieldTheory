/* 
   lattice simulation code of EFT for normal/super-fluid hydro/thermodynamics: 
   needed input: EoS; specifically,  p - T*dp/dT  as a function of (s,mu,xi^2: 
   entropy density, chemical potential, normal/super-fluid squared relative motion); 
   after hydro/thermo matching, this is the Lagrangian density F(b,X,y), where b=s, 
   X=xi^2+y^2, and y=mu 

   normal fluid fields (comoving coordinates): phi^I (--> b,y) 
   conserved charge (possibly superfluid) U(1) phase: psi (--> X,y) 

   [actually, the phi^I and psi fields herein differ from those in the refs below 
   by a simple shift: phi^I --> phi^I - x^I ; psi --> psi - mu*t] 

   see the following articles: 

   D.T.Son, arXiv:hep-ph/0204199 (relativistic superfluids) 
   S.Endlich et al., arXiv:1011.6396 (relativistic perfect fluids) 
   S.Dubovsky et al., arXiv:1107.0731 (relativistic fluids w/conserved charges) 
   A.Nicolis, arXiv:1108.2513 (relativistic two-fluid approach; this version is 
   coded herein) ... 

   ... and refs therein for further details 

   see G.Torrieri, arXiv:1112.4086 for a perturbative calculation of quantum 
   viscosity of a perfect fluid 
   see V.Cirigliano etal., arXiv:1102.5379 for possible extensions with "solids" 
   (NS crust, supersolids) 
   see P.Huovinen etal., arXiv:0912.2541 & 1202.3104 for example parametrizations 
   of (lattice) QCD EoS 

   other related articles: 
   P.Romatschke, arXiv:0902.3663 ; P.Kovtun, arXiv:1205.5040 ; 
   T.Kalaydzhyan, arXiv:1208.0012 

   initial coding and testing phase: 2011-2012 T.Burch (special thanks and apologies 
   to his family for the late nights spent catching up with his old friend, astrophysics) 
   (needed files: main.c, lattice.h, macros.h, setup.c, ltools.c, io.c, local_meas.c, 
   glob_meas.c, action.c, update.c; also ranlxd.* for Luscher's random number generator) 

   03/2013 (T.B.): OpenMP parallelization added, along with correlator output and more 
   observables (e.g., vorticity) 

   04/2013 (T.B.): correlator analysis code created (see read_corrs.c, 
   scalar_corr_space.c, etc. in corr_analysis/ subdirectory) 
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAIN 

#include "lattice.h"
#include "macros.h"
#include "ranlxd.h"
#include "ltools.h"
#include "setup.h"
#include "io.h"
#include "local_meas.h"
#include "action.h"
#include "update.h"
#include "glob_meas.h"


/*  main program  */
int main( ){
  register int is,n,dmet,met,ns,traj,t,nd,mu,nu,munu,js;
  register site *s;
  int dr[4];
  Ttensor Tmn;
  flow Om;
  double *AvgLocObs,dJ,temp,v2,gamma,uu_L[4];
  Ttensor_corr *TmnTmn;
  flow_corr *OmOm;
  phon_corr *PiPi;
  scal_corrs *scsc;

  /*  read input data from stdin  */
  scanf("%*s%d",&(params.Nx));
  scanf("%*s%d",&(params.Ny));
  scanf("%*s%d",&(params.Nz));
  scanf("%*s%d",&(params.Nt));
  params.Vol2 = params.Nx * params.Ny;
  params.Vol3 = params.Vol2 * params.Nz;
  params.Vol4 = params.Vol3 * params.Nt;
  Nsite = params.Vol4;
  Nspace = params.Vol3;
  scanf("%*s%d",&(params.iseed));
  scanf("%*s%d",&(params.upd_phi));
  scanf("%*s%le",&(params.S1));
  scanf("%*s%le",&(params.S2));
  scanf("%*s%le",&(params.S3));
  scanf("%*s%le",&(params.chem_pot));
  scanf("%*s%le",&(params.dt));
  scanf("%*s%d",&(params.Met_step));
  scanf("%*s%d",&(params.N_therm));
  scanf("%*s%d",&(params.N_meas));
  scanf("%*s%d",&(params.N_step));
  scanf("%*s%d",&(params.N_dump));
  scanf("%*s%d",&(params.start));  /* 0=cold, 1=file, 2=hot */
  scanf("%*s%d",&(params.out));  /* 0=forget_lattice, 1=dump_to_file_after_meas */
  scanf("%*s%s",&(params.infile));
  scanf("%*s%s",&(params.outfile));
  scanf("%*s%s",&(params.corrfile));

  /*  initialize random # generator  */
  rlxd_init( 2 , params.iseed );
  //  srand48( params.iseed );

  /*  set initial field parameters based upon the action  */
  action_init();

  /*  allocate memory for lattice sites  */
  lattice = (site *) malloc( Nsite * sizeof(site) );
  n = 2 * ( Nspace + 3 * params.Nt * params.Vol2 );
  if( params.bc != PERIODIC ) 
    latt_bound = (site *) malloc( n * sizeof(site) );

  /*  allocate memory for averaged, local observables  */
  AvgLocObs = (double *) malloc( NUM_LOCAL_OBS * sizeof(double) );

  /*  allocate memory for tensor correlators  */
  //  n = params.Nx/2 + 1 + params.Nt;
  //  n = (params.Nx/2 + 1) * (params.Nt/2 + 1);
  n = params.Nt/2 + 1 + params.Nt;
  TmnTmn = (Ttensor_corr *) malloc( n * sizeof(Ttensor_corr) );
  OmOm = (flow_corr *) malloc( n * sizeof(flow_corr) );

  /*  allocate memory for scalar correlators  */
  n = 3 * (params.Nx/2 + 1) * (params.Nt/2 + 1) + params.Nt;
  PiPi = (phon_corr *) malloc( n * sizeof(phon_corr) );
  scsc = (scal_corrs *) malloc( n * sizeof(scal_corrs) );

  /*  setup lattice and initialize lattice fields  */
  setup();
  //  check_setup(); fflush(stdout);
  init_fields();
  Saction = action();
  random_conj_mom();
  Hamilto = hamiltonian( Saction );

  /*  thermalization updates  */
  //  fprintf(stderr,"thermalization updates:\n");
  traj = 0;
  met = 0;
  for( n=0 ; n<params.N_therm ; n++, traj++ ){
    //    dmet = Met_update();
    dmet = HMC_update();
    met += dmet;
    //    fprintf(stderr," traj=%08d met=%d met_tot=%d S=%le H=%le\n", 
    //	    traj,dmet,met,Saction,Hamilto);
  }

  /*  measurement updates  */
  //  fprintf(stderr,"measurement updates:\n");
  met = 0; nd = 0;
  for( n=0 ; n<params.N_meas ; n++ ){

    /*  N_step trajectories between measurements  */
    for( ns=0 ; ns<params.N_step ; ns++, traj++ ){
      //      dmet = Met_update();
      dmet = HMC_update();
      met += dmet;
      //      fprintf(stderr," traj=%08d met=%d met_tot=%d S=%le H=%le\n", 
      //	      traj,dmet,met,Saction,Hamilto);
    }

    printf(" traj= %d  met= %d  S= %le  H= %le\n",traj,met,Saction,Hamilto);

    /*  averaged, local observables: e.g., <b>, <dF/db>, etc.  */
    dS_dbdXdy( AvgLocObs );

    is = 0;
    s = &(lattice[is]);
    dJ = divJ(s);
    printf("div.J(x=%d) = %le\n",is,dJ);


    eucl_to_lrnz_vel( s , uu_L );
    gamma = gammaL( s );
    v2 = 1.0 - 1.0 / (gamma*gamma);
    printf("v^2(x=%d) = %le  ;  gamma(x) = %le\n",is,v2,gamma);
    printf("u_mu_E(x) = ( %le , %le , %le , %le )\n",(*s).uu[0],(*s).uu[1],(*s).uu[2],(*s).uu[3]);
    printf("u^mu_L(x) = ( %le , %le , %le , %le )\n",uu_L[3],uu_L[0],uu_L[1],uu_L[2]);


    /*
    Tmn_match( s , &Tmn );
    printf("T_mu,nu(x=%d):\n",is);
    temp = 0.0;
    for( mu=0,munu=0 ; mu<4 ; mu++ ) 
      for( nu=mu ; nu<4 ; nu++,munu++ ){
	if( nu == mu ) temp += Tmn.dir12[munu];
	printf("%1d %1d %le\n",mu,nu,Tmn.dir12[munu]);
      }
    if ( fabs(temp) > 1.e-5 ) printf("T_mu,mu(x=%d) = %le !!!\n",is,temp);
    //    printf("T_mu,mu(x=%d) = %le\n",is,temp);
    dr[0] = params.Nx/2; dr[1] = params.Ny/2; dr[2] = params.Nz/2; dr[3] = 0;
    js = new_index_from_dr(is,dr);
    s = &(lattice[js]);
    Tmn_match( s , &Tmn );
    printf("T_mu,nu(x=%d):\n",js);
    temp = 0.0;
    for( mu=0,munu=0 ; mu<4 ; mu++ ) 
      for( nu=mu ; nu<4 ; nu++,munu++ ){
	if( nu == mu ) temp += Tmn.dir12[munu];
	printf("%1d %1d %le\n",mu,nu,Tmn.dir12[munu]);
      }
    if ( fabs(temp) > 1.e-5 ) printf("T_mu,mu(x=%d) = %le !!!\n",is,temp);
    */

    //    printf("u_mu(x) = ( %le , %le , %le , %le )\n",(*s).uu[0],(*s).uu[1],(*s).uu[2],(*s).uu[3]);
    //    temp = (*s).uu[0] * (*s).uu[0] + (*s).uu[1] * (*s).uu[1] + (*s).uu[2] * (*s).uu[2] + (*s).uu[3] * (*s).uu[3];
    //    printf("u_mu(x).u_mu(x) = %le\n",temp);

    /*  stress-energy and flow (tensor,vector) averages and their two-point correlations 
	(including scalars)  */
    Correlators( &Tmn , TmnTmn , &Om , OmOm , PiPi , scsc );

    /*  write out observables  */
    write_corrs( traj , AvgLocObs , &Tmn , TmnTmn , &Om , OmOm , PiPi , scsc );

    fflush(stdout);

    /*  save present {phi^I(x),psi(x)} config?  */
    nd += 1;
    if( params.out && nd == params.N_dump ){
      dump_lattice( traj );
      nd = 0;
    }

  }

  /*  free allocated memory  */
  free(AvgLocObs);
  free(OmOm);
  if( params.bc != PERIODIC ) free(latt_bound);
  free(lattice);

  return(0);
}


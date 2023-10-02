

typedef struct {
  int Nx,Ny,Nz,Nt;
  int Vol2,Vol3,Vol4;
  int iseed;
  int upd_phi;
  double S1,S2,S3;
  double chem_pot;
  double dt;
  int Met_step,N_therm,N_meas,N_step,N_dump;
  int start,out;
  char infile[256];
  char outfile[256];
  char corrfile[256];
  double b0,y0;
  int bc;
} latt_params ;


typedef struct {

  //  int xyzt[4];
  int index;
  char *neib[8];

#ifdef DUAL_LATTICE
  char *dneib[16];
  char *dneibBC[16];
#endif

  /*  U(1) link phases ?  */
  //  double phase[4] ;

  /*  the phi^I and psi fields herein differ from those in the refs by a 
      simple shift: phi^I --> pi^I = phi^I - x^I ; psi --> pi^3 = psi - mu*t ; 
      the derivatives take these shifts into account, however  */

  /*  superfluid phase and 1st derivs  */
  double psi , dpsi[4] ;
  double psi_old ;
  /*  conjugate momentum variable (for MD update)  */
  double psiCM ;

  /*  normal fluid fields (comoving coords.) and 1st derivs  */
  double phi[3] , dphi[3][4] ;
  double AImu[3][4] ;
  //  double Binv[3][3] ;
  double phi_old[3] ;
  /*  conjugate momentum variables  */
  double phiCM[3] ;

  /*  normal-fluid velocity: u_mu(x)  */
  double uu[4] ;

  /*  local scalar invariants  */
  double b , X , y ;

} site ;


typedef struct {
  /*  pi^I = phi^I - x^I ; pi^3 = psi - mu*t ; 
      phonon field correlator: <pi^a(x) pi^b(x')>  */
  double dir12[4][4];
} phon_corr;

typedef struct {
  /*  stress-energy tensor: T_mu,nu(x)  */
  double dir12[10];
} Ttensor;

typedef struct {
  /*  stress-energy tensor correlator: <T_mu,nu(x) T_mu',nu'(x')>  */
  double dir1234[10][10];
} Ttensor_corr;

typedef struct {
  /*  normal-fluid flow velocities: <u_mu(x)> , <J_mu(x)>  */
  double u1[4],J1[4],uL1[4],JL1[4];
  /*  normal-fluid flow tensor: Omega_mu,nu(x)  */
  double ndir12[10];
  /*  normal/superfluid relative motion: xi_mu(x)  */
  double sdir[4];
  /*  normal-fluid vorticity  */
  double nvort[6];
  /*  superfluid vorticity  */
  double svort[6];
  /*  phonon fields pi^a and derivs d_mu pi^a  */
  double pi[4];
  //  double dpi[4][4];
} flow;

typedef struct {
  /*  normal-fluid flow velocity correlators: <u_mu(x) u_nu(x')> , <J_mu(x) J_nu(x')>  */
  double uu12[4][4],JJ12[4][4],uuL12[4][4],JJL12[4][4];
  /*  normal-fluid flow tensor correlator: <Omega_mu,nu(x) Omega_mu',nu'(x')>  */
  double ndir1234[10][10];
  /*  relative flow correlator: <xi_mu(x) xi_nu(x')>  */
  double sdir12[4][4];
  /*  normal-fluid vorticity  */
  double nvort12[6][6];
  /*  superfluid vorticity  */
  double svort12[6][6];
  /*  phonon field correlator: <pi^a(x) pi^b(x')>  */
  //  phon_corr pipi;
  //  phon_corr dpidpi[10];
  /*  scalar correlators: <b(x) b(x')> , etc.  */
  //  double bb,XX,yy,xixi,nn,rhorho,pp,TT;
} flow_corr;

typedef struct {
  double bb,XX,yy,xixi,nn,rhorho,pp,TT;
} scal_corrs;


#define NUM_LOCAL_OBS 31

#ifdef MAIN
#define EXTERN 
#else
#define EXTERN extern
#endif

EXTERN int Nsite , Nspace ;
EXTERN double Saction , Hamilto ;
EXTERN latt_params params ;

#ifdef DUAL_LATTICE
EXTERN int dsign[4][16] ;
#endif

EXTERN site *lattice , *latt_bound ;



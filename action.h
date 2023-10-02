
void action_init( );
double dLagrangian( site *s );
void dF_dbXy( double *dF, site *s );
void d2F_dbXy2( double *d2F, site *s );
void EoS_match( site *s , double *n_rho_p_T );
void Tmn_match( site *s , Ttensor *Tmn );
double dS_ddmuphi( site *s, int i, int mu );
double dS_ddmupsi( site *s, int mu );
void init_bound( );


void field_derivs( );

double deriv( field_offset fo , site *s , int mu );
double derivBC( field_offset fo , site *s , int mu );
double avgBC( field_offset fo , site *s );
double divJ( site *s );
void vorticity( field_offset velFO , site *s , double *vort );

double gammaL( site *s );
void eucl_to_lrnz_vel( site *s , double *uu_L );
void lrnz_flow_match( site *s , flow *Om );

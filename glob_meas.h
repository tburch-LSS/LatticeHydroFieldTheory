
double action( );
double hamiltonian( double S );
void dS_dbdXdy( double *output );
void Correlators( Ttensor *Tmn , Ttensor_corr *TmnTmn , 
		  flow *Om , flow_corr *OmOm , 
		  phon_corr *PiPi , scal_corrs *scsc );

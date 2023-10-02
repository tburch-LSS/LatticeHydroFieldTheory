void read_lattice( );
void dump_lattice( int traj );
void write_corrs( int traj , double *AvgLocObs , 
		  Ttensor *Tmn , Ttensor_corr *TmnTmn , 
		  flow *Om , flow_corr *OmOm , 
		  phon_corr *PiPi , scal_corrs *scsc );

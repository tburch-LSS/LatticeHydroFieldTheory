#!/bin/csh

#set exec = "IG_DL_OMP.x"
#set exec = "IG_DL_OMP_v2.x"
#set exec = "IG_DL_OMP_v3.x"
#set exec = "IG_DL_OMP_v3_tcorr.x"
#set exec = "IG_DL_OMP_v3_xtcorr.x"
#set exec = "IG_DL_OMP_v4_xtcorr.x"
#set exec = "IG_DL_OMP_v5L_xtcorr.x"
#set exec = "IG_DL_OMP_v6E_xtcorr.x"
#set exec = "IG_DL_OMP_v7E_xtcorr.x"
#set exec = "IG_DL_OMP_v8E_xtcorr.x"
#set exec = "IG_DL_OMP_v9E_xtcorr.x"
#set exec = "IG_DL_OMP_v10_xtcorr.x"
set exec = "IG_DL_OMP_v11.x"
#set exec = "RDFG_DL_OMP_v8E_xtcorr.x"
#set exec = "IGwCP_DL_OMP_v2.x"
#set exec = "blah.x"
#set exec = "test.x"
#set exec = "testDL.x"
#set exec = $1

#gcc -O3 -c *.c
#gcc -O3 -DDUAL_LATTICE -c *.c
gcc -fopenmp -O3 -DDUAL_LATTICE -c *.c

gcc *.o -fopenmp -lm -o $exec


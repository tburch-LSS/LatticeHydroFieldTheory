
#define XUP 0
#define YUP 1
#define ZUP 2
#define TUP 3
#define TDOWN 4
#define ZDOWN 5
#define YDOWN 6
#define XDOWN 7
#define OPP_DIR(dir) (7-dir)

#define PERIODIC 0
#define TIME_NeuBC 1
#define TIME_MixBC 2
#define ALL_NeuBC 3

typedef int field_offset;
#define F_OFFSET(a) \
  ((field_offset)(((char *)&(lattice[0]. a ))-((char *)&(lattice[0])))) 
#define F_PT(site,fo)  ((char *)(site) + (fo)) 

#define FORALLSITES(i,s)  for(i=0;i<Nsite;i++){s=&(lattice[i]); 
#define ENDALLSITES  } 

#define FORALLSPACE(i,s,tt)  for(i=tt*Nspace,s=&(lattice[i]);i<(tt+1)*Nspace;i++,s++) 

#define FORALLXBOUND(i,s,xx)  for(i=xx,s=&(lattice[i]);i<Nsite;i+=params.Nx,s+=params.Nx) 

#define FORALLYBOUND(i,s,j,yy)  for(i=yy*params.Nx,s=&(lattice[i]);i<Nsite;i+=params.Vol2,s+=params.Vol2)for(j=0;j<params.Nx;j++,i++,s++) 

#define FORALLZBOUND(i,s,j,zz)  for(i=zz*params.Vol2,s=&(lattice[i]);i<Nsite;i+=Nspace,s+=Nspace)for(j=0;j<params.Vol2;j++,i++,s++) 


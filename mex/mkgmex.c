#include <mex.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "mkgmex.h"
void maxmin(double *p, int N, double *max, double *min){
    double a,b;
    int i;
    a=p[0];
    b=p[0];
    for(i=1; i<N; i++){
	if(p[i]>a) a=p[i];
	if(p[i]<b) b=p[i];
    }
    *max=a;
    *min=b;
}
dsp * spnew(long nx, long ny, long nzmax){
    /*a nx*ny dsp matrix with nmax max elements.*/
    dsp *sp;
    sp = calloc(sizeof(dsp),1);
    sp->p=malloc((ny+1)*sizeof(long));
    sp->i=malloc(nzmax*sizeof(long));
    sp->x=malloc(nzmax*sizeof(double));
    sp->m=nx;
    sp->n=ny;
    sp->nzmax=nzmax;
    sp->nz=-1;
    sp->shmkey=0;
    return sp;
}
void loc_free_map(LOC_T *loc){
    if(loc->map && loc->map->loc==loc){
	free(loc->map->p);
	loc->map->p=NULL;
	loc->map->loc=NULL;
	free(loc->map);
	loc->map=NULL;
    }
}
void loc_create_map(LOC_T *loc){
#define NPAD 0 /*We don't need pad. This makes wrapping possible*/
    if(loc->map){
	if(loc->map->loc==loc)/*already have a map that is valid.*/
	    return;
	else
	    loc->map=NULL;
    }
    loc->map = calloc(1,sizeof(LOCMAP_T));
    loc->map->loc = loc;
    double xmin,xmax,ymin,ymax;
    maxmin(loc->locx, loc->nloc, &xmax, &xmin);
    maxmin(loc->locy, loc->nloc, &ymax, &ymin);
    int map_nx, map_ny;
    map_nx=(int) ceil((xmax-xmin)/loc->dx)+2+NPAD*2;
    map_ny=(int) ceil((ymax-ymin)/loc->dx)+2+NPAD*2;
    loc->map->p=calloc(map_nx*map_ny,sizeof(long));
    loc->map->nx=map_nx;
    loc->map->ny=map_ny;
    loc->map->ox=xmin;
    loc->map->oy=ymin;
    loc->map->npad=NPAD;
    int iloc;
    long (*map)[map_nx]=(long(*)[map_nx])loc->map->p;
    const double *locx,*locy;
    locx=loc->locx;
    locy=loc->locy;
    double dx_in1=1./loc->dx;
    int ix,iy;
    for(iloc=0; iloc<loc->nloc; iloc++){
	ix=(int)round((locx[iloc]-xmin)*dx_in1)+NPAD;
	iy=(int)round((locy[iloc]-ymin)*dx_in1)+NPAD;
	map[iy][ix]=iloc+1;/*start from 1.*/
    }
}
static void __attribute__((constructor)) init(void){
    fprintf(stderr, "%s: Compiled on %s %s by %s ", 
	    __BASE_FILE__, __DATE__, __TIME__, __VERSION__);
#ifdef __OPTIMIZE__
    fprintf(stderr, "with optimization.\n");
#else
    fprintf(stderr, "without optimization!!!\n");
#endif
}
dsp *spcat(const dsp *A, const dsp *B, int type){
    dsp *C=NULL;
    if(type==0){
	error("Not implemented\n");
	/*
	  |A|
	  |B|
	 */
    }else if(type==1){
	/*|AB|*/
	if(A->m != B->m){
	    error("Dsp matrix doesn't match\n");
	}
	const long nzmax=A->nzmax+B->nzmax;
	C=spnew(A->m, A->n+B->n, nzmax);
	memcpy(C->p, A->p, A->n*sizeof(long));
	memcpy(C->i, A->i, A->nzmax*sizeof(long));
	memcpy(C->x, A->x, A->nzmax*sizeof(double));
	memcpy(C->i+A->nzmax, B->i, B->nzmax*sizeof(long));
	memcpy(C->x+A->nzmax, B->x, B->nzmax*sizeof(long));
	const long Anzmax=A->nzmax;
	long i;
	for(i=0; i<B->n+1; i++){
	    C->p[i+A->n]=Anzmax+B->p[i];
	}
    }else{
	error("Wrong type\n");
    }
    return C;
}
void spsetnzmax(dsp *sp, long nzmax){
    if(sp->nzmax!=nzmax){
	sp->i=realloc(sp->i, sizeof(long)*nzmax);
	sp->x=realloc(sp->x, sizeof(double)*nzmax);
	sp->nzmax=nzmax;
    }
}
void spfree(dsp *sp){
    if(sp){
	if(sp->shmkey==0){
	    free(sp->x);
	    free(sp->p);
	    free(sp->i);
	}
	free(sp);
    }
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    enum{
	P_XLOC,/*XLOC is where phase points are defined. XLOC==PLOC if on same plan*/
	P_DX,
	P_PLOC,/*PLOC is where AMP is defined.*/
	P_DP,
	P_AMP,
	P_SALOC,
	P_SAORC, /*saorc: SALOC is subaperture origin or center.
		   1: origin (lower left corner), 0: center.*/
	P_DSA,
	P_SCALE,
	P_DISPLACE,
	P_DOPARTIAL,
	P_TOT,
    };
    int do_partial;
    int saorc;
    LOC_T xloc,ploc,saloc;
    double *amp;
    double scale, *displace;
    if(nrhs!=P_TOT){
	printf("Number of inputs expected: %d\n", P_TOT);
	mexErrMsgTxt("Wrong number of inputs!");
    }
    if(mxGetNumberOfElements(prhs[P_DISPLACE])!=2){
	mexErrMsgTxt("Displace needs to be a two-vector.\n");
    }
#define LOC_FROM_MATLAB(A,B)						\
    {A.nloc=mxGetM(B); A.locx=mxGetPr(B); A.locy=A.locx+A.nloc;A.map=NULL;}
    LOC_FROM_MATLAB(xloc,prhs[P_XLOC]);
    LOC_FROM_MATLAB(ploc,prhs[P_PLOC]);
    LOC_FROM_MATLAB(saloc,prhs[P_SALOC]);
#undef LOC_FROM_MATLAB
    amp=mxGetPr(prhs[P_AMP]);
    xloc.dx=mxGetScalar(prhs[P_DX]);
    ploc.dx=mxGetScalar(prhs[P_DP]);
    saloc.dx=mxGetScalar(prhs[P_DSA]);
    
    /*scale and displace is between XLOC and PLOC->*/
    scale=mxGetScalar(prhs[P_SCALE]);
    displace=mxGetPr(prhs[P_DISPLACE]);
    saorc=(int)mxGetScalar(prhs[P_SAORC]);
    do_partial=(int)mxGetScalar(prhs[P_DOPARTIAL]);
    dsp *GS0T=mkgt(&xloc, &ploc, amp, &saloc, 
		       saorc, scale, displace, do_partial);
    mxArray *GS0T2=mxCreateSparse(GS0T->m, GS0T->n, GS0T->nzmax, mxREAL);
    memcpy(mxGetPr(GS0T2), GS0T->x, GS0T->nzmax*sizeof(double));
    memcpy(mxGetIr(GS0T2), GS0T->i, GS0T->nzmax*sizeof(long));
    memcpy(mxGetJc(GS0T2), GS0T->p, (GS0T->n+1)*sizeof(long));
    spfree(GS0T);
    mexCallMATLAB(1,plhs,1, &GS0T2,"tanspose");
}

#ifdef __INTEL_COMPILER
#undef _GNU_SOURCE /*avoid compiling problem*/
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "proj.h"
#define SPLIT(A,B,C) {C=(int)floor(A); B=A-C;}
const double pi=3.1415926535897932384626433832795;
static __inline double cosangle(double a[3], double b[3]){
    return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
	/sqrt((a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
	      *(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]));
}
static void 
proj_rect_grid(rectmap_t *mapin, double thetax, double thetay,
	       const loc_t *locout,const double ratiox, const double ratioy,
	       const double *ampout, double* phiout, 
	       double sc, double hs, 
	       double betax, double betay){
    /*
      input parameters:
      mapin: M3 surface map. 
      thetax, thetay: tilt of surface map along x or y. one has to be pi/2
      locout: Pupil plan opd grid.
      ratio: scaling of pupil plan to exit pupil plane.
      sc: scaling of opd. don't include projection
      hs: distance from exit pupil to m3.
      betax,betay: beam angle
     */
    double ht=mapin->h;
    const double (*phiin)[mapin->nx]=(void*)mapin->p;
    const int wrapx1 = mapin->nx;
    const int wrapy1 = mapin->ny;
    const int wrapx = wrapx1-1;
    const int wrapy = wrapy1-1;
    double offx=-hs*betax;
    double offy=-hs*betay;
    if(ratiox<0) offx=-offx;
    if(ratioy<0) offy=-offy;
    const double dx_in1 = 1./mapin->dx;
    const double dy_in1 = 1./mapin->dy;
    double a0x=(-offx/sin(thetax)-mapin->ox)*dx_in1;
    double a0y=(-offy/sin(thetay)-mapin->oy)*dy_in1;
    double ddx=(hs-ht)*dx_in1;
    double ddy=(hs-ht)*dy_in1;
  
    int nplocx,nplocy,nplocx1,nplocy1;
    if(fabs(thetax-pi*0.5)>pi*.45 || fabs(thetay-pi*0.5)>pi*0.45){
	fprintf(stderr,"Tilting angle is too much\n");
	return;
    }
    
    double vm3[3];
    vm3[0]=-sin(pi/2-thetax);
    vm3[1]=-sin(pi/2-thetay);
    vm3[2]=-sqrt(1.-vm3[0]*vm3[0]-vm3[1]*vm3[1]);
    double vi[3];
    vi[2]=-hs;
    double sc2;
    int iloc;
    for(iloc=0; iloc<locout->nloc; iloc++){
	if(ampout && fabs(ampout[iloc])<1.e-10)
	    continue;/*skip points that has zero amplitude*/
	double alx=atan2(locout->locx[iloc]*ratiox+offx,hs);
	double aly=atan2(locout->locy[iloc]*ratioy+offy,hs);
	double btx=thetax-alx;
	double bty=thetay-aly;
	double dplocx=ddx*sin(alx)/sin(btx)+a0x;
	double dplocy=ddy*sin(aly)/sin(bty)+a0y;
	vi[0]=locout->locx[iloc]*ratiox+offx;
	vi[1]=locout->locy[iloc]*ratioy+offy;
	
	SPLIT(dplocx,dplocx,nplocx);
	SPLIT(dplocy,dplocy,nplocy);
	
	if(nplocx<0||nplocx>=wrapx||nplocy<0||nplocy>=wrapy){
	    continue;
	}else{
	    nplocx1=nplocx+1;
	    nplocy1=nplocy+1;
	}
	sc2=sc*cosangle(vi,vm3);
	
	phiout[iloc]+=sc2*(phiin[nplocy][nplocx]*(1.-dplocx)*(1.-dplocy)
			  +phiin[nplocy][nplocx1]*(dplocx)*(1.-dplocy)
			  +phiin[nplocy1][nplocx]*(1.-dplocx)*(dplocy)
			  +phiin[nplocy1][nplocx1]*(dplocx)*(dplocy));
    }
}
#ifdef MATLAB_MEX_FILE
#include <mex.h>
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    enum{
	P_SURF,
	P_X,
	P_Y,
	P_ALX,/*M3 tilt X . pi/2 is no tilt*/
	P_ALY,/*M3 tilt Y */
	P_THETAX,/*guide star offset*/
	P_THETAY,
	P_LOC,
	P_AMP,
	P_TOT,
    };
    enum{
	PL_OPD,
	PL_TOT,
    };
    if(P_TOT!=nrhs) mexErrMsgTxt("Incorrect input arguments");
    if(PL_TOT!=nlhs) mexErrMsgTxt("Incorrect output arguments");
    rectmap_t *mapin=calloc(1, sizeof(rectmap_t));
    mapin->p=mxGetPr(prhs[P_SURF]);
    mapin->nx=mxGetM(prhs[P_SURF]);
    mapin->ny=mxGetN(prhs[P_SURF]);
    double alx=mxGetScalar(prhs[P_ALX]);
    double aly=mxGetScalar(prhs[P_ALY]);
    double *X=mxGetPr(prhs[P_X]);
    double *Y=mxGetPr(prhs[P_Y]);
    mapin->dx=X[1]-X[0];
    mapin->dy=Y[1]-Y[0];
    mapin->ox=X[0];
    mapin->oy=Y[0];
    double bx=mxGetScalar(prhs[P_THETAX]);
    double by=mxGetScalar(prhs[P_THETAY]);;
    if (X[mapin->nx-1]-X[0]-(mapin->nx-1)*mapin->dx > 1.e-10)
	mexErrMsgTxt("X has non even spacing\n");
    if (Y[mapin->ny-1]-Y[0]-(mapin->ny-1)*mapin->dy > 1.e-10)
	mexErrMsgTxt("Y has non even spacing\n");
    
    double d_m3_f=20.;/*from m3 to focus*/
    double d_exitpupil_f=46.38661051;
    double d_exitpupil_m3=d_exitpupil_f-d_m3_f;
    double r_exitpupil=1.546220350;
    double r_pupil=15;
    mapin->h=d_exitpupil_m3;
    loc_t* loc2=calloc(1, sizeof(loc_t));
    loc2->locx=mxGetPr(prhs[P_LOC]);
    loc2->nloc=mxGetM(prhs[P_LOC]);
    loc2->locy=loc2->locx+loc2->nloc;
    loc2->dx=loc2->locx[1]-loc2->locx[0];
    
    double *amp=mxGetPr(prhs[P_AMP]);
    plhs[PL_OPD]=mxCreateDoubleMatrix(loc2->nloc,1,mxREAL);
    double *phi2=mxGetPr(plhs[PL_OPD]);
    proj_rect_grid(mapin,alx,aly,
		   loc2,-r_exitpupil/r_pupil,r_exitpupil/r_pupil,
		   amp,phi2,-2,d_exitpupil_f,bx,by);
    free(loc2);
}
#endif

/*Unused routines buried here*/

/*
  y=A'*x where A is sparse. x, y are vectors
*/
__global__ void cusptmul_do(Real *y, int icol, cusp *A, Real *x, Real alpha){
    __shared__ Real val;
    if(threadIdx.x==0) val=0;
    int i=blockIdx.x * blockDim.x + threadIdx.x;
    int j=i+A->p[icol];
    atomicAdd(&val, A->x[j]*x[A->i[j]]);
    if(threadIdx.x==0) y[icol]+=val*alpha;
}

__global__ void cuspmul_do(Real *y, cusp *A, Real *x, Real alpha){
    int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<A->ny; i+=step){
	for(int j=A->p[i]; j<A->p[i+1]; j++){
	    atomicAdd(&y[A->i[j]], A->x[j]*x[i]*alpha);
	}
    }
}
in gpu_prop_grid_do()
	int match2=fabs(xratio-0.5)<EPS && fabs(yratio-0.5)<EPS
		    && fabs(dispx)<EPS && fabs(dispy)<EPS;
		if(match2 && trans!='t' && 0){
		    Real fxv[4], fxm[4];
		    Real x;
		    x=0;
		    fxv[0]=(1.f-x)*(1.f-x)*(cc[3]+cc[4]*(1.f-x));			
		    fxv[1]=cc[0]+x*x*(cc[1]+cc[2]*x);			
		    fxv[2]=cc[0]+(1.f-x)*(1.f-x)*(cc[1]+cc[2]*(1.f-x));			
		    fxv[3]=x*x*(cc[3]+cc[4]*x);		
		    x=0.5;
		    fxm[0]=(1.f-x)*(1.f-x)*(cc[3]+cc[4]*(1.f-x));			
		    fxm[1]=cc[0]+x*x*(cc[1]+cc[2]*x);			
		    fxm[2]=cc[0]+(1.f-x)*(1.f-x)*(cc[1]+cc[2]*(1.f-x));			
		    fxm[3]=x*x*(cc[3]+cc[4]*x);	
		    Real *fx, *fy;
		    int pps_xmin=-datai->offpsx-1;
		    int pps_ymin=-datai->offpsy-1;
		    int pps_xmax=datai->nxps-datai->offpsx;
		    int pps_ymax=datai->nxps-datai->offpsy;
		    if(trans=='t'){
		

		    }else{//32 us without test. 45 us with test.
			for(int my=iy0; my<ny; my+=stepy){
			    int iy=my/2;
			    if((my & 1) == 0){
				fy=fxv;
			    }else{
				fy=fxm;
			    }
			    for(int mx=ix0; mx<nx; mx+=stepx){
				int ix=mx/2;
				if((mx & 1) == 0){
				    fx=fxv;
				}else{
				    fx=fxm;
				}
				Real sum=0;
#pragma unroll
				for(int ky=-1; ky<3; ky++){
				    int kyi=ky+iy;
				    if((kyi)>pps_ymin && (kyi)<pps_ymax){
#pragma unroll
					for(int kx=-1; kx<3; kx++){
					    int kxi=kx+ix;
					    if((kxi)>pps_xmin && (kxi)<pps_xmax){
						sum+=fx[kx+1]*fy[ky+1]*pps[(iy+ky)*nxps+kx+ix];
					    }
					}
				    }
				}
				pdir[mx+my*nxdir]+=sum*alpha;
			    }
			}
		    }
		}else 





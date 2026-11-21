/*
  Copyright 2009-2026 Lianqi Wang

  This file is part of Multithreaded Adaptive Optics Simulator (MAOS).

  MAOS is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  MAOS is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  MAOS.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
  Tests of the building blocks of the Fourier domain preconditioner (FDPCG):
  aperture selection, permutation, gradient and propagation operators in the
  Fourier domain, and assembly and block inversion of the preconditioner with
  csp operations. Ported from test_fd.c, which read saloc.bin and SAA.bin and
  wrote the results to files. The helper functions are local to the test.

  Fixes to the helpers relative to test_fd.c:
  - fdpcg_sa: frequencies whose magnitude equals that of the zero frequency
    (the aliases) are kept with a relative tolerance instead of an exact float
    comparison which could drop them due to round off.
  - fdpcg_g: the gradient was identically zero for os==1 because the phase
    offset was initialized to 0 instead of 1.
  - fdpcg_prop: the coarse layer interpolation used a kernel with only the
    diagonal neighbors evaluated at the layer (aliased) frequency. It now uses
    the bilinear kernel at the ground frequency, as in maos/fdpcg.c.
  - fdpcg_prop: VLAs are replaced by allocated arrays.
*/
#include "unittest.h"

#if CPU_SINGLE
#define TOL 1e-4
#else
#define TOL 1e-10
#endif

/**
   Create aperture selection function that selects the gradients for valid
   subapertures from ground layer xloc (or ploc) in Fourier domain.
*/
static csp *fdpcg_sa(loc_t *xloc, loc_t *saloc, real *saa){
	const real threas=1;
	long nxloc=xloc->nloc;
	long nx=(long)sqrt((real)nxloc);
	long ny=nx;
	if(nxloc!=nx*ny){
		warning("xloc is not square. The code may not work\n");
	}
	cmat *xsel=cnew(nx, ny);
	real dx1=1./xloc->dx;
	long offx=(long)round(-xloc->locx[0]*dx1+saloc->dx*0.5*dx1);
	long offy=(long)round(-xloc->locy[0]*dx1+saloc->dx*0.5*dx1);
	for(long isa=0; isa<saloc->nloc; isa++){
		if(saa[isa]>0.9){
			long ix=(long)round(saloc->locx[isa]*dx1)+offx;/*subaperture center. */
			long iy=(long)round(saloc->locy[isa]*dx1)+offy;
			P(xsel, ix, iy)=1;
		}
	}
	cfftshift(xsel);
	cfft2(xsel, -1);
	cfftshift(xsel);
	cscale(xsel, 1./(real)(nx*ny));
	real xselc=creal(P(xsel, nx/2, ny/2))*threas;/*Fourier center */
	for(long ix=0; ix<nx; ix++){
		for(long iy=0; iy<ny; iy++){
			if(cabs(P(xsel, ix, iy))<xselc*(1-1e-10)){
				P(xsel, ix, iy)=0;
			}
		}
	}
	csp *sel=cspconvolvop(xsel);
	cfree(xsel);
	return sel;
}
/**
   Create a permulation acting on xhat (fft of x on xloc) so that points of
   the same spatial frequency in all layers are grouped together.
*/
static long *fdpcg_perm(long *nperm, loc_t **xloc, int nps, loc_t *saloc){
	long *xdim=mycalloc(nps, long);
	long *os=mycalloc(nps, long);
	long *fxlim=mycalloc(nps, long);
	long *noff=mycalloc(nps, long);
	long xloctot=0;
	for(long ips=0; ips<nps; ips++){
		xdim[ips]=(long)sqrt((real)xloc[ips]->nloc);
		fxlim[ips]=xdim[ips]/2;
		noff[ips]=xloctot;
		xloctot+=xloc[ips]->nloc;
		if(xloc[ips]->nloc!=xdim[ips]*xdim[ips]){
			error("xloc must be square\n");
		}
		os[ips]=(long)round(saloc->dx/xloc[ips]->dx);
		if(os[ips]>os[0]){
			error("Layer %ld oversampling ratio is greater than ground layer\n", ips);
		}
	}
	long *perm=mycalloc(xloctot, long);
	long use_os=os[0];
	long adim=xdim[0]/use_os;
	long osx=xdim[0]/2;
	long count=0;
	for(long iy=-adim/2; iy<adim/2; iy++){
		for(long ix=-adim/2; ix<adim/2; ix++){
			/*(ix,iy) is the frequency in SALOC grid. */
			for(long juse_os=0; juse_os<use_os; juse_os++){
				long jy=(iy+adim*juse_os);
				if(jy>=osx) jy-=xdim[0];
				for(long iuse_os=0; iuse_os<use_os; iuse_os++){
					long jx=(ix+adim*iuse_os);
					if(jx>=osx) jx-=xdim[0];
					/*jx, jy is the frequency in XLOC grid. */
					for(long ips=0; ips<nps; ips++){
						if(jy>=-fxlim[ips]&&jy<fxlim[ips] /*this layer has such freq. */
							&&jx>=-fxlim[ips]&&jx<fxlim[ips]){
							perm[count]=noff[ips]+(jx+fxlim[ips])+(jy+fxlim[ips])*xdim[ips];
							count++;
						}
					}
				}
			}
		}
	}
	free(fxlim);
	free(noff);
	free(xdim);
	free(os);
	*nperm=count;
	return perm;
}
/**
   Compute gradient operator in Fourier domain (zero frequency at the center).
*/
static void fdpcg_g(cmat **gx, cmat **gy, long nx, long ny, real dx, real dsa){
	long os=(long)round(dsa/dx);
	if(fabs(dsa-dx*os)>1.e-10){
		error("dsa must be multiple of dx");
	}
	real *wt=mymalloc(os+1, real);
	real *st=mymalloc(os+1, real);
	/*Trapzoidal weights for averaging. */
	wt[os]=wt[0]=0.5/(real)os/dsa;
	for(long ios=1; ios<os; ios++){
		wt[ios]=1./(real)os/dsa;
	}
	for(long ios=0; ios<os+1; ios++){
		st[ios]=ios*dx;
	}
	long ny2=ny/2;
	long nx2=nx/2;
	comp cf=COMPLEX(0, 2*M_PI);
	real dfy=1/(ny*dx);
	real dfx=1/(nx*dx);
	*gx=cnew(nx*ny, 1);
	*gy=cnew(nx*ny, 1);
	comp *pgx=P(*gx);
	comp *pgy=P(*gy);
	real dsa2=dsa*0.5;
	for(long iy=0; iy<ny; iy++){
		real fy=(real)(iy-ny2)*dfy;
		for(long ix=0; ix<nx; ix++){
			real fx=(real)(ix-nx2)*dfx;
			comp tx=0;
			comp ty=0;
			comp offset=1;
			if(os>1){
				offset=cexp(-cf*(fx+fy)*dsa2);/*shift by half a subaperture */
			}
			for(int ios=0; ios<os+1; ios++){
				tx+=wt[ios]*(cexp(cf*(fx*dsa+fy*st[ios]))-cexp(cf*(fy*st[ios])));
				ty+=wt[ios]*(cexp(cf*(fy*dsa+fx*st[ios]))-cexp(cf*(fx*st[ios])));
			}
			pgx[ix+iy*nx]=offset*tx;
			pgy[ix+iy*nx]=offset*ty;
		}
	}
	free(wt);
	free(st);
}
/**
   Propagate operator for nlayer screens to ground of size nxg*nxg, sampling
   dx, with displacement of dispx, dispy. Layers with os[ips]=os[0]/2 are
   sampled at 2*dx and are interpolated bilinearly. Assumes unitary FFT.
*/
static csp *fdpcg_prop(long nps, const long *os, long nxg, real dx, const real *dispx, const real *dispy){
	long *nxi=mycalloc(nps, long);
	long *nxi2=mycalloc(nps, long);
	long *nxi3=mycalloc(nps, long);
	long *noff=mycalloc(nps, long);
	long nxtot=0;
	real dk=1./(nxg*dx);
	for(long ips=0; ips<nps; ips++){
		nxi[ips]=nxg/os[0]*os[ips];
		nxi2[ips]=nxi[ips]/2;
		nxi3[ips]=nxi[ips]/2+nxi[ips];
		noff[ips]=nxtot;
		nxtot+=nxi[ips]*nxi[ips];
	}
	long nxg2=nxg/2;
	csp *propt=cspnew(nxtot, nxg*nxg, nxg*nxg*nps);
	spint *pp=propt->pp;
	spint *pi=propt->pi;
	comp *px=propt->px;
	long count=0;
	comp cf=COMPLEX(0, 2*M_PI);
	for(long iy=0; iy<nxg; iy++){
		real fyg=(iy-nxg2)*dk;
		for(long ix=0; ix<nxg; ix++){
			real fxg=(ix-nxg2)*dk;
			long icol=ix+iy*nxg;
			pp[icol]=count;
			for(long ips=0; ips<nps; ips++){
				long jx=((ix-nxg2)+nxi3[ips])%nxi[ips];/*map to layer ips. */
				long jy=((iy-nxg2)+nxi3[ips])%nxi[ips];/*map to layer ips. */
				real fx=(jx-nxi2[ips])*dk;
				real fy=(jy-nxi2[ips])*dk;
				pi[count]=jx+jy*nxi[ips]+noff[ips];
				comp shift=cexp(cf*(fx*dispx[ips]+fy*dispy[ips]));
				switch(os[0]/os[ips]){
				case 1:
					px[count]=conj(shift);
					break;
				case 2:
					//bilinear interpolation. 0.5 accounts for the unitary FFT scaling.
					px[count]=conj(shift)*((1.+cos(2*M_PI*(fxg*dx)))*(1.+cos(2*M_PI*(fyg*dx)))*0.5);
					break;
				default:
					error("Invalid\n");
				}
				count++;
			}
		}
	}
	pp[nxg*nxg]=count;
	/*we put conj above because csptrans applies a conjugation. */
	csp *prop0=csptrans(propt);
	cspfree(propt);
	free(nxi); free(nxi2); free(nxi3); free(noff);
	return prop0;
}
/*Unitary 2d FFT with zero frequency at the center.*/
static cmat *fft_centered(const dmat *in){
	cmat *out=cnew(NX(in), NY(in));
	for(long i=0; i<PN(in); i++){
		P(out, i)=P(in, i);
	}
	cfft2(out, -1);
	cfftshift(out);
	cscale(out, 1./sqrt((real)PN(in)));
	return out;
}
/*Maximum absolute value of a complex sparse matrix*/
static real csp_maxabs(const csp *A){
	real maxv=0;
	for(long i=0; i<A->pp[NY(A)]; i++){
		if(cabs(A->px[i])>maxv) maxv=cabs(A->px[i]);
	}
	return maxv;
}
/*Count of nonzero values in column icol of A*/
static long csp_colnnz(const csp *A, long icol){
	return A->pp[icol+1]-A->pp[icol];
}
/*Find value of A(irow, icol). Returns 0 if not found.*/
static comp csp_get(const csp *A, long irow, long icol){
	for(long k=A->pp[icol]; k<A->pp[icol+1]; k++){
		if(A->pi[k]==irow) return A->px[k];
	}
	return 0;
}
/*Geometry: N*N xloc grid with spacing dx, subapertures of 2*dx in a circular pupil.*/
#define N 32
#define DX 0.25
#define DSA 0.5
static void make_geom(loc_t **xloc, loc_t **saloc, dmat **saa, real D){
	*xloc=mksqloc_auto(N, N, DX, DX);
	const long nsa=(long)(N*DX/DSA);
	*saloc=mksqloc_auto(nsa, nsa, DSA, DSA);
	*saa=dnew((*saloc)->nloc, 1);
	for(long isa=0; isa<(*saloc)->nloc; isa++){
		real cx=(*saloc)->locx[isa]+DSA*0.5;
		real cy=(*saloc)->locy[isa]+DSA*0.5;
		P(*saa, isa)=(sqrt(cx*cx+cy*cy)<D/2)?1:0;
	}
}

/** The aperture selection operator is a convolution with the zero frequency and its aliases of equal magnitude. */
static void test_sa(void **state){
	(void)state;
	for(int ifull=0; ifull<2; ifull++){
		loc_t *xloc, *saloc;
		dmat *saa;
		make_geom(&xloc, &saloc, &saa, ifull?100:7);
		const real nvalid=dsum(saa);
		csp *sel=fdpcg_sa(xloc, saloc, P(saa));
		assert_int_equal(NX(sel), N*N);
		assert_int_equal(NY(sel), N*N);
		const real dc=nvalid/(N*N);
		//The zero frequency of the subaperture centers sampled at 2*dx is aliased to N/2.
		const long offs[4][2]={{0,0},{N/2,0},{0,N/2},{N/2,N/2}};
		for(long icol=0; icol<N*N; icol++){
			assert_int_equal(csp_colnnz(sel, icol), 4);
			long ix=icol%N, iy=icol/N;
			for(int io=0; io<4; io++){
				long irow=(ix+offs[io][0])%N+((iy+offs[io][1])%N)*N;
				comp val=csp_get(sel, irow, icol);
				assert_float_equal(cabs(val), dc, TOL*dc);
				if(io==0){
					assert_float_equal(creal(val), dc, TOL*dc);//diagonal is the fill factor
				}
				assert_float_equal(fabs(creal(val)), dc, TOL*dc);//aliases are real (+-1)
			}
		}
		if(ifull){
			assert_float_equal(dc, 0.25, 0);//all subapertures are valid
		}
		cspfree(sel);
		dfree(saa);
		locfree(saloc);
		locfree(xloc);
	}
}
/** The permutation is a valid permutation that groups each frequency with its aliases in all layers. */
static void test_perm(void **state){
	(void)state;
	loc_t *xloc, *saloc;
	dmat *saa;
	make_geom(&xloc, &saloc, &saa, 7);
	loc_t *xloc2=mksqloc_auto(N/2, N/2, DX*2, DX*2);//coarse layer
	for(int imix=0; imix<2; imix++){
		const int nps=3;
		loc_t *xlocs[3]={xloc, imix?xloc2:xloc, xloc};
		long nperm=0;
		long *perm=fdpcg_perm(&nperm, xlocs, nps, saloc);
		long ntot=0;
		for(int ips=0; ips<nps; ips++) ntot+=xlocs[ips]->nloc;
		assert_int_equal(nperm, ntot);
		int *seen=mycalloc(ntot, int);
		for(long i=0; i<nperm; i++){
			assert_true(perm[i]>=0&&perm[i]<ntot);
			seen[perm[i]]++;
		}
		for(long i=0; i<ntot; i++){
			assert_int_equal(seen[i], 1);
		}
		free(seen);
		//block size: 4 aliases for each fine layer, 1 for coarse layer.
		const long bs=imix?(4+1+4):(4*nps);
		assert_int_equal(nperm%bs, 0);
		const long adim=N/2;
		for(long ib=0; ib<nperm/bs; ib++){
			long fx0=0, fy0=0;
			for(long k=0; k<bs; k++){
				long ind=perm[ib*bs+k];
				long ips=0;
				while(ind>=xlocs[ips]->nloc){
					ind-=xlocs[ips]->nloc;
					ips++;
				}
				long nd=(long)sqrt((real)xlocs[ips]->nloc);
				long fx=(ind%nd)-nd/2;
				long fy=(ind/nd)-nd/2;
				//frequency modulo the subaperture sampling frequency
				fx=((fx%adim)+adim)%adim;
				fy=((fy%adim)+adim)%adim;
				if(k==0){
					fx0=fx; fy0=fy;
				}else{
					assert_int_equal(fx, fx0);
					assert_int_equal(fy, fy0);
				}
			}
		}
		free(perm);
	}
	locfree(xloc2);
	dfree(saa);
	locfree(saloc);
	locfree(xloc);
}
/** Fourier domain gradient operators vanish at zero frequency and match the closed form of averaged gradient. */
static void test_g(void **state){
	(void)state;
	for(int ios=0; ios<2; ios++){
		const real dsa=ios?DSA:DX;//os=2 or os=1
		cmat *gx, *gy;
		fdpcg_g(&gx, &gy, N, N, DX, dsa);
		const real dk=1./(N*DX);
		assert_float_equal(cabs(P(gx, N/2+N/2*N)), 0, TOL);
		assert_float_equal(cabs(P(gy, N/2+N/2*N)), 0, TOL);
		for(long iy=0; iy<N; iy++){
			real fy=(iy-N/2)*dk;
			for(long ix=0; ix<N; ix++){
				real fx=(ix-N/2)*dk;
				//symmetry between x and y
				assert_float_equal(cabs(P(gx, ix+iy*N)-P(gy, iy+ix*N)), 0, TOL);
				//conjugate symmetry of a real operator
				if(ix>0&&iy>0){
					assert_float_equal(cabs(P(gx, ix+iy*N)-conj(P(gx, (N-ix)+(N-iy)*N))), 0, TOL);
				}
				//along fy=0: difference of the two edges divided by dsa
				if(iy==N/2){
					comp expect=ios?COMPLEX(0, 2*sin(M_PI*fx*dsa)/dsa):(cexp(COMPLEX(0, 2*M_PI*fx*dsa))-1.)/dsa;
					assert_float_equal(cabs(P(gx, ix+iy*N)-expect), 0, TOL);
				}
				(void)fy;
			}
		}
		//Low frequency limit: |g|=2*pi*f (os=1 is not centered and has an extra phase of pi*f*dsa)
		comp g1=P(gx, N/2+1+N/2*N);
		assert_float_equal(cabs(g1)/(2*M_PI*dk), 1, 0.01);
		if(ios){
			assert_float_equal(cabs(g1/COMPLEX(0, 2*M_PI*dk)-1.), 0, 0.01);
		}
		cfree(gx);
		cfree(gy);
	}
}
/** The Fourier propagation operator matches shifting and bilinear upsampling done in the spatial domain. */
static void test_prop(void **state){
	(void)state;
	const long nps=2;
	const long os[2]={2, 1};//layer 1 is sampled at 2*DX
	//zero displacement: identity for the fine layer
	{
		const real disp[2]={0, 0};
		csp *prop=fdpcg_prop(nps, os, N, DX, disp, disp);
		assert_int_equal(NX(prop), N*N);
		assert_int_equal(NY(prop), N*N+N*N/4);
		for(long i=0; i<N*N; i++){
			assert_float_equal(cabs(csp_get(prop, i, i)-1.), 0, TOL);
		}
		cspfree(prop);
	}
	rand_t rstat;
	seed_rand(&rstat, 1);
	dmat *l0=dnew(N, N);
	dmat *l1=dnew(N/2, N/2);
	drandn(l0, 1, &rstat);
	drandn(l1, 1, &rstat);
	const long k0x=3, k0y=-2, k1x=-1, k1y=2;//integer shifts (in pixels of each layer)
	const real dispx[2]={k0x*DX, k1x*2*DX};
	const real dispy[2]={k0y*DX, k1y*2*DX};
	//spatial domain: g(x)=l0(x+d0)+interp(l1)(x+d1) with periodic boundary
	dmat *up=dnew(N, N);
	for(long iy=0; iy<N; iy++){
		for(long ix=0; ix<N; ix++){
			long mx=ix/2, my=iy/2;
			long mx1=(mx+1)%(N/2), my1=(my+1)%(N/2);
			real wx=(ix%2)*0.5, wy=(iy%2)*0.5;
			P(up, ix, iy)=(1-wx)*(1-wy)*P(l1, mx, my)+wx*(1-wy)*P(l1, mx1, my)
				+(1-wx)*wy*P(l1, mx, my1)+wx*wy*P(l1, mx1, my1);
		}
	}
	dmat *g=dnew(N, N);
	for(long iy=0; iy<N; iy++){
		for(long ix=0; ix<N; ix++){
			P(g, ix, iy)=P(l0, (ix+k0x+N)%N, (iy+k0y+N)%N)
				+P(up, (ix+2*k1x+N)%N, (iy+2*k1y+N)%N);
		}
	}
	cmat *ghat=fft_centered(g);
	//Fourier domain
	cmat *l0hat=fft_centered(l0);
	cmat *l1hat=fft_centered(l1);
	cmat *lhat=cnew(N*N+N*N/4, 1);
	memcpy(P(lhat), P(l0hat), sizeof(comp)*N*N);
	memcpy(P(lhat)+N*N, P(l1hat), sizeof(comp)*N*N/4);
	csp *prop=fdpcg_prop(nps, os, N, DX, dispx, dispy);
	cmat *ghat2=cnew(N*N, 1);
	cspmv(ghat2, prop, lhat, 'n', 1);
	real maxd=0;
	for(long i=0; i<N*N; i++){
		real d=cabs(P(ghat2, i)-P(ghat, i));
		if(d>maxd) maxd=d;
	}
	assert_float_equal(maxd, 0, TOL);
	cspfree(prop);
	cfree(ghat2); cfree(lhat); cfree(l0hat); cfree(l1hat); cfree(ghat);
	dfree(g); dfree(up); dfree(l0); dfree(l1);
}
/** The assembled preconditioner is Hermitian and block diagonal after permutation and cspinvbdiag inverts each block. */
static void test_mhat(void **state){
	(void)state;
	loc_t *xloc, *saloc;
	dmat *saa;
	make_geom(&xloc, &saloc, &saa, 7);
	const long nps=3;
	const long nwfs=3;
	loc_t *xlocs[3]={xloc, xloc, xloc};
	const long os[3]={2, 2, 2};
	const real ht[3]={0, 4000, 10000};
	const real wt[3]={0.5, 0.3, 0.2};
	const real hs=90000;
	const real thetax[3]={0, -0.0001613797784, 9.973820086e-05};
	const real thetay[3]={0.000169684629, 5.243545924e-05, -0.0001372777737};
	const real NEA=1.712908313e-07;
	const real r0=0.198749305619780;
	const long nxtot=N*N*nps;
	csp *sel=fdpcg_sa(xloc, saloc, P(saa));
	long nperm=0;
	long *perm=fdpcg_perm(&nperm, xlocs, nps, saloc);
	assert_int_equal(nperm, nxtot);
	cmat *gx, *gy;
	fdpcg_g(&gx, &gy, N, N, DX, DSA);
	//inverse PSD (zero frequency at the center)
	comp *invpsd=mycalloc(nxtot, comp);
	for(long ips=0; ips<nps; ips++){
		dmat *tmp=turbpsd(N, N, DX*(1-ht[ips]/hs), r0*pow(wt[ips], -3./5.), 30, 0, -1);
		P(tmp, 0)=P(tmp, 1);//turbpsd sets piston to 0. Use a finite value to make Mhat invertible.
		dfftshift(tmp);
		dscale(tmp, N*N);
		for(long i=0; i<N*N; i++){
			invpsd[ips*N*N+i]=P(tmp, i);
		}
		dfree(tmp);
	}
	csp *Mpsd=cspnewdiag(nxtot, invpsd, 1);
	free(invpsd);
	/*Compute gx'*sel'*sel*gx+gy'*sel'*sel*gy */
	csp *Mmid=NULL;
	for(int i=0; i<2; i++){
		csp *tmp=cspdup(sel);
		cspmuldiag(tmp, P(i==0?gx:gy), 1);
		csp *tmp2=cspmulsp(tmp, tmp, "cn");
		cspadd(&Mmid, 1, tmp2, 1);
		cspfree(tmp);
		cspfree(tmp2);
	}
	csp *Mhat=NULL;
	real *dispx=mycalloc(nps, real);
	real *dispy=mycalloc(nps, real);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
		for(long ips=0; ips<nps; ips++){
			dispx[ips]=ht[ips]*thetax[iwfs]/(1-ht[ips]/hs);
			dispy[ips]=ht[ips]*thetay[iwfs]/(1-ht[ips]/hs);
		}
		csp *prop1=fdpcg_prop(nps, os, N, DX, dispx, dispy);
		cspscale(prop1, 1./NEA);
		csp *tmp=cspmulsp(Mmid, prop1, "nn");
		csp *tmp2=cspmulsp(prop1, tmp, "cn");
		cspadd(&Mhat, 1, tmp2, 1);
		cspfree(tmp);
		cspfree(tmp2);
		cspfree(prop1);
	}
	free(dispx);
	free(dispy);
	cspadd(&Mhat, 1, Mpsd, 1);
	cspfree(Mpsd);
	assert_int_equal(NX(Mhat), nxtot);
	assert_int_equal(NY(Mhat), nxtot);
	const real mmax=csp_maxabs(Mhat);
	//Hermitian
	{
		csp *MhatH=csptrans(Mhat);
		csp *diff=cspadd2(Mhat, 1, MhatH, -1);
		assert_float_equal(csp_maxabs(diff), 0, 1e-12*mmax);
		cspfree(diff);
		cspfree(MhatH);
	}
	cspdroptol(Mhat, 1e-15);
	cspsym(&Mhat);
	csp *Mhatp=cspperm(Mhat, 0, perm, perm);
	//Mhatp(i,j)=Mhat(perm[i],perm[j])
	{
		cmat *x=cnew(nxtot, 1);
		cmat *xp=cnew(nxtot, 1);
		for(long i=0; i<nxtot; i++){
			P(x, i)=COMPLEX(sin(0.1*i), cos(0.37*i));
		}
		for(long i=0; i<nxtot; i++){
			P(xp, i)=P(x, perm[i]);
		}
		cmat *y=cnew(nxtot, 1);
		cmat *yp=cnew(nxtot, 1);
		cspmv(y, Mhat, x, 'n', 1);
		cspmv(yp, Mhatp, xp, 'n', 1);
		real maxd=0;
		for(long i=0; i<nxtot; i++){
			real d=cabs(P(yp, i)-P(y, perm[i]));
			if(d>maxd) maxd=d;
		}
		assert_float_equal(maxd, 0, 1e-12*mmax*nxtot);
		cfree(x); cfree(xp); cfree(y); cfree(yp);
	}
	//Block diagonal with block size of 4 aliases times nps layers
	const long bs=4*nps;
	for(long icol=0; icol<nxtot; icol++){
		for(long k=Mhatp->pp[icol]; k<Mhatp->pp[icol+1]; k++){
			assert_int_equal(Mhatp->pi[k]/bs, icol/bs);
		}
	}
	//Each block is Hermitian positive definite
	ccell *blocks=cspblockextract(Mhatp, bs);
	assert_int_equal(PN(blocks), nxtot/bs);
	for(long ib=0; ib<PN(blocks); ib++){
		const cmat *bk=P(blocks, ib);
		for(long i=0; i<bs; i++){
			assert_true(creal(P(bk, i, i))>0);
			assert_float_equal(cimag(P(bk, i, i)), 0, 1e-12*mmax);
		}
	}
	//Block inverse
	csp *Minv=cspinvbdiag(Mhatp, bs);
	csp *MM=cspmulsp(Minv, Mhatp, "nn");
	real maxd=0;
	for(long icol=0; icol<nxtot; icol++){
		real dd=cabs(csp_get(MM, icol, icol)-1.);
		if(dd>maxd) maxd=dd;
		for(long k=MM->pp[icol]; k<MM->pp[icol+1]; k++){
			if(MM->pi[k]!=icol && cabs(MM->px[k])>maxd) maxd=cabs(MM->px[k]);
		}
	}
	assert_float_equal(maxd, 0, 1e-8);
	cspfree(MM);
	cspfree(Minv);
	cellfree(blocks);
	cspfree(Mhatp);
	cspfree(Mhat);
	cspfree(Mmid);
	cfree(gx);
	cfree(gy);
	free(perm);
	cspfree(sel);
	dfree(saa);
	locfree(saloc);
	locfree(xloc);
}

UNITTEST_MAIN(
	cmocka_unit_test(test_sa),
	cmocka_unit_test(test_perm),
	cmocka_unit_test(test_g),
	cmocka_unit_test(test_prop),
	cmocka_unit_test(test_mhat)
)

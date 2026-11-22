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
  Tests of the Fourier domain preconditioner (FDPCG) in lib/fdpcg.c: aperture
  selection, permutation, gradient and propagation operators in the Fourier
  domain, assembly of the preconditioner matrix, and the block inversion and
  application of the preconditioner. Ported from test_fd.c, which read
  saloc.bin and SAA.bin and wrote the results to files.
*/
#include "unittest.h"

#if CPU_SINGLE
#define TOL 1e-4
#else
#define TOL 1e-10
#endif

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
		csp *sel=fdpcg_saselect(N, N, DX, saloc, P(saa));
		assert_int_equal(NX(sel), N*N);
		assert_int_equal(NY(sel), N*N);
		const real dc=nvalid/(N*N);
		//The zero frequency of the subapertures sampled at 2*dx is aliased to N/2.
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
	for(int imix=0; imix<2; imix++){
		const int nps=3;
		long nx[3]={N, imix?N/2:N, N};//layer 1 is coarse in the mixed case
		const real os[3]={2, (real)(imix?1:2), 2};
		const int bs=(int)(os[0]*os[0]+os[1]*os[1]+os[2]*os[2]);
		long ntot=0;
		long noff[3];
		for(int ips=0; ips<nps; ips++){
			noff[ips]=ntot;
			ntot+=nx[ips]*nx[ips];
		}
		lmat *perm=fdpcg_perm(nx, nx, os, bs, nps, 0, 0);
		assert_int_equal(PN(perm), ntot);
		int *seen=mycalloc(ntot, int);
		for(long i=0; i<ntot; i++){
			assert_true(P(perm, i)>=0&&P(perm, i)<ntot);
			seen[P(perm, i)]++;
		}
		for(long i=0; i<ntot; i++){
			assert_int_equal(seen[i], 1);
		}
		free(seen);
		//block size: 4 aliases for each fine layer, 1 for coarse layer.
		assert_int_equal(bs, imix?(4+1+4):(4*nps));
		assert_int_equal(ntot%bs, 0);
		const long adim=N/2;
		for(long ib=0; ib<ntot/bs; ib++){
			long fx0=0, fy0=0;
			for(long k=0; k<bs; k++){
				long ind=P(perm, ib*bs+k);
				long ips=nps-1;
				while(ind<noff[ips]) ips--;
				ind-=noff[ips];
				long nd=nx[ips];
				long fx=(ind%nd)-nd/2;//zero frequency at the center (shift=0)
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
		lfree(perm);
		//with fftshift embedded, it is still a permutation
		perm=fdpcg_perm(nx, nx, os, bs, nps, 1, 0);
		seen=mycalloc(ntot, int);
		for(long i=0; i<ntot; i++){
			assert_true(P(perm, i)>=0&&P(perm, i)<ntot);
			seen[P(perm, i)]++;
		}
		for(long i=0; i<ntot; i++){
			assert_int_equal(seen[i], 1);
		}
		free(seen);
		lfree(perm);
		//about half of the frequencies for real to complex FFT. Negative values indicate conjugate.
		perm=fdpcg_perm(nx, nx, os, bs, nps, 1, 1);
		assert_int_equal(PN(perm), (adim/2+1)*adim*bs);
		for(long i=0; i<PN(perm); i++){
			assert_true(labs(P(perm, i))<ntot);
		}
		lfree(perm);
	}
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
			for(long ix=0; ix<N; ix++){
				real fx=(ix-N/2)*dk;
				//symmetry between x and y
				assert_float_equal(cabs(P(gx, ix+iy*N)-P(gy, iy+ix*N)), 0, TOL);
				//conjugate symmetry of a real operator
				if(ix>0&&iy>0){
					assert_float_equal(cabs(P(gx, ix+iy*N)-conj(P(gx, (N-ix)+(N-iy)*N))), 0, TOL);
				}
				//along fy=0: difference of the two edges of the subaperture (at its lower left corner) divided by dsa
				if(iy==N/2){
					comp expect=(cexp(COMPLEX(0, 2*M_PI*fx*dsa))-1.)/dsa;
					assert_float_equal(cabs(P(gx, ix+iy*N)-expect), 0, TOL);
				}
			}
		}
		//Low frequency limit: |g|=2*pi*f
		comp g1=P(gx, N/2+1+N/2*N);
		assert_float_equal(cabs(g1)/(2*M_PI*dk), 1, 0.01);
		cfree(gx);
		cfree(gy);
	}
}
/** The Fourier propagation operator matches shifting and bilinear upsampling done in the spatial domain. */
static void test_prop(void **state){
	(void)state;
	const long nps=2;
	long nx[2]={N, N/2};//layer 1 is sampled at 2*DX
	//zero displacement: identity for the fine layer
	{
		real disp[2]={0, 0};
		csp *propt=fdpcg_prop_trans(nps, N, N, nx, nx, DX, disp, disp);
		csp *prop=csptrans(propt);
		cspfree(propt);
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
	real dispx[2]={k0x*DX, k1x*2*DX};
	real dispy[2]={k0y*DX, k1y*2*DX};
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
	csp *propt=fdpcg_prop_trans(nps, N, N, nx, nx, DX, dispx, dispy);
	csp *prop=csptrans(propt);//csptrans also conjugates
	cspfree(propt);
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
/**
   The assembled preconditioner is Hermitian and block diagonal after
   permutation, and the preconditioner applied in the spatial domain inverts it.
*/
static void test_mhat(void **state){
	(void)state;
	loc_t *xloc, *saloc;
	dmat *saa;
	make_geom(&xloc, &saloc, &saa, 7);
	const long nps=3;
	const long nwfs=3;
	long nx[3]={N, N, N};
	const real os[3]={2, 2, 2};
	const real ht[3]={0, 4000, 10000};
	const real wt[3]={0.5, 0.3, 0.2};
	const real hs=90000;
	const real thetax[3]={0, -0.0001613797784, 9.973820086e-05};
	const real thetay[3]={0.000169684629, 5.243545924e-05, -0.0001372777737};
	const real neai[3]={1.712908313e-07, 1.712908313e-07, 1.712908313e-07};
	const real r0=0.198749305619780;
	const long nxtot=N*N*nps;
	//inverse PSD (zero frequency at the center)
	cmat *invpsd=cnew(nxtot, 1);
	for(long ips=0; ips<nps; ips++){
		dmat *tmp=turbpsd(N, N, DX*(1-ht[ips]/hs), r0*pow(wt[ips], -3./5.), 30, 0, -1);
		P(tmp, 0)=P(tmp, 1);//turbpsd sets piston to 0. Use a finite value to make Mhat invertible.
		dfftshift(tmp);
		dscale(tmp, N*N);
		for(long i=0; i<N*N; i++){
			P(invpsd, ips*N*N+i)=P(tmp, i);
		}
		dfree(tmp);
	}
	dmat *dispx=dnew(nps, nwfs);
	dmat *dispy=dnew(nps, nwfs);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
		for(long ips=0; ips<nps; ips++){
			P(dispx, ips, iwfs)=ht[ips]*thetax[iwfs]/(1-ht[ips]/hs);
			P(dispy, ips, iwfs)=ht[ips]*thetay[iwfs]/(1-ht[ips]/hs);
		}
	}
	csp *Mhat=fdpcg_mhat(nps, nx, nx, N, N, DX, saloc, P(saa), invpsd, nwfs, dispx, dispy, neai, NULL, 0);
	cfree(invpsd);
	dfree(dispx);
	dfree(dispy);
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
	//the preconditioner as created by maos (fdpcg_new frees its input)
	fdpcg_t *fdpcg=fdpcg_new(cspdup(Mhat), nps, nx, nx, os, 1, NULL, 1, 0);
	const long bs=4*nps;//4 aliases times nps layers
	assert_int_equal(fdpcg->bs, bs);
	assert_int_equal(fdpcg->nxtot, nxtot);
	assert_int_equal(fdpcg->scale, 0);
	assert_int_equal(PN(fdpcg->Mbinv), nxtot/bs);
	assert_non_null(fdpcg->permhf);

	csp *Msym=cspdup(Mhat);
	cspdroptol(Msym, 1e-15);
	cspsym(&Msym);
	lmat *perm=fdpcg_perm(nx, nx, os, (int)bs, (int)nps, 0, 0);
	csp *Mhatp=cspperm(Msym, 0, P(perm), P(perm));
	//Mhatp(i,j)=Mhat(perm[i],perm[j])
	{
		cmat *x=cnew(nxtot, 1);
		cmat *xp=cnew(nxtot, 1);
		for(long i=0; i<nxtot; i++){
			P(x, i)=COMPLEX(sin(0.1*i), cos(0.37*i));
		}
		for(long i=0; i<nxtot; i++){
			P(xp, i)=P(x, P(perm, i));
		}
		cmat *y=cnew(nxtot, 1);
		cmat *yp=cnew(nxtot, 1);
		cspmv(y, Msym, x, 'n', 1);
		cspmv(yp, Mhatp, xp, 'n', 1);
		real maxd=0;
		for(long i=0; i<nxtot; i++){
			real d=cabs(P(yp, i)-P(y, P(perm, i)));
			if(d>maxd) maxd=d;
		}
		assert_float_equal(maxd, 0, 1e-12*mmax*nxtot);
		cfree(x); cfree(xp); cfree(y); cfree(yp);
	}
	//Block diagonal
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

	//The preconditioner applied to x gives y with Mhat*U*y=U*x, where U is the unitary FFT with zero frequency at the center.
	{
		rand_t rstat;
		seed_rand(&rstat, 2);
		dcell *x=dcellnew_same(nps, 1, N, N);
		for(long ips=0; ips<nps; ips++){
			drandn(P(x, ips), 1, &rstat);
		}
		dcell *y=NULL;
		fdpcg_apply(&y, fdpcg, x, 1);
		cmat *xhat=cnew(nxtot, 1);
		cmat *yhat=cnew(nxtot, 1);
		for(long ips=0; ips<nps; ips++){
			cmat *tmp=fft_centered(P(x, ips));
			memcpy(P(xhat)+ips*N*N, P(tmp), sizeof(comp)*N*N);
			cfree(tmp);
			tmp=fft_centered(P(y, ips));
			memcpy(P(yhat)+ips*N*N, P(tmp), sizeof(comp)*N*N);
			cfree(tmp);
		}
		cmat *r=cnew(nxtot, 1);
		cspmv(r, Mhat, yhat, 'n', 1);
		real maxr=0, maxx=0;
		for(long i=0; i<nxtot; i++){
			if(cabs(P(r, i)-P(xhat, i))>maxr) maxr=cabs(P(r, i)-P(xhat, i));
			if(cabs(P(xhat, i))>maxx) maxx=cabs(P(xhat, i));
		}
		assert_float_equal(maxr, 0, 1e-6*maxx);
		dcell *y2=NULL;//linear
		dcellscale(x, 2);
		fdpcg_apply(&y2, fdpcg, x, 2);
		dcelladd(&y2, 1, y, -2);
		real maxd2=0, maxy=0;
		for(long ips=0; ips<nps; ips++){
			maxd2=MAX(maxd2, dmaxabs(P(y2, ips)));
			maxy=MAX(maxy, dmaxabs(P(y, ips)));
		}
		assert_float_equal(maxd2, 0, 1e-10*maxy);
		dcellfree(y2);
		cfree(r); cfree(xhat); cfree(yhat);
		dcellfree(y);
		dcellfree(x);
	}
	fdpcg_free(fdpcg);
	cspfree(MM);
	cspfree(Minv);
	cellfree(blocks);
	cspfree(Mhatp);
	cspfree(Msym);
	lfree(perm);
	cspfree(Mhat);
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

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
#include "../math/mathdef.h"
#include "psd.h"
#include "libmisc.h"
//#define W_J(i,N2) (1-pow((real)(i-N2)/(real)N2, 2))
#define W_J(i,N2) 1

/**
   Compute the 1-d PSD from a 1-d sequence.

   @param[in] v The data sequence. A row vector, or one sequence per column.
   @param[in] nseg Number of overlapping segments.
   @return Newly allocated one sided PSD, integrated in each frequency bin, one column per sequence; the caller frees it.
*/
dmat* psd1d(const dmat* v, long nseg){
	long nx;
	long ncol;
	if(NX(v)==1){
		nx=NY(v);
		ncol=1;
	} else{
		nx=NX(v);
		ncol=NY(v);
	}
	if(nseg<=1) nseg=1;
	const long lseg2=nx/(nseg+1);
	const long lseg=lseg2*2;
	dmat* psd=dnew(lseg2+1, ncol);
	cmat* hat=cnew(lseg, 1);
	//cfft2plan(hat, -1);
	for(long icol=0; icol<ncol; icol++){
		real* ppsd=PCOL(psd, icol);
		for(int iseg=0; iseg<nseg; iseg++){
			real* p=P(v)+icol*nx+iseg*lseg2;
			for(long ix=0; ix<lseg; ix++){
				P(hat,ix)=p[ix]*W_J(ix, lseg2);
			}
			cfft2(hat, -1);
			ppsd[0]+=cabs2(P(hat,0));
			for(long ix=1; ix<lseg2; ix++){
				ppsd[ix]+=cabs2(P(hat,ix))+cabs2(P(hat,lseg-ix));
			}
			ppsd[lseg2]+=cabs2(P(hat,lseg2));
		}
	}
	real sumwt=0;
	for(long ix=0; ix<lseg; ix++){
		sumwt+=pow(W_J(ix, lseg2), 2);
	}
	sumwt*=lseg*nseg;
	dscale(psd, 1./sumwt);
	cfree(hat);
	return psd;
}
/**
   Wrap of psd1d to put the frequency along the first column.

   @param[in] v The data sequence. A row vector, or one sequence per column.
   @param[in] nseg Number of overlapping segments.
   @param[in] dt Sampling period (s).
   @return Newly allocated PSD with frequency (Hz) in the first column and PSD density in the remaining; the caller frees it.
*/
dmat* psd1dt(const dmat* v, long nseg, real dt){
	dmat* psd=psd1d(v, nseg);
	dmat* psd2=dnew(NX(psd), psd->ny+1);
	long N=(psd->nx-1)*2;
	real df=1./(N*dt);
	for(long i=0; i<NX(psd); i++){
		P(psd2,i)=df*i;
	}
	dscale(psd, 1./df);//divide so the value is point, not integrated in a bin.
	memcpy(PCOL(psd2, 1), P(psd), PN(psd)*sizeof(real));
	dfree(psd);
	return psd2;
}

/**
   Interpolate psd onto new f. We interpolate in log space which is usually more linear.

   @param[in] psdin The PSD with frequency in the first column and PSD in the second column.
   @param[in] fnew The new frequency.
   @param[in] uselog Interpolate in log space of the PSD.
   @return Newly allocated PSD (without the frequency column) on fnew; the caller frees it.
*/
dmat* psd_interp1(const dmat* psdin, const dmat* fnew, int uselog){
	if(NY(psdin)!=2){
		error("psdin should have 2 columns for frequency and psd.\n");
		return NULL;
	}
	dmat* f1=drefcols(psdin, 0, 1);
	dmat* psd1=dsub(psdin, 0, 0, 1, 1);//copy
	dmat* f2=dref(fnew);
	real t1=dtrapz(f1, psd1);
	real ydefault=1e-40;
	if(uselog){
		dcwlog(psd1);
		ydefault=log(ydefault);
	}
	dmat* psd2=dinterp1(f1, psd1, f2, ydefault);
	if(uselog){
		dcwexp(psd2, 1);
	}
	real t2=dtrapz(f2, psd2);
	if(fabs(t1-t2)>fabs(0.5*(t1+t2)*2)){
		warning("psd interpolation failed. int_orig=%g, int_interp=%g\n", t1, t2);
	}
	//Don't scale psd2 as it may have overlapping frequency regions
	dfree(f1); dfree(f2); dfree(psd1);
	return psd2;
}
/**
   Find vibration peaks in the PSD by comparing the PSD against a LPF version plus noise.

   @param[in] psdin The PSD with frequency in the first column and PSD in the second column.
   @return Newly allocated 4*npeak array with each column containing the centroid frequency, integrated power, starting index, and ending index of a peak; the caller frees it.
*/
dmat* psd_vibid(const dmat* psdin){
	int nmaxp=100;
	if(!psdin||NY(psdin)<2||NX(psdin)<3){
		warning("psdin should have at least 3 rows and 2 columns.\n");
		return dnew(4, 0);
	}
	real* f=P(psdin);
	real* psd=PCOL(psdin, 1);
	dmat* y=dsub(psdin, 0, 0, 1, 1);
	const real gain=0.1;
	const real gain2=0.05;
	int inpeak=0;
	real ylpf0=P(y,1);
	real dylpf0=fabs(P(y,1)-P(y,0));
	real ylpf=ylpf0, dylpf=dylpf0;//start the second stage from the first stage.
	dmat* res=dnew(4, nmaxp);
	real thres=25e-18;/*threshold: 5 nm*/
	real sumxy=0, sumy=0, sum=0;
	int count=0;
	for(long i=1; i<psdin->nx-1; i++){
		if(!inpeak){
			//second order LPF
			ylpf0=(1.-gain)*ylpf0+P(y,i)*gain;
			ylpf=(1.-gain)*ylpf+ylpf0*gain;
			real diff=P(y,i)-P(y,i-1);
			if(diff>0){
				dylpf0=(1.-gain2)*dylpf0+diff*gain2;
				dylpf=(1.-gain2)*dylpf+dylpf0*gain2;
			}
			if(P(y,i+1)>ylpf+dylpf*5&&f[i]>1){//beginning of peak
				inpeak=1;
				if(count>0&&f[i]<f[(int)P(res, 3, count-1)]+0.1){
					//combine with last peak if within 1 Hz.
					count--;
				} else{
					P(res, 2, count)=i;
					sumxy=f[i]*psd[i];//for CoG
					sumy=psd[i];//for CoG
					sum=0;//integration
				}
			}
		} else{
			//continuation of peak
			sumxy+=f[i]*psd[i];
			sumy+=psd[i];
			sum+=(f[i]-f[i-1])*(psd[i]+psd[i-1]);
			if(P(y,i)<ylpf+dylpf&&P(y,i+1)<ylpf+dylpf){//end of peak
				inpeak=0;
				if(sum*0.5>thres){
					P(res, 0, count)=sumxy/sumy;
					P(res, 1, count)=sum*0.5;
					P(res, 3, count)=i;
					count++;
					if(count==nmaxp){
						nmaxp*=2;
						dresize(res, 4, nmaxp);
					}
				}
			}
		}
	}
	if(inpeak&&sum*0.5>thres){//peak extends to the end. Also keeps a peak that has been combined.
		P(res, 0, count)=sumxy/sumy;
		P(res, 1, count)=sum*0.5;
		P(res, 3, count)=NX(psdin)-1;
		count++;
	}
	dfree(y);
	dresize(res, 4, count);
	return res;
}

/**
   Convert temporal PSD to spatial

   @param[in] psdt The temporal PSD as a row vector [alpha, beta, fmin, fmax] of power law beta*f^alpha. fmin and fmax are optional.
   @param[in] vmean The mean wind speed (m/s).
   @return Newly allocated spatial PSD in the same format; the caller frees it.
*/
dmat* psd_t2s(const dmat* psdt, real vmean){
	if(NX(psdt)!=1||NY(psdt)>4){
		error("psdt needs to be 1 row and less than 4 cols\n");
		return NULL;
	}
	real alpha=P(psdt,0);//power
	real beta=P(psdt,1);//strength
	real alpha2=alpha-1;
	//n is -alpha2
	real bfun=tgamma(0.5)*tgamma((-alpha2-1)*0.5)/tgamma(-alpha2*0.5);
	real beta2=beta/bfun*pow(vmean, 2+alpha2);
	dmat* psds=dnew(NX(psdt), NY(psdt));
	P(psds,0)=alpha2;
	P(psds,1)=beta2;
	if(NY(psds)>2){
		P(psds,2)=P(psdt,2)/vmean;
	}
	if(NY(psds)>3){
		P(psds,3)=P(psdt,3)/vmean;
	}
	return psds;
}

/**
   Convert spatial PSD to temporal

   @param[in] psds The spatial PSD as a row vector [alpha, beta, fmin, fmax] of power law beta*f^alpha. fmin and fmax are optional.
   @param[in] vmean The mean wind speed (m/s).
   @return Newly allocated temporal PSD in the same format; the caller frees it.
*/
dmat* psd_s2t(const dmat* psds, real vmean){
	if(NX(psds)!=1||NY(psds)>4){
		error("psds needs to be 1 row and less than 4 cols\n");
		return NULL;
	}
	real alpha=P(psds,0);//power
	real beta=P(psds,1);//strength
	real alpha2=alpha+1;
	//n is -alpha
	real bfun=tgamma(0.5)*tgamma((-alpha-1)*0.5)/tgamma(-alpha*0.5);
	real beta2=beta*bfun*pow(vmean, -2-alpha);
	dmat* psdt=dnew(NX(psds), NY(psds));
	P(psdt,0)=alpha2;
	P(psdt,1)=beta2;
	if(NY(psds)>2){
		P(psdt,2)=P(psds,2)*vmean;
	}
	if(NY(psds)>3){
		P(psdt,3)=P(psds,3)*vmean;
	}
	return psdt;
}


/**
   Integrated a PSD that defines on linear or logrithmically spaced grid nu.

   @param[in] nu The frequency grid.
   @param[in] psd The PSD.
   @param[in] n Number of points.
   @return The integral, or 0 if the grid is neither linear nor logarithmic.
*/
real psd_inte(const real* nu, const real* psd, long n){
	real dnu=(nu[n-1]-nu[0])/(n-1);
	real dlognu=(log(nu[n-1])-log(nu[0]))/(n-1);
	real res_sig=0;
	if(fabs(nu[1]-nu[0]-dnu)<dnu*1.e-4){
		for(long i=0; i<n; i++){
			res_sig+=psd[i];
		}
		res_sig*=dnu;
	} else if(fabs(log(nu[1])-log(nu[0])-dlognu)<dlognu*1.e-4){
		for(long i=0; i<n; i++){
			res_sig+=psd[i]*nu[i];
		}
		res_sig*=dlognu;
	}
	return res_sig;
}
/**
   wraps psd_inte with a single input variable containing both frequency and psd.

   @param[in] psdin The PSD with frequency in the first column and PSD in the second column.
   @return The integral of the PSD.
*/
real psd_inte2(const dmat* psdin){
	if(NY(psdin)!=2){
		error("psdin should have two columns\n");
		return NAN;
	}
	real* nu=P(psdin);
	long n=NX(psdin);
	real* psd=nu+n;
	return psd_inte(nu, psd, n);
}
/**
   Reverse cumulative sum of PSD

   @param[in] psdin The input PSD with frequency in the first column and PSD in the remaining columns.
   @param[in] scale Scaling applied to the result.
   @return Newly allocated array with frequency in the first column and sqrt of the reverse integration of PSD in the remaining; the caller frees it.
*/
dmat *psd_reverse_cumu(const dmat *psdin, real scale){
	long nx=NX(psdin);
	dmat *cumu=dnew(NX(psdin), NY(psdin));
	const real *nu=P(psdin);
	real dlogf=(log(nu[nx-1])-log(nu[0]))/(nx-1);	
	if(fabs(log(nu[1])-log(nu[0])-dlogf)<dlogf*1.e-4){//log spacing
		info("log spacing\n");
		for(int iy=1; iy<NY(psdin); iy++){
			long ix=nx-1;
			P(cumu, ix, iy)=P(psdin, ix, iy)*P(psdin, ix, 0)*dlogf;
			for(ix=nx-2; ix>-1; ix--){
				P(cumu, ix, iy)=P(cumu, ix+1, iy)+P(psdin, ix, iy)*P(psdin, ix, 0)*dlogf;
			}
		}
	} else{
		for(int iy=1; iy<NY(psdin); iy++){
			long ix=nx-1;
			P(cumu, ix, iy)=P(psdin, ix, iy)*(P(psdin, ix, 0)-P(psdin, ix-1, 0));
			for(ix=nx-2; ix>-1; ix--){
				P(cumu, ix, iy)=P(cumu, ix+1, iy)+P(psdin, ix, iy)*(P(psdin, ix+1, 0)-P(psdin, ix, 0));
			}
		}
	}
	for(int iy=1; iy<NY(psdin); iy++){
		for(int ix=0; ix<NX(cumu); ix++){
			P(cumu, ix, iy)=sqrt(P(cumu, ix, iy))*scale;
		}
	}
	for(int ix=0; ix<NX(cumu); ix++){
		P(cumu, ix, 0)=P(psdin, ix, 0);
	}
	//real inte=psd_inte2(psdin);
	//info("psd integral is %g, cum[0] is %g\n", sqrt(inte)*1e9, P(cumu, 0, 1));
	return cumu;
}

/**
   Convert PSD into time series. wraps psd2ts with a seed instead of rand_t as input.

   @param[in] psdin The PSD. See psd2ts().
   @param[in] seed The random seed.
   @param[in] dt The sampling period (s).
   @param[in] nstepin Number of time steps.
   @return Newly allocated time series of length nstepin; the caller frees it.
*/
dmat *psd2ts2(const dmat *psdin, int seed, real dt, int nstepin){
	rand_t stat0;
	seed_rand(&stat0, seed);
	return psd2ts(psdin, &stat0, dt, nstepin);
}
/**
   Convert PSD into time series.

   @param[in] psdin The PSD. Either a row vector [alpha, beta, fmin, fmax] of power law beta*f^alpha (fmin and fmax are optional), or two columns of frequency and PSD.
   @param[in,out] rstat The random number generator state.
   @param[in] dt The sampling period (s).
   @param[in] nstepin Number of time steps.
   @return Newly allocated time series of length nstepin; the caller frees it.
*/
dmat* psd2ts(const dmat* psdin, rand_t* rstat, real dt, int nstepin){
	if(!psdin){
		error("psdin cannot be null\n");
		return NULL;
	}
	long nstep=nextpow2(nstepin);
	real df=1./(dt*nstep);
	dmat* fs=dlinspace(0, df, nstep);
	dmat* psd=NULL;
	if(NY(psdin)==1){//[alpha, beta, fmin, fmax] discribes power law with cut on/off freq.
		psd=dnew(nstep, 1);
		real alpha=P(psdin,0);
		real beta=P(psdin,1);
		long i0=1, imax=nstep;
		if(NX(psdin)>2){
			i0=(long)round(P(psdin,2)/df);
			if(i0<1) i0=1;
		}
		if(NX(psdin)>3){
			imax=(long)round(P(psdin,3)/df);
			if(imax>nstep) imax=nstep;//psd has nstep elements
		}
		dbg("fmin=%g, fmax=%g, df=%g, i0=%ld, imax=%ld\n",
			P(psdin,2), P(psdin,3), df, i0, imax);
		for(long i=i0; i<imax; i++){
			P(psd,i)=beta*pow(i*df, alpha);
		}
	} else if(NY(psdin)==2){
		if(NX(psdin)<2){
			error("Invalid PSD. nx should be more than 1\n");
			dfree(fs);
			return NULL;
		}
		psd=dinterp1(psdin, 0, fs, 1e-40);
		P(psd,0)=0;/*disable pistion. */
	} else{
		error("psdin is invalid format. shape is %ldx%ld\n", NX(psdin), NY(psdin));
		dfree(fs);
		return NULL;
	}
	cmat* wshat=cnew(nstep, 1);
	//cfft2plan(wshat, -1);
	for(long i=0; i<nstep; i++){
		P(wshat,i)=sqrt(P(psd,i)*df)*COMPLEX(randn(rstat), randn(rstat));
	}
	cfft2(wshat, -1);
	dmat* out=NULL;
	creal2d(&out, 0, wshat, 1);
	cfree(wshat);
	dfree(psd);
	dfree(fs);
	dresize(out, nstepin, 1);
	return out;
}

/**
   Check whether two PSDs has the same frequency. the first column of each
   dmat is the frequency nu, and the second column is PSD.

   @param[in] psd1 The first PSD.
   @param[in] psd2 The second PSD.
   @return 1 if the frequencies match, 0 otherwise.
*/
static int check_psd_match(const dmat* psd1, const dmat* psd2){
	real ans=1;
	if(NX(psd1)!=NX(psd2)){
		ans=0;
	} else{
		for(long i=0; i<NX(psd1); i++){
			if(fabs(P(psd1,i)-P(psd2,i))>fabs(P(psd1,i)+P(psd2,i))*1.e-6){
				ans=0;
				break;
			}
		}
	}
	return ans;
}
/**
   Add two PSDs. The first column of each dmat is the frequency nu, and the
   second column is PSD

   @param[in] psd1 The first PSD.
   @param[in] psd2 The second PSD, interpolated onto the frequency of psd1 if needed.
   @param[in] scale2 Scaling of psd2.
   @return Newly allocated psd1+psd2*scale2; the caller frees it.
*/
dmat* add_psd(const dmat* psd1, const dmat* psd2, real scale2){
	dmat* out=ddup(psd1);
	add_psd2(&out, psd2, scale2);
	return out;
}

/**
   Add a PSD scaled by scale to another. The first column of each dmat is the
   frequency nu, and the second column is PSD.

   @param[in,out] pout The PSD to add to. If NULL, it is set to a copy of in times scale.
   @param[in] in The PSD to add, interpolated onto the frequency of pout if needed.
   @param[in] scale Scaling of in.
*/
void add_psd2(dmat** pout, const dmat* in, real scale){
	if(!*pout){
		*pout=ddup(in);
		if(scale!=1){
			for(long i=0; i<NX(in); i++){
				P(*pout, i, 1)*=scale;
			}
		}
	} else{
		dmat* out=*pout;
		real* p1=PCOL(out, 1);
		dmat* p2new=0;
		const long nx=NX(out);
		const real* p2=0;
		if(check_psd_match(out, in)){
			p2=PCOL(in, 1);
		} else{
			dmat* nu1=dsub(out, 0, nx, 0, 1);
			p2new=dinterp1(in, 0, nu1, 1e-40);
			p2=PCOL(p2new, 0);
			dfree(nu1);
		}

		for(long i=0; i<nx; i++){
			p1[i]+=p2[i]*scale;
		}
		dfree(p2new);
	}
}
/**
   Choose and sum modes from the PSD. The first column is frequency.

   @param[in] psd The input PSD with frequency in the first column.
   @param[in] im Index of the first mode, 0 to ny(psd)-2. -1: sum all modes.
   @param[in] nm Number of modes. At least 1.
   @param[in] keepdc Do not remove DC component (f=0).
   @param[in] scale Scaling factor.
   @return Newly allocated two column array of frequency and summed PSD; the caller frees it.
*/
dmat* psd_select(const dmat* psd, int im, int nm, int keepdc, real scale){
	int skip=(keepdc||P(psd,0)>0)?0:1;
	dmat *psd2=dnew(psd->nx-skip, 2);
	int jm=0;
	if(im<0){//all modes
		im=1;
		jm=(int)NY(psd);
	}else{
		im++;
		if(nm<=0){
			nm=1;
		}
		jm=(int)MIN(im+nm, NY(psd));
	}
	for(int ix=0; ix<NX(psd2); ix++){
		real tmp=0;
		for(int iy=im; iy<jm; iy++){
			tmp+=P(psd, ix+skip, iy);
		}
		P(psd2, ix, 0)=P(psd, ix+skip, 0);
		P(psd2, ix, 1)=tmp*scale;
	}
	return psd2;
}
/**
   from a 2-d screen, compute 1-d PSD by radially averaging the 2-D psd assuming isotropy.

   @param[in] screen The square 2-d screen.
   @param[in] dx The spatial sampling (m).
   @return Newly allocated two column array of spatial frequency (1/m) and PSD; the caller frees it.
*/
dmat *psd2d_aniso(const dmat *screen, real dx){
	if(!screen || NX(screen)!=NY(screen)){
		error("psd2d_aniso: expects a square screen\n");
		return NULL;
	}
	cmat *hat=NULL;
	ccpd(&hat, screen);
	cfft2(hat, -1);
	dmat *psd=NULL;
	cfftshift(hat);
	cabs22d(&psd, 1, hat, 1./(PN(screen)*PN(screen)));
	cfree(hat);
	long npsd=MIN(NX(screen), NY(screen))/2+1;
	dmat *rvec=dnew(npsd,1);
	for(long i=0; i<npsd; i++){
		P(rvec,i)=i;
	}
	//dbg("npsd=%ld\n", npsd);
	dmat *psd1d=denc(psd, rvec, -1, NTHREAD);
	dfree(rvec);
	dfree(psd);
	dresize(psd1d, NX(psd1d), 2);
	real df=1./(NX(screen)*dx);
	real df1=1./df;
	for(long i=0; i<npsd; i++){
		P(psd1d, i,1)=P(psd1d,i,0)*df1;
		P(psd1d, i,0)=df*i;
	}
	return psd1d;
}
/**
   From a 2-d screen, compute 1-d PSD by averaging the 1-D psd of each
   row/column, and then scale back to the 2-D psd using the technique of computing
   temporay PSD of turbulence.

   @param[out] extra Extra output [r0 equivalent at 0.5 micron (m), power law slope]. Can be NULL.
   @param[in] screens_ 2-d square screen or cell of 2-d screens of the same size.
   @param[in] dx Sampling (m).
   @return Newly allocated two column array of spatial frequency (1/m) and PSD, or NULL if screens_ is empty; the caller frees it.
*/
dmat *psd2d(dmat **extra, const_anyarray screens_, real dx){
	const cell* screens=screens_.c;
	if(!screens){
		dbg("screen is empty\n");
		return NULL;
	}
	long nc=iscell(screens)?PN(screens):1;
	dmat *psd2=NULL;
	real df=0;
	for(long ic=0; ic<nc; ic++){
		dmat *screen=dmat_cast(iscell(screens)?P(screens,ic):screens);
		if(NX(screen)!=NY(screen)){
			error("psd2d: expects a square screen\n");
			return NULL;
		}
		dmat *p1=psd1d(screen, 1);
		dmat *screen2=dtrans(screen);
		dmat *p2=psd1d(screen2, 1);
		dfree(screen2);
		if(!psd2){
			psd2=dnew(NX(p1),2);
		}else if(NX(psd2)!=NX(p1)){
			error("psd2d: screens shall have the same size\n");
			continue;
		}
		if(ic==0){
			df=1/(NX(screen)*dx);
			for(long ix=0; ix<NX(psd2); ix++){
				P(psd2, ix, 0)=df*ix;
			}
		}
		real scale=1./(NY(p1)+NY(p1))/df/nc;
		for(long ix=0; ix<NX(psd2); ix++){
			real tmp=0;
			for(long iy=0; iy<NY(p1); iy++){
				tmp+=P(p1, ix, iy)+P(p2,ix,iy);
			}
			P(psd2, ix, 1)+=tmp*scale/P(psd2, ix, 0);
		}
		dfree(p1);
		dfree(p2);
	}
	P(psd2, 0, 1)=P(psd2, 1, 1);//remove the singular value.
	
	/*
		1-D psd is converted to 2-D psd. B(1/2, (n-1)/2) is accounted for.
	*/
	//Fit power law, then determine the coefficient of scaling
	dmat *x=dnew(NX(psd2)-1,1);//skip DC term
	dmat *y=dnew(NX(psd2)-1,1);
	for(long i=0; i<NX(x); i++){
		P(x, i)=log(P(psd2, i+1, 0));
		P(y, i)=log(P(psd2, i+1, 1));
	}
	dmat *coeff=polyfit(x,y,1);
	dfree(x);
	dfree(y);
	real ss=(-P(coeff,0)-1)*0.5;
	//beta(a,b)=gamma(a)gamma(b)/gamma(a+b)
	real bt=1./(tgamma(0.5)*tgamma(ss)/tgamma(0.5+ss));//correction factor
	//dbg("PSD scaled by %g\n", bt);
	for(long i=0; i<NX(psd2); i++){
		P(psd2, i, 1)*=bt;
	}
	if(extra){
		*extra=dnew(2,1);//r0, slope
		real ckp2=0.0229*pow(2.*M_PI/0.5e-6,-2);
		P(*extra,0)=pow(exp(P(coeff,1))*bt/ckp2, -3./5.);//r0 equivalent
		P(*extra,1)=P(coeff,0); //slope
	}
	dfree(coeff);
	return psd2;
}

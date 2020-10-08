/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>

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
#ifndef FRACTAL
#error "do not use fractal_do.c along."
#endif

/**
   Define the following:
   FRACTAL: function name
   INVERSE: invert or now. 1: inverse, 0: not
   TRANPOSE: transpose the operation or not.
   F(r): The structure function
*/

/*derived parameter: FORWARD: direction, 1: forward, 0: backward */
#if (INVERSE==0 && TRANSPOSE ==0) || (INVERSE == 1 && TRANSPOSE == 1)
#define FORWARD 1
#else
#define FORWARD 0
#endif


void FRACTAL(dmat* p0, real dx, real r0, real L0, long ninit){
	const long nx=p0->nx;
	const long ny=p0->ny;
	assert(nx==ny);
	LOCK(mutex_cov);
	long step0=(nx-1)/(ninit-1);
#ifndef NDEBUG
	if(ninit<1||(ninit-1)*step0!=(nx-1)||(step0&(step0-1))!=0){
		dbg("nx=%ld, ninit=%ld, step0=%ld\n", nx, ninit, step0);
	}
#endif
	assert(ninit>1);
	assert((ninit-1)*step0==(nx-1));/*must evenly divide. */
	assert((step0&(step0-1))==0);/*step0 must be power of 2. */
	vkcov_t* node=vkcov_calc(r0, L0, dx, nx, ninit);
	dmat* cov=node->cov;
	UNLOCK(mutex_cov);
	dmat* pcov=cov;
	const long nx1=nx-1;
	const long ny1=ny-1;
	const long norder=mylog2(step0);
	const real c0=P(pcov, 0, 0);
#if FORWARD == 1
	{
		dmat* pi=dnew(ninit, ninit);
		for(long iy=0; iy<ninit;iy++){
			for(long ix=0; ix<ninit; ix++){
				P(pi, ix, iy)=P(p0, ix*step0, iy*step0);
			}
		}
		/*reshape pi; */
		pi->nx=ninit*ninit;
		pi->ny=1;
		dmat* qi=dnew(ninit*ninit, 1);
#if INVERSE ==0
		dmm(&qi, 0, node->K, pi, "nn", 1);
#else
		dmm(&qi, 0, node->KI, pi, "tn", 1);
#endif
	/*reshape to square. */
		qi->nx=ninit;
		qi->ny=ninit;
		dmat* pqi=qi;
		for(long iy=0; iy<ninit;iy++){
			for(long ix=0; ix<ninit; ix++){
				P(p0, ix*step0, iy*step0)=P(pqi, ix, iy);
			}
		}
		dfree(pi);
		dfree(qi);
	}
#endif
#if TRANSPOSE == 0
#if INVERSE==0
#define QUA(p0,p1,p2,p3,p4) p0=qua1*(p1+p2+p3+p4)+qua0*p0
#define DIA(p0,p1,p2,p3,p4) p0=dia1*(p1+p2+p3+p4)+dia0*p0
#define TRI(p0,p1,p2,p3)    p0=(p1+p2)*tri1+p3*tri3+tri0*p0
#else /*INVERSE==0 */
#define QUA(p0,p1,p2,p3,p4) p0=(p0-qua1*(p1+p2+p3+p4))/qua0;
#define DIA(p0,p1,p2,p3,p4) p0=(p0-dia1*(p1+p2+p3+p4))/dia0;
#define TRI(p0,p1,p2,p3)    p0=(p0-((p1+p2)*tri1+p3*tri3))/tri0;
#endif
#else /*TRANPOSE == 0 */
	real tmp;
#if INVERSE==0
#define QUA(p0,p1,p2,p3,p4) tmp=qua1*p0; p1+=tmp; p2+=tmp; p3+=tmp; p4+=tmp; p0=qua0*p0;
#define DIA(p0,p1,p2,p3,p4) tmp=dia1*p0; p1+=tmp; p2+=tmp; p3+=tmp; p4+=tmp; p0=dia0*p0;
#define TRI(p0,p1,p2,p3) p1+=p0*tri1; p2+=p0*tri1; p3+=p0*tri3; p0*=tri0;
#else /*INVERSE==0 */
#define QUA(p0,p1,p2,p3,p4) tmp=-qua1/qua0*p0; p1+=tmp; p2+=tmp; p3+=tmp; p4+=tmp; p0=p0/qua0;
#define DIA(p0,p1,p2,p3,p4) tmp=-dia1/dia0*p0; p1+=tmp; p2+=tmp; p3+=tmp; p4+=tmp; p0=p0/dia0;
#define TRI(p0,p1,p2,p3) p1-=tri1/tri0*p0; p2-=tri1/tri0*p0; p3-=tri3/tri0*p0; p0/=tri0;
#endif
#endif
#if FORWARD == 1 /*Define LOOP to avoid wrong indentation */
#define LOOP long order=norder; order>0; order--
#else
#define LOOP long order=1; order<=norder; order++
#endif
	for(LOOP){
#undef LOOP
		long step=1<<order;/*step of the parent grid. */
		real c1=P(pcov, 0, order);
		real c2=P(pcov, 1, order);
		real c3=P(pcov, 0, order+1);
		real c4=P(pcov, 1, order+1);
		/*for case 1: square case */
		real qua1=c2/(c0+2*c3+c4);
		real qua0=sqrt(c0-4.*c2*qua1);
		/*for case 3: Diamond configuration */
		real dia1=c1/(c0+2*c2+c3);
		real dia0=sqrt(c0-4.*c1*dia1);
		/*for case 2: triangular case */
		real trii=1./(c0*(c0+c3)-2*c2*c2);
		real tri1=c1*(c0-c2)*trii;
		real tri3=c1*(c0-2*c2+c3)*trii;
		real tri0=sqrt(c0-c1*c1*(3*c0-4*c2+c3)*trii);
		long ny2=ny1-step;
		long nx2=nx1-step;
		long step2=step>>1;
#if FORWARD == 1
	/*finish all the black ones before we do white ones, which depend on this */
		for(long offy=0; offy<ny1; offy+=step){
			for(long offx=0; offx<nx1; offx+=step){
			/*case 1: square configuration */
				QUA(P(p0, offx+step2, offy+step2),
					P(p0, offx, offy),
					P(p0, offx, offy+step),
					P(p0, offx+step, offy+step),
					P(p0, offx+step, offy));
			}
		}
#endif
	/*do all the white ones */
		for(long offy=0; offy<ny1; offy+=step){
			for(long offx=0; offx<nx1; offx+=step){
				if(offx==0){/*do the top one, triangular  */
					TRI(P(p0, offx, offy+step2),
						P(p0, offx, offy),
						P(p0, offx, offy+step),
						P(p0, offx+step2, offy+step2));
				}

				/*bottom one.  */
				if(offx!=nx2){/*diamond */
					DIA(P(p0, offx+step, offy+step2),
						P(p0, offx+step, offy),
						P(p0, offx+step2, offy+step2),
						P(p0, offx+step+step2, offy+step2),
						P(p0, offx+step, offy+step));
				} else{/*do triangular case */
					TRI(P(p0, offx+step, offy+step2),
						P(p0, offx+step, offy),
						P(p0, offx+step, offy+step),
						P(p0, offx+step2, offy+step2));

				}

				if(offy==0){/*do the left one, triangular  */
					TRI(P(p0, offx+step2, offy),
						P(p0, offx, offy),
						P(p0, offx+step, offy),
						P(p0, offx+step2, offy+step2));
				}

				/*right one. */
				if(offy!=ny2){/*diamond */
					DIA(P(p0, offx+step2, offy+step),
						P(p0, offx+step2, offy+step2),
						P(p0, offx, offy+step),
						P(p0, offx+step, offy+step),
						P(p0, offx+step2, offy+step+step2));
				} else{/*do triangular case */
					TRI(P(p0, offx+step2, offy+step),
						P(p0, offx, offy+step),
						P(p0, offx+step, offy+step),
						P(p0, offx+step2, offy+step2));
				}
			}
		}
#if FORWARD == 0
	/*finish all the black ones before we do white ones, which depend on this */
		for(long offy=0; offy<ny1; offy+=step){
			for(long offx=0; offx<nx1; offx+=step){
			/*case 1: square configuration */
				QUA(P(p0, offx+step2, offy+step2),
					P(p0, offx, offy),
					P(p0, offx, offy+step),
					P(p0, offx+step, offy+step),
					P(p0, offx+step, offy));
			}
		}
#endif
	}
#if FORWARD == 0
	{
		dmat* pi=dnew(ninit, ninit);
		dmat* ppi=pi;
		for(long iy=0; iy<ninit;iy++){
			for(long ix=0; ix<ninit; ix++){
				P(ppi, ix, iy)=P(p0, ix*step0, iy*step0);
			}
		}
		/*reshape pi; */
		pi->nx=ninit*ninit;
		pi->ny=1;
		dmat* qi=dnew(ninit*ninit, 1);
#if INVERSE ==0
		dmm(&qi, 0, node->K, pi, "tn", 1);
#else
		dmm(&qi, 0, node->KI, pi, "nn", 1);
#endif
	/*reshape to square. */
		qi->nx=ninit;
		qi->ny=ninit;
		dmat* pqi=qi;
		for(long iy=0; iy<ninit;iy++){
			for(long ix=0; ix<ninit; ix++){
				P(p0, ix*step0, iy*step0)=P(pqi, ix, iy);
			}
		}
		dfree(pi);
		dfree(qi);
	}
#endif
#undef QUA
#undef DIA
#undef TRI
}
#undef FORWARD

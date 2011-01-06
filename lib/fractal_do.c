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
//derived parameter: FORWARD: direction, 1: forward, 0: backward
#if (INVERSE==0 && TRANSPOSE ==0) || (INVERSE == 1 && TRANSPOSE == 1)
#define FORWARD 1
#else
#define FORWARD 0
#endif


void FRACTAL(double *p0, long nx, long ny, double dx, double r0, double L0){
    if(((nx-1) & (nx-2)) != 0  || ((ny-1) & (ny-2)) !=0 || nx != ny){
	error("nx=%ld, ny=%ld: they need to be 1+power of 2, and equal\n", nx, ny);
    }
    LOCK(mutex_cov);
    dmat *cov=vkcov_calc(r0, L0, dx, nx);
    UNLOCK(mutex_cov);
    PDMAT(cov, pcov);
    const long nx1=nx-1;
    const long ny1=ny-1;
    const long norder=cov->ny-2;
    const double c0=pcov[0][0];
    double (*p)[nx]=(void*)p0;
#if FORWARD == 1
    {
	//First generate four outmost values.
	double c1=pcov[norder+1][0];
	double c2=pcov[norder+1][1];
#if INVERSE == 0
	double a=sqrt(c0+2*c1+c2);
	double b=sqrt(c0-2*c1+c2);
	double c=sqrt(2*(c0-c2));
#else //INVERSE
	double a=1./sqrt(c0+2*c1+c2);
	double b=1./sqrt(c0-2*c1+c2);
	double c=2./sqrt(2*(c0-c2));
#endif
	assert(c0+2*c1+c2>=0 && c0-2*c1+c2>=0 && c0-c2>=0);
	double *p1=&p[0][0];
	double *p2=&p[0][nx1];
	double *p3=&p[ny1][nx1];
	double *p4=&p[ny1][0];

	double q1=0.5*(a**p1-b**p2-c**p3      );
	double q2=0.5*(a**p1+b**p2      -c**p4);
	double q3=0.5*(a**p1-b**p2+c**p3      );
	double q4=0.5*(a**p1+b**p2      +c**p4);

	*p1=q1;
	*p2=q2;
	*p3=q3;
	*p4=q4;
    }
#endif
#if TRANSPOSE == 0
#if INVERSE==0
#define QUA(p0,p1,p2,p3,p4) p0=qua1*(p1+p2+p3+p4)+qua0*p0
#define DIA(p0,p1,p2,p3,p4) p0=dia1*(p1+p2+p3+p4)+dia0*p0
#define TRI(p0,p1,p2,p3)    p0=(p1+p2)*tri1+p3*tri3+tri0*p0
#else //INVERSE==0
#define QUA(p0,p1,p2,p3,p4) p0=(p0-qua1*(p1+p2+p3+p4))/qua0;
#define DIA(p0,p1,p2,p3,p4) p0=(p0-dia1*(p1+p2+p3+p4))/dia0;
#define TRI(p0,p1,p2,p3)    p0=(p0-((p1+p2)*tri1+p3*tri3))/tri0;
#endif
#else //TRANPOSE == 0
    double tmp;
#if INVERSE==0
#define QUA(p0,p1,p2,p3,p4) tmp=qua1*p0; p1+=tmp; p2+=tmp; p3+=tmp; p4+=tmp; p0=qua0*p0;
#define DIA(p0,p1,p2,p3,p4) tmp=dia1*p0; p1+=tmp; p2+=tmp; p3+=tmp; p4+=tmp; p0=dia0*p0;
#define TRI(p0,p1,p2,p3) p1+=p0*tri1; p2+=p0*tri1; p3+=p0*tri3; p0*=tri0;
#else //INVERSE==0
#define QUA(p0,p1,p2,p3,p4) tmp=-qua1/qua0*p0; p1+=tmp; p2+=tmp; p3+=tmp; p4+=tmp; p0=p0/qua0;
#define DIA(p0,p1,p2,p3,p4) tmp=-dia1/dia0*p0; p1+=tmp; p2+=tmp; p3+=tmp; p4+=tmp; p0=p0/dia0;
#define TRI(p0,p1,p2,p3) p1-=tri1/tri0*p0; p2-=tri1/tri0*p0; p3-=tri3/tri0*p0; p0/=tri0;
#endif
#endif
#if FORWARD == 1 //Define LOOP to avoid wrong indentation
#define LOOP long order=norder; order>0; order--
#else
#define LOOP long order=1; order<=norder; order++
#endif
    for(LOOP){
#undef LOOP
	long step=1<<order;
	double c1=pcov[order][0];
	double c2=pcov[order][1];
	double c3=pcov[order+1][0];
	double c4=pcov[order+1][1];
	//for case 1: square case
	double qua1=c2/(c0+2*c3+c4);
	double qua0=sqrt(c0-4.*c2*qua1);
	//for case 3: Diamond configuration
	double dia1=c1/(c0+2*c2+c3);
	double dia0=sqrt(c0-4.*c1*dia1);
	//for case 2: triangular case
	double trii=1./(c0*(c0+c3)-2*c2*c2);
	double tri1=c1*(c0-c2)*trii;
	double tri3=c1*(c0-2*c2+c3)*trii;
	double tri0=sqrt(c0-c1*c1*(3*c0-4*c2+c3)*trii);
	long ny2=ny1-step;
	long nx2=nx1-step;
	long step2=step>>1;
#if FORWARD == 1
	//finish all the black ones before we do white ones, which depend on this
	for(long offy=0; offy<ny1; offy+=step){
	    for(long offx=0; offx<nx1; offx+=step){
		//case 1: square configuration
		QUA(p[offy+step2][offx+step2],
		    p[offy][offx],
		    p[offy+step][offx],
		    p[offy+step][offx+step],
		    p[offy][offx+step]);
	    }
	}
#endif
	//do all the white ones
	for(long offy=0; offy<ny1; offy+=step){
	    for(long offx=0; offx<nx1; offx+=step){
		if(offx==0){//do the top one, triangular 
		    TRI(p[offy+step2][offx],
			p[offy][offx],
			p[offy+step][offx],
			p[offy+step2][offx+step2]);
		}
		
		//bottom one. 
		if(offx!=nx2){//diamond
		    DIA(p[offy+step2][offx+step],
			p[offy][offx+step],
			p[offy+step2][offx+step2],
			p[offy+step2][offx+step+step2],
			p[offy+step][offx+step]);
		}else{//do triangular case
		    TRI(p[offy+step2][offx+step],
			p[offy][offx+step],
			p[offy+step][offx+step],
			p[offy+step2][offx+step2]);

		}

		if(offy==0){//do the left one, triangular 
		    TRI(p[offy][offx+step2],
			p[offy][offx],
			p[offy][offx+step],
			p[offy+step2][offx+step2]);
		}

		//right one.
		if(offy!=ny2){//diamond
		    DIA(p[offy+step][offx+step2],
			p[offy+step2][offx+step2],
			p[offy+step][offx],
			p[offy+step][offx+step],
			p[offy+step+step2][offx+step2]);
		}else{//do triangular case
		    TRI(p[offy+step][offx+step2],
			p[offy+step][offx],
			p[offy+step][offx+step],
			p[offy+step2][offx+step2]);
		}
	    }
	}
#if FORWARD == 0
	//finish all the black ones before we do white ones, which depend on this
	for(long offy=0; offy<ny1; offy+=step){
	    for(long offx=0; offx<nx1; offx+=step){
		//case 1: square configuration
		QUA(p[offy+step2][offx+step2],
		    p[offy][offx],
		    p[offy+step][offx],
		    p[offy+step][offx+step],
		    p[offy][offx+step]);
	    }
	}
#endif
    }
#if FORWARD == 0
    {
	//generate four outmost values. 
	double c1=pcov[norder+1][0];
	double c2=pcov[norder+1][1];
#if INVERSE == 0
	double a=sqrt(c0+2*c1+c2);
	double b=sqrt(c0-2*c1+c2);
	double c=sqrt(2*(c0-c2));
#else //INVERSE
	double a=1./sqrt(c0+2*c1+c2);
	double b=1./sqrt(c0-2*c1+c2);
	double c=2./sqrt(2*(c0-c2));
#endif
	assert(a>0 && b>0 && b>0);
	double *p1=&p[0][0];
	double *p2=&p[0][nx1];
	double *p3=&p[ny1][nx1];
	double *p4=&p[ny1][0];

	double q1=0.5*( a**p1+a**p2+a**p3+a**p4);
	double q2=0.5*(-b**p1+b**p2-b**p3+b**p4);
	double q3=0.5*(-c**p1      +c**p3      );
	double q4=0.5*(      -c**p2      +c**p4);
	*p1=q1;
	*p2=q2;
	*p3=q3;
	*p4=q4;
    }
#endif
#undef QUA
#undef DIA
#undef TRI
}
#undef FORWARD

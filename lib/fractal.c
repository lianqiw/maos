#include "fractal.h"

/**
   \file fractal.c

   Implementation of the fractal operation for atmospheric turbulence screen
   generation and reconstruction.

*/


#define F(r) coeff*pow(r, power)/*<Compute structure function of separation
				     r, and Fried parameter of r0.*/
/**
   Apply forward fractal operator recursively.

   To confirm that fractal_trans is right, compare u'Av with v'A'u, where ' is transpose. 
   u,v are two vectors (reshape to 2-d arrays when applying A, or A')
*/
void fractal(double *p0, long nx, long ny, double dx, double r0){
    if((nx-1 & (nx-2)) != 0  || (ny-1 & (ny-2)) !=0 || nx != ny){
	error("nx=%ld, ny=%ld: they need to be 1+power of 2, and equal\n", nx, ny);
    }
    double power=5./3.;
    double coeff=6.88*pow(2*M_PI/0.5e-6, -2) * pow(r0, -power);
    double (*p)[nx]=(void*)p0;
    double sqrt2=sqrt(2.);
    double sqrt2i=1./sqrt2;
    long nx1=nx-1;
    long ny1=ny-1;
    double D=(nx1)*dx;
    double sigma2=0.5*F(sqrt2*D);//changed from sqrt(2) to sqrt(5);
    double c0=sigma2;
    {
	//First generate four outmost values.
	double fsqrt2d=F(sqrt2*D);
	double fd=F(D);
	double a=sqrt(4*sigma2-fd-0.5*fsqrt2d);
	double b=sqrt(fd-0.5*fsqrt2d);
	double c=sqrt(fsqrt2d);
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
#define QUA(p0,p1,p2,p3,p4) p0=qua1*(p1+p2+p3+p4)+qua0*p0
#define DIA(p0,p1,p2,p3,p4) p0=dia1*(p1+p2+p3+p4)+dia0*p0
#define TRI(p0,p1,p2,p3)    p0=(p1+p2)*tri1+p3*tri3+tri0*p0

    for(long step=nx1; step>1; step>>=1){
	double r=step*dx;
	double c1=sigma2-F(r*0.5)*0.5;
	double c2=sigma2-F(r*sqrt2i)*0.5;
	double c3=sigma2-F(r)*0.5;
	double c4=sigma2-F(r*sqrt2)*0.5;
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
    }
#undef QUA
#undef DIA
#undef TRI
}
/**
   Apply inverse fractal operator recursively. We simply compute the random
numbers from the turbulence in this step by just revert the computation in QUA
and TRI. Having reversed order as fractal().  */
void fractal_inv(double *p0, long nx, long ny, double dx, double r0){
    if((nx-1 & (nx-2)) != 0  || (ny-1 & (ny-2)) !=0 || nx != ny){
	error("nx=%ld, ny=%ld: they need to be 1+power of 2, and equal\n", nx, ny);
    }
    double power=5./3.;
    double coeff=6.88*pow(2*M_PI/0.5e-6, -2) * pow(r0, -power);
    double (*p)[nx]=(void*)p0;
    double sqrt2=sqrt(2.);
    double sqrt2i=1./sqrt2;
    long nx1=nx-1;
    long ny1=ny-1;
    double D=(nx1)*dx;
    double sigma2=0.5*F(sqrt2*D);//changed from sqrt(2) to sqrt(5);
    double c0=sigma2;
#define QUA(p0,p1,p2,p3,p4) p0=(p0-alpha1*(p1+p2+p3+p4))/alpha0;
#define TRI(p0,p1,p2,p3)    p0=(p0-((p1+p2)*tri1+p3*tri3))/tri0;
    long nreg=nx1>>1;//number of regions to do.
    for(long step=2; step<nx; step<<=1){
	double r=step*dx;
	double c1=sigma2-F(r*0.5)*0.5;
	double c2=sigma2-F(r*sqrt2i)*0.5;
	double c3=sigma2-F(r)*0.5;
	double c4=sigma2-F(r*sqrt2)*0.5;
	long step2=step>>1;
	for(long my=0; my<nreg; my++){
	    long offy=my*step;
	    for(long mx=0; mx<nreg; mx++){
		long offx=mx*step;
		//for case 3: Diamond configuration
		double alpha1=c1/(c0+2*c2+c3);
		double alpha0=sqrt(c0-4.*c1*alpha1);
		//for case 2: triangular case
		double trii=1/(c0*(c0+c3)-2*c2*c2);
		double tri1=c1*(c0-c2)*trii;
		double tri3=c1*(c0-2*c2+c3)*trii;
		double tri0=sqrt(c0-c1*c1*(3*c0-4*c2+c3)*trii);
		if(mx==0){//do the top one, triangular 
		    TRI(p[offy+step2][offx],
			p[offy][offx],
			p[offy+step][offx],
			p[offy+step2][offx+step2]);
		}
		
		//bottom one. 
		if(mx!=nreg-1){//diamond
		    QUA(p[offy+step2][offx+step],
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

		if(my==0){//do the left one, triangular 
		    TRI(p[offy][offx+step2],
			p[offy][offx],
			p[offy][offx+step],
			p[offy+step2][offx+step2]);
		}
		//right one.
		if(my!=nreg-1){//triangular
		    QUA(p[offy+step][offx+step2],
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
	for(long my=0; my<nreg; my++){
	    long offy=my*step;
	    for(long mx=0; mx<nreg; mx++){
		long offx=mx*step;
		//case 1: square configuration
		double alpha1=c2/(c0+2*c3+c4);
		double alpha0=sqrt(c0-4.*c2*alpha1);
		QUA(p[offy+step2][offx+step2],
		    p[offy][offx],
		    p[offy][offx+step],
		    p[offy+step][offx],
		    p[offy+step][offx+step]);
	    }
	}
	nreg>>=1;
    }
#undef QUA
#undef TRI
  {
	//Finally generate four outmost values.
	double fsqrt2d=F(sqrt2*D);
	double fd=F(D);
	double a=1./sqrt(4*sigma2-fd-0.5*fsqrt2d);
	double b=1./sqrt(fd-0.5*fsqrt2d);
	double c=2./sqrt(fsqrt2d);
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
}/**
   Apply transpose fractal operator recursively. The order are the same as in fractal_inv(). The only change being the definition of QUA and TRI
*/
void fractal_trans(double *p0, long nx, long ny, double dx, double r0){
    if((nx-1 & (nx-2)) != 0  || (ny-1 & (ny-2)) !=0 || nx != ny){
	error("nx=%ld, ny=%ld: they need to be 1+power of 2, and equal\n", nx, ny);
    }
    double power=5./3.;
    double coeff=6.88*pow(2*M_PI/0.5e-6, -2) * pow(r0, -power);
    double (*p)[nx]=(void*)p0;
    double sqrt2=sqrt(2.);
    double sqrt2i=1./sqrt2;
    long nx1=nx-1;
    long ny1=ny-1;
    double D=(nx1)*dx;
    double sigma2=0.5*F(sqrt2*D);//changed from sqrt(2) to sqrt(5);
    double c0=sigma2;
#define QUA(p0,p1,p2,p3,p4)			\
    double tmp=alpha1*p0;			\
    p1+=tmp;					\
    p2+=tmp;					\
    p3+=tmp;					\
    p4+=tmp;					\
    p0=alpha0*p0;

#define TRI(p0,p1,p2,p3)			\
    p1+=p0*tri1;				\
    p2+=p0*tri1;				\
    p3+=p0*tri3;				\
    p0=tri0*p0;

    long nreg=nx1>>1;//number of regions to do.
    for(long step=2; step<nx; step<<=1){
	double r=step*dx;
	double c1=sigma2-F(r*0.5)*0.5;
	double c2=sigma2-F(r*sqrt2i)*0.5;
	double c3=sigma2-F(r)*0.5;
	double c4=sigma2-F(r*sqrt2)*0.5;
	long step2=step>>1;
	for(long my=0; my<nreg; my++){
	    long offy=my*step;
	    for(long mx=0; mx<nreg; mx++){
		long offx=mx*step;
		//for case 3: Diamond configuration
		double alpha1=c1/(c0+2*c2+c3);
		double alpha0=sqrt(c0-4.*c1*alpha1);
		//for case 2: triangular case
		double trii=1/(c0*(c0+c3)-2*c2*c2);
		double tri1=c1*(c0-c2)*trii;
		double tri3=c1*(c0-2*c2+c3)*trii;
		double tri0=sqrt(c0-c1*c1*(3*c0-4*c2+c3)*trii);
		if(mx==0){//do the top one, triangular 
		    TRI(p[offy+step2][offx],
			p[offy][offx],
			p[offy+step][offx],
			p[offy+step2][offx+step2]);
		}
		
		//bottom one. 
		if(mx!=nreg-1){//diamond
		    QUA(p[offy+step2][offx+step],
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

		if(my==0){//do the left one, triangular 
		    TRI(p[offy][offx+step2],
			p[offy][offx],
			p[offy][offx+step],
			p[offy+step2][offx+step2]);
		}
		//right one.
		if(my!=nreg-1){//triangular
		    QUA(p[offy+step][offx+step2],
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
	for(long my=0; my<nreg; my++){
	    long offy=my*step;
	    for(long mx=0; mx<nreg; mx++){
		long offx=mx*step;
		//case 1: square configuration
		double alpha1=c2/(c0+2*c3+c4);
		double alpha0=sqrt(c0-4.*c2*alpha1);
		QUA(p[offy+step2][offx+step2],
		    p[offy][offx],
		    p[offy][offx+step],
		    p[offy+step][offx],
		    p[offy+step][offx+step]);
	    }
	}
	nreg>>=1;
    }
#undef QUA
#undef TRI
  {
	//Finally generate four outmost values.
	double fsqrt2d=F(sqrt2*D);
	double fd=F(D);
	double a=sqrt(4*sigma2-fd-0.5*fsqrt2d);
	double b=sqrt(fd-0.5*fsqrt2d);
	double c=sqrt(fsqrt2d);
	double *p1=&p[0][0];
	double *p2=&p[0][nx1];
	double *p3=&p[ny1][nx1];
	double *p4=&p[ny1][0];
	//takes the same form as fractal_inv(), except that a, b, c are defined like in fractal()
	double q1=0.5*( a**p1+a**p2+a**p3+a**p4);
	double q2=0.5*(-b**p1+b**p2-b**p3+b**p4);
	double q3=0.5*(-c**p1      +c**p3      );
	double q4=0.5*(      -c**p2      +c**p4);
	*p1=q1;
	*p2=q2;
	*p3=q3;
	*p4=q4;
    }
}
/**
   Apply transpose of inverse fractal operator recursively. Having the same order as fractal(). The QUA and TRI are changed to be like transpose of the operations in fractal_inv.
*/
void fractal_inv_trans(double *p0, long nx, long ny, double dx, double r0){
    if((nx-1 & (nx-2)) != 0  || (ny-1 & (ny-2)) !=0 || nx != ny){
	error("nx=%ld, ny=%ld: they need to be 1+power of 2, and equal\n", nx, ny);
    }
    double power=5./3.;
    double coeff=6.88*pow(2*M_PI/0.5e-6, -2) * pow(r0, -power);
    double (*p)[nx]=(void*)p0;
    double sqrt2=sqrt(2.);
    double sqrt2i=1./sqrt2;
    long nx1=nx-1;
    long ny1=ny-1;
    double D=(nx1)*dx;
    double sigma2=0.5*F(sqrt2*D);//changed from sqrt(2) to sqrt(5);
    double c0=sigma2;
    {
	//First generate four outmost values.
	double fsqrt2d=F(sqrt2*D);
	double fd=F(D);
	double a=1./sqrt(4*sigma2-fd-0.5*fsqrt2d);
	double b=1./sqrt(fd-0.5*fsqrt2d);
	double c=2./sqrt(fsqrt2d);
	double *p1=&p[0][0];
	double *p2=&p[0][nx1];
	double *p3=&p[ny1][nx1];
	double *p4=&p[ny1][0];
	//take the same form as in fractal(), except that a,b,c are defined like in fractal_inv().
	double q1=0.5*(a**p1-b**p2-c**p3);
	double q2=0.5*(a**p1+b**p2-c**p4);
	double q3=0.5*(a**p1-b**p2+c**p3);
	double q4=0.5*(a**p1+b**p2+c**p4);
	*p1=q1;
	*p2=q2;
	*p3=q3;
	*p4=q4;
    }
#define QUA(p0,p1,p2,p3,p4)			\
    double tmp=-alpha1/alpha0*p0;		\
    p1+=tmp; p2+=tmp; p3+=tmp; p4+=tmp;		\
    p0=p0/alpha0;
#define TRI(p0,p1,p2,p3)			\
    p1-=tri1/tri0*p0;				\
    p2-=tri1/tri0*p0;				\
    p3-=tri3/tri0*p0;				\
    p0=p0/tri0;
    long nreg=1;//number of regions to do
    for(long step=nx1; step>1; step>>=1){
	double r=step*dx;
	double c1=sigma2-F(r*0.5)*0.5;
	double c2=sigma2-F(r*sqrt2i)*0.5;
	double c3=sigma2-F(r)*0.5;
	double c4=sigma2-F(r*sqrt2)*0.5;
	long step2=step>>1;
	//finish all the black ones before we do white ones, which depend on this
	for(long my=0; my<nreg; my++){
	    long offy=my*step;
	    for(long mx=0; mx<nreg; mx++){
		long offx=mx*step;
		//case 1: square configuration
		double alpha1=c2/(c0+2*c3+c4);
		double alpha0=sqrt(c0-4.*c2*alpha1);
		QUA(p[offy+step2][offx+step2],
		    p[offy][offx],
		    p[offy][offx+step],
		    p[offy+step][offx],
		    p[offy+step][offx+step]);
	    }
	}
	//do all the white ones
	for(long my=0; my<nreg; my++){
	    long offy=my*step;
	    for(long mx=0; mx<nreg; mx++){
		long offx=mx*step;
		//for case 3: Diamond configuration
		double alpha1=c1/(c0+2*c2+c3);
		double alpha0=sqrt(c0-4.*c1*alpha1);
		//for case 2: triangular case
		double trii=1/(c0*(c0+c3)-2*c2*c2);
		double tri1=c1*(c0-c2)*trii;
		double tri3=c1*(c0-2*c2+c3)*trii;
		double tri0=sqrt(c0-c1*c1*(3*c0-4*c2+c3)*trii);
		if(mx==0){//do the top one, triangular 
		    TRI(p[offy+step2][offx],
			p[offy][offx],
			p[offy+step][offx],
			p[offy+step2][offx+step2]);
		}
		
		//bottom one. 
		if(mx!=nreg-1){//diamond
		    QUA(p[offy+step2][offx+step],
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

		if(my==0){//do the left one, triangular 
		    TRI(p[offy][offx+step2],
			p[offy][offx],
			p[offy][offx+step],
			p[offy+step2][offx+step2]);
		}
		//right one.
		if(my!=nreg-1){//triangular
		    QUA(p[offy+step][offx+step2],
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
	nreg<<=1;
    }
#undef QUA
#undef TRI
}

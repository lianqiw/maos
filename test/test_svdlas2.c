/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#include "../lib/aos.h"
#include "errno.h"
/*
  Implement and test SVD of sparse matrix based in las2.c in SVDPACK
*/
static void dspsvd(dmat **Sdiag, dmat **U, dmat **VT, const dsp *A);

int main(){
    if(!zfexist("FLM.bin")){
	info("FLM.bin does not exist\n");
	return 0;
    }
    dsp *A=dspread("FLM");
    dmat *s, *u, *vt;
    dspsvd(&s, &u, &vt, A);
}

/**************************************************************
 * Sparse svd via eigensystem of A'A matrix   		      *
 * The equivalent symmetric eigenvalue problem:               *
 *                                                            *
 *  B x = lambda x, where x' = (u',v'), lambda = sigma**2,    *
 *                                                            *
 *  B = A'A, and A is m (nrow) by n (ncol) (nrow >> ncol),    *
 *							      *
 *  so that {u, sqrt(lambda), v} is a singular triplet of A.  *
 *  (A' = transpose of A)				      *
 *                                                            *
 * global variables and common areas used by las2 and its     *
 * procedures.                                                *
 **************************************************************/

#define LMTNW   600000  /* max. size of working area allowed  */
#define NMAX    3000    /* bound on ncol, order of A          */
#define NZMAX   100000  /* bound on number of nonzeros in a   */

#define STORQ 1
#define RETRQ 2
#define STORP 3
#define RETRP 4
#define TRUE  1
#define FALSE 0
#define MAXLL 2

static 
long    ierr,           /* error flag                         */
	j,              /* number of lanczos steps taken      */
	neig,           /* number of ritz values stabilized   */
	nsig,           /* number of accepted ritz values     *
			 * based on kappa (relative accuracy) */
    	ncol,           /* number of columns of A             */
    	nrow,           /* number of rows of A                */
	mxvcount = 0;

/**************************************************************
 * pointers to areas holding input matrix which is stored in  *
 * harwell-boeing format.                                     *
 **************************************************************/
static 
long    *pointr = NULL, /* pointer to column start array      */
   	*rowind = NULL; /* pointer to row indices array       */
static 
real  *value = NULL;  /* pointer to nonzero values array    */

static 
real  rnm,            /* norm of the next residual vector   */
	anorm,
	tol,
	eps,            /* positive machine epsilon           */
	eps1,           /* roundoff estimate for dot product  *
			 * of two unit vector                 */
 	reps,
	eps34;

static 
real  *xv1 = NULL,    /* temp arrays needed for computing   */
	*xv2 = NULL,    /* singular vectors                   */
	*ztemp = NULL,

        *a = NULL;      /* pointer to area used by user-      *
			 * supplied procedure store and holds *
		  	 * lanczos vectors                    */

static
const char	*error[10] = {  /* error messages used by function    *
			 * check_parameters                   */
	    NULL,
	  " SORRY, YOUR MATRIX IS TOO BIG ",
	  " ***** ENDL MUST BE LESS THAN ENDR *****",
	  " ***** MAXPRS CANNOT EXCEED LANMAX *****",
	  " ***** N = NROW + NCOL MUST BE GREATER THAN ZERO *****",
	  " ***** LANMAX (NUMBER OF LANCZOS STEPS) IS INVALID *****",
	  " ***** MAXPRS (NUMBER OF IEGENPAIRS DESIRED) IS INVALID *****",
	  " ***** 6*N+4*LANMAX+1 + LANMAX*LANMAX CANNOT EXCEED NW *****",
	  " ***** 6*N+4*LANMAX+1 CANNOT EXCEED NW *****",
	    NULL};
static long check_parameters(long maxprs, long lanmax, long n, 
		      real endl, real endr, long vectors, 
			     long nnzero) ;
static long landr(long n, long lanmax, long maxprs, long nnzero, real endl, 
	   real endr, long vectors, real kappa, real *ritz, 
		  real *bnd, real *r);

static void ritvec(long n, real kappa, real *ritz, real *bnd, real *alf,
		   real *bet, real *w1, real *w2);

static void lanso(long n, long lanmax, long maxprs, real endl,
		  real endr, real *ritz, real *bnd, real *wptr[]);

static void lanczos_step(long n, long first, long last, real *wptr[],
		  real *alf, real *eta, real *oldeta,
			 real *bet, long *ll, long *enough);
static void ortbnd(real *alf, real *eta, real *oldeta, real *bet);

static void purge(long n, long ll, real *r, real *q, real *ra,  
		  real *qa, real *wrk, real *eta, real *oldeta);

static void stpone(long n, real *wrkptr[]);

static real startv(long n, real *wptr[]);
static real pythag(real a, real b);
static void error_bound(long *enough, real endl, real endr, 
		     real *ritz, real *bnd);

static void imtqlb(long n, real d[], real e[], real bnd[]);
static void imtql2(long nm, long n, real d[], real e[], real z[]);
static void machar(long *ibeta, long *it, long *irnd, long *machep, long *negep);
static void store(long n, long isw, long j, real *s);
static real fsign(real a,real b);
static real ddmax(real a, real b);
static real ddmin(real a, real b);
static long imin(long a, long b);
static long imax(long a,long b);
static void daxpy (long n,real da,real *dx,long incx,real *dy,long incy);
static real ddot(long n,real *dx,long incx,real *dy,long incy);
static void datx(long n,real da,real *dx,long incx,real *dy,long incy);
static void dsort2(long igap,long n,real *array1,real *array2);
static void dswap(long n,real *dx,long incx,real *dy,long incy);
static long idamax(long n,real *dx,long incx);
static void dscal(long n,real da,real *dx,long incx);
static void dcopy(long n,real *dx,long incx,real *dy,long incy);
static const dsp* A2=NULL;
/**************************************************************
 * multiplication of matrix B by vector x, where B = A'A,     *
 * and A is nrow by ncol (nrow >> ncol). Hence, B is of order *
 * n = ncol (y stores product vector).		              *
 **************************************************************/
static void opb(long n, real *x, real *y){
    (void)n;
    dmat *tmp=dnew(A2->nx, 1);
    dspmulvec(P(tmp), A2, x, 'n',1);
    memset(y, 0, sizeof(real)*A2->ny);
    dspmulvec(y, A2, P(tmp), 't',1);
    dfree(tmp);
    mxvcount +=2;
}
static void dspsvd(dmat **Sdiag, dmat **U, dmat **VT, const dsp *A){
    /* make a lanczos run; see landr for meaning of parameters */
    A2=A;/*Save it to global. */
    nrow=A->nx;
    ncol=A->ny;
    long i;
    long n=ncol;      /*Dimension of B=A'A; */
    long nn=ncol+nrow;
    long lanmax=ncol;    /*upper limit of desired number of Lanczos steps */
    long maxprs=ncol; /*upper limit of desired number of eigenpairs */
    long nnzero=A->nzmax;
    real endl=-1e-30;
    real endr=1e-30;
    long vectors=1;    /*1: indicates both eigenvalues and eigen vectors are wanted. */
    real kappa=1e-6; /*relative accuracy of ritz values acceptable as
			 eigenvalues of B (singular values of A)*/
   /*******************************************************************
     * allocate memory                                                 *
     * r      - work array                                        (n)  *
     * ritz   - array of ritz values                              (n)  *
     * bnd    - array of error bounds                             (n)  *
     * d      - array of approximate singular values of matrix A  (n)  *
     * ztemp  - work array for user function opb	       (nrow)  *
     * a      - storage area for Lanczos vectors     (n * (lanmax+2))  *
     *******************************************************************/
    long size1 = sizeof(real) * (6 * n + nrow + nnzero + (n * lanmax));

    if (!(value  = (real *) malloc(size1))){
	perror("MALLOC FAILED in MAIN()");
	exit(errno);
    }
    real* tptr1 = value;
    real *r,*ritz,*bnd,*d;
    rowind = pointr + ncol + 1;
    tptr1 += nnzero;
    r      = tptr1;
    tptr1 += n;
    ritz   = tptr1;
    tptr1 += n;
    bnd    = tptr1;
    tptr1 += n;
    d      = tptr1;
    tptr1 += n;
    ztemp  = tptr1;
    tptr1 += nrow;
    a      = tptr1;
    *Sdiag=dnew(maxprs,1);
    dmat *V=dnew(nrow, maxprs); dmat* pV=V;
    TIC;tic;
    if(landr(n, lanmax, maxprs, nnzero, endl, endr, vectors, kappa, ritz, bnd, r)){ 
	warning("landr failed\n");
	/*clean up */
    }
    toc("landr");
    {
	/* print ritz values and error bounds */
	dbg("\n");
	dbg(" ...... NUMBER OF LANCZOS STEPS = %3ld       NEIG = %3ld\n", j+1, neig);
	dbg(" ...... \n");
	dbg(" ......         COMPUTED RITZ VALUES  (ERROR BNDS)\n");
	dbg(" ...... \n");
	for (i = 0; i <= j; i++)
	    dbg("...... %3ld   %22.14E  (%11.2E)\n", i + 1, ritz[i], bnd[i]);
    }

    if(vectors){/*we do want eigen vectors */
	size1 = sizeof(real) * nrow;
	
	long id = 0;

	for (i = 0; i < nsig; i++) {

	    /* multiply by matrix B first */
	    opb(n, &xv1[id], xv2);
	    real tmp0 = ddot(n, &xv1[id], 1, xv2, 1);
	    daxpy(n, -tmp0, &xv1[id], 1, xv2, 1);
	    tmp0 = sqrt(tmp0);
	    real xnorm = sqrt(ddot(n, xv2, 1, xv2, 1));
	    long ida = id + ncol;

	    /* multiply by matrix A to get (scaled) left s-vector */
	    /*opa(&xv1[id], &xv1[ida]); */
	    mxvcount+=1;
	    memset(&xv1[ida], 0, sizeof(real)*nrow);
	    dspmulvec(&xv1[ida], A, &xv1[id],'n', 1);
	    real tmp1 = 1.0 / tmp0;
	    dscal(nrow, tmp1, &xv1[ida], 1);
	    xnorm *= tmp1;
	    bnd[i] = xnorm;
	    d[i] = tmp0;
	    P(*Sdiag,i)=tmp0;
	    /* write left s-vector to output file */
	    /*write(fp_out2, (char *)&xv1[ida], size1); */
	    memcpy(&pV[i], &xv1[ida], size1);
	    id += nn;
	}
	toc("done");
	long count1=(mxvcount-nsig)/2 + nsig;
	long count2=(mxvcount-nsig)/2;
	dbg(" ...... \n");
	dbg(" ...... NO. MULTIPLICATIONS BY A  =%10ld\n", count1);
	dbg(" ...... NO. MULT. BY TRANSPOSE(A) =%10ld\n", count2);
	dbg("\n");
	dbg(" ...... \n");
	dbg(" ......        NSIG = %4ld\n", nsig);
	dbg(" ...... \n");
	dbg(" ......         COMPUTED S-VALUES     (RES. NORMS)\n");
	dbg(" ...... \n");
	for (i = 0; i < nsig; i++)
	    dbg(" ...... %3ld   %22.14E  (%11.2E)\n",
		    i + 1, d[i], bnd[i]);
    }
    else {
	for (i = j; i >= 0; i--)
	    if (bnd[i] > kappa * fabs(ritz[i])) break;
	nsig = j - i;

	long count1=(mxvcount-nsig)/2 + nsig;
	long count2=(mxvcount-nsig)/2;
	dbg(" ...... \n");
	dbg(" ...... NO. MULTIPLICATIONS BY A  =%10ld\n", count1);
	dbg(" ...... NO. MULT. BY TRANSPOSE(A) =%10ld\n", count2);
	dbg("\n");
	dbg(" ...... \n");
	dbg(" ......         NSIG = %4ld\n" , nsig);
	dbg(" ...... \n");
	dbg(" ......         COMPUTED S-VALUES   (ERROR BNDS)\n");
	dbg(" ...... \n");

	long k = j + 1 - nsig;
	for (i = 1 ; i <= nsig; i++) {
	    dbg(" ...... %3ld   %22.14E  (%11.2E)\n", 
		    i, sqrt(ritz[k]), bnd[k]);
	    k++;
	    P(*Sdiag,i)=sqrt(ritz[k]);
	}
    }
    if (vectors) {
	free(xv1);
	free(xv2);
	*VT=dtrans(V);
	*U=NULL;
	dspmm(U, A, V, "nn",1);
	dmat *SdiagI=ddup(*Sdiag); 
	dcwpow(SdiagI, -1);
	dmuldiag(*U, SdiagI);
	dfree(V);
	dfree(SdiagI);
    }
    writebin(*U,"U");
    writebin(*VT,"VT");
    writebin(*Sdiag,"Sdiag");
    free(value);
}
/*
  The following code are taken from las2.c in SVDPACK
*/
   
/*************************************************************************
                           (c) Copyright 1993
                        University of Tennessee
                          All Rights Reserved                          
 *************************************************************************/

/***********************************************************************
 *								       *
 *		      check_parameters()			       *
 *								       *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------
   Function validates input parameters and returns error code (long)  

   Parameters 
   ----------
  (input)
   maxprs   upper limit of desired number of eigenpairs of B           
   lanmax   upper limit of desired number of lanczos steps             
   n        dimension of the eigenproblem for matrix B               
   endl     left end of interval containing unwanted eigenvalues of B
   endr     right end of interval containing unwanted eigenvalues of B
   vectors  1 indicates both eigenvalues and eigenvectors are wanted 
            and they can be found in lav2; 0 indicates eigenvalues only
   nnzero   number of nonzero elements in input matrix (matrix A)      
                                                                      
***********************************************************************/

static long check_parameters(long maxprs, long lanmax, long n, 
		      real endl, real endr, long vectors, 
		      long nnzero) 
{
    long error_inbex, ncells;
    error_inbex = 0;

    /* assuming that nrow >= ncol... */
    if (ncol >= NMAX || nnzero > NZMAX) error_inbex = 1;
    else if (endl >= endr)  error_inbex = 2;
    else if (maxprs > lanmax)  error_inbex = 3;
    else if (n <= 0)  error_inbex = 4;
    else if (lanmax <= 0 || lanmax > n)  error_inbex = 5;
    else if (maxprs <= 0 || maxprs > lanmax)  error_inbex = 6;
    else {
	if (vectors) {
	    ncells = 6 * n + 4 * lanmax + 1 + lanmax * lanmax;
	    if (ncells > LMTNW) error_inbex = 7;
	}
	else {
	    ncells = 6 * n + 4 * lanmax + 1;
	    if (ncells > LMTNW) error_inbex = 8;
	}
    }
    if (error_inbex) dbg("%s\n", error[error_inbex]);
    return(error_inbex);
}

/***********************************************************************
 *                                                                     *
 *				landr()				       *
 *        Lanczos algorithm with selective orthogonalization           *
 *                    Using Simon's Recurrence                         *
 *                       (real precision)                            *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   landr() is the LAS2 driver routine that, upon entry,
     (1)  checks for the validity of input parameters of the 
	  B-eigenproblem 
     (2)  determines several machine constants
     (3)  makes a Lanczos run
     (4)  calculates B-eigenvectors (singular vectors of A) if requested 
	  by user


   arguments
   ---------

   (input)
   n        dimension of the eigenproblem for A'A
   lanmax   upper limit of desired number of Lanczos steps
   maxprs   upper limit of desired number of eigenpairs
   nnzero   number of nonzeros in matrix A
   endl     left end of interval containing unwanted eigenvalues of B
   endr     right end of interval containing unwanted eigenvalues of B
   vectors  1 indicates both eigenvalues and eigenvectors are wanted
              and they can be found in output file lav2; 
	    0 indicates only eigenvalues are wanted
   kappa    relative accuracy of ritz values acceptable as eigenvalues
	      of B (singular values of A)
   r        work array

   (output)
   j        number of Lanczos steps actually taken                     
   neig     number of ritz values stabilized                           
   ritz     array to hold the ritz values                              
   bnd      array to hold the error bounds


   External parameters
   -------------------

   Defined and documented in las2.h


   local parameters
   -------------------

   ibeta    radix for the floating-point representation
   it       number of base ibeta digits in the floating-point significand
   irnd     floating-point addition rounded or chopped
   machep   machine relative precision or round-off error
   negeps   largest negative integer
   wptr	    array of pointers each pointing to a work space


   Functions used
   --------------

   MISC         ddmax, machar, check_parameters
   LAS2         ritvec, lanso

***********************************************************************/

static long landr(long n, long lanmax, long maxprs, long nnzero, real endl, 
	   real endr, long vectors, real kappa, real *ritz, 
	   real *bnd, real *r)

{
    long i, size, ibeta, it, irnd, machep, negep;
    real *wptr[10], *tptr, *tptr2;

    /* data validation */
    if (check_parameters(maxprs, lanmax, n, endl, endr, vectors, nnzero))
	return(-1);

    /* Compute machine precision */ 
    machar(&ibeta, &it, &irnd, &machep, &negep);

    eps1 = eps * sqrt( (real) n );
    reps = sqrt(eps);
    eps34 = reps * sqrt(reps);

    /* allocate work area and initialize pointers         *
     * ptr             symbolic name         size         *
     * wptr[0]             r                  n           *
     * wptr[1]             q		     n           *
     * wptr[2]             q_previous         n           *
     * wptr[3]             p		     n           *
     * wptr[4]             p_previous         n           *
     * wptr[5]             wrk                n           *
     * wptr[6]             alf              lanmax        *
     * wptr[7]             eta              lanmax        *
     * wptr[8]             oldeta           lanmax        *
     * wptr[9]             bet              lanmax+1      */

    size = 5 * n + (lanmax * 4 + 1);
    tptr = NULL;
    if (!(tptr = (real*)malloc(size * sizeof(real)))){
	perror("FIRST MALLOC FAILED in LANDR()");
	raise(errno);
    }
    tptr2 = tptr;
    wptr[0] = r;
    for (i = 1; i <= 5; i++) {
	wptr[i] = tptr;
	tptr += n;
    }
    for (i = 6; i <= 9; i++) {
	wptr[i] = tptr;
	tptr += lanmax;
    }
    lanso(n, lanmax, maxprs, endl, endr, ritz, bnd, wptr);

    /* compute eigenvectors */
    if (vectors) {
	dbg("j=%ld\n",j);
	xv1 = mymalloc((nrow+ncol)*(j+1) ,real);
	xv2 = mymalloc(ncol ,real);

	kappa = ddmax(fabs(kappa), eps34);
	ritvec(n, kappa, ritz, bnd, wptr[6], wptr[9], wptr[4], wptr[5]);
    }
    free(tptr2);
    warning("Change above\n");
    return(0);
}

/***********************************************************************
 *                                                                     *
 *                        ritvec()                                     *
 * 	    Function computes the singular vectors of matrix A	       *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   This function is invoked by landr() only if eigenvectors of the A'A
   eigenproblem are desired.  When called, ritvec() computes the 
   singular vectors of A and writes the result to an unformatted file.


   Parameters
   ----------

   (input)
   nrow       number of rows of A
   j	      number of Lanczos iterations performed
   fp_out2    pointer to unformatted output file
   n	      dimension of matrix A
   kappa      relative accuracy of ritz values acceptable as 
		eigenvalues of A'A
   ritz       array of ritz values
   bnd        array of error bounds
   alf        array of diagonal elements of the tridiagonal matrix T
   bet        array of off-diagonal elements of T
   w1, w2     work space

   (output)
   xv1        array of eigenvectors of A'A (right singular vectors of A)
   ierr	      error code
              0 for normal return from imtql2()
	      k if convergence did not occur for k-th eigenvalue in
	        imtql2()
   nsig       number of accepted ritz values based on kappa

   (local)
   s	      work array which is initialized to the identity matrix
	      of order (j + 1) upon calling imtql2().  After the call,
	      s contains the orthonormal eigenvectors of the symmetric 
	      tridiagonal matrix T

   Functions used
   --------------

   BLAS		dscal, dcopy, daxpy
   USER		store
   		imtql2

***********************************************************************/

static void ritvec(long n, real kappa, real *ritz, real *bnd, real *alf,
	    real *bet, real *w1, real *w2)

{
    long js, jsq, i, k, id, id2, tmp;
    real *s;

    js = j + 1;
    jsq = js * js;
    //size = sizeof(real) * n;

    if(!(s = (real *) malloc (jsq * sizeof(real)))) {
	perror("MALLOC FAILED in RITVEC()");
	raise(errno);
    }

    /* initialize s to an identity matrix */
    for (i = 0; i < jsq; i++) s[i] = 0.0;
    for (i = 0; i < jsq; i+= (js+1)) s[i] = 1.0;
    dcopy(js, alf, 1, w1, -1);
    dcopy(j, &bet[1], 1, &w2[1], -1);

    /* on return from imtql2(), w1 contains eigenvalues in ascending 
     * order and s contains the corresponding eigenvectors */
    imtql2(js, js, w1, w2, s);
    if (ierr) return;

    /*write(fp_out2, (char *)&n, sizeof(n));
    write(fp_out2, (char *)&js, sizeof(js));
    write(fp_out2, (char *)&kappa, sizeof(kappa));
    */
    id = 0;
    nsig = 0;
    id2 = jsq - js;
    for (k = 0; k < js; k++) {
	tmp = id2;
	if (bnd[k] <= kappa * fabs(ritz[k]) && k > js-neig-1) {
	    for (i = 0; i < n; i++) w1[i] = 0.0;
	    for (i = 0; i < js; i++) {
		store(n, RETRQ, i, w2);
		daxpy(n, s[tmp], w2, 1, w1, 1);
		tmp -= js;
	    }
	    /*write(fp_out2, (char *)w1, size); */

	    /* store the w1 vector row-wise in array xv1;   
	     * size of xv1 is (j+1) * (nrow+ncol) elements 
	     * and each vector, even though only ncol long,
	     * will have (nrow+ncol) elements in xv1.      
	     * It is as if xv1 is a 2-d array (j+1) by     
	     * (nrow+ncol) and each vector occupies a row  */

	    for (i = 0; i < n; i++) {
		xv1[id++] = w1[i];
	    }
	    id += nrow;
	    nsig += 1;
	}
	id2++;
    }
    free(s);
    return;
}


/***********************************************************************
 *                                                                     *
 *                          lanso()                                    *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   Function determines when the restart of the Lanczos algorithm should 
   occur and when it should terminate.

   Arguments 
   ---------

   (input)
   n         dimension of the eigenproblem for matrix B
   lanmax    upper limit of desired number of lanczos steps           
   maxprs    upper limit of desired number of eigenpairs             
   endl      left end of interval containing unwanted eigenvalues
   endr      right end of interval containing unwanted eigenvalues
   ritz      array to hold the ritz values                       
   bnd       array to hold the error bounds                          
   wptr      array of pointers that point to work space:            
  	       wptr[0]-wptr[5]  six vectors of length n		
  	       wptr[6] array to hold diagonal of the tridiagonal matrix T
  	       wptr[9] array to hold off-diagonal of T	
  	       wptr[7] orthogonality estimate of Lanczos vectors at 
		 step j
 	       wptr[8] orthogonality estimate of Lanczos vectors at 
		 step j-1

   (output)
   j         number of Lanczos steps actually taken
   neig      number of ritz values stabilized
   ritz      array to hold the ritz values
   bnd       array to hold the error bounds
   ierr      (globally declared) error flag
	     ierr = 8192 if stpone() fails to find a starting vector
	     ierr = k if convergence did not occur for k-th eigenvalue
		    in imtqlb()
	     ierr = 0 otherwise


   Functions used
   --------------

   LAS		stpone, error_bound, lanczos_step
   MISC		dsort2
   UTILITY	imin, imax

***********************************************************************/

static void lanso(long n, long lanmax, long maxprs, real endl,
	   real endr, real *ritz, real *bnd, real *wptr[])

{
    real *alf, *eta, *oldeta, *bet, *wrk;
    long ll, first, last, ENOUGH, id1, id2, id3, i, l;

    //real *r = wptr[0];
    alf = wptr[6];
    eta = wptr[7];
    oldeta = wptr[8];
    bet = wptr[9];
    wrk = wptr[5];
    j = 0;

    /* take the first step */
    stpone(n, wptr);
    if (!rnm || ierr) return;
    eta[0] = eps1;
    oldeta[0] = eps1;
    ll = 0;
    first = 1;
    last = imin(maxprs + imax(8,maxprs), lanmax);
    ENOUGH = FALSE;
    id1 = 0;
    while (id1 < maxprs && !ENOUGH) {
	if (rnm <= tol) rnm = 0.0;

	/* the actual lanczos loop */
	lanczos_step(n, first, last, wptr, alf, eta, oldeta, bet, &ll,
		     &ENOUGH);
	if (ENOUGH) j = j - 1;
	else j = last - 1;
	first = j + 1;
	bet[j+1] = rnm;

	/* analyze T */
	l = 0;
	for (id2 = 0; id2 < j; id2++) {
	    if (l > j) break;
	    for (i = l; i <= j; i++) if (!bet[i+1]) break;
	    if (i > j) i = j;

	    /* now i is at the end of an unreduced submatrix */
	    dcopy(i-l+1, &alf[l],   1, &ritz[l],  -1);
	    dcopy(i-l,   &bet[l+1], 1, &wrk[l+1], -1);

	    imtqlb(i-l+1, &ritz[l], &wrk[l], &bnd[l]);

	    if (ierr) {
		printf("IMTQLB FAILED TO CONVERGE (IERR = %ld)\n", ierr);
		printf("L = %ld    I = %ld\n", l, i);
		for (id3 = l; id3 <= i; id3++) 
		    printf("%ld  %lg  %lg  %lg\n",
			   id3, ritz[id3], wrk[id3], bnd[id3]);
	    }
	    for (id3 = l; id3 <= i; id3++) 
		bnd[id3] = rnm * fabs(bnd[id3]);
	    l = i + 1;
	}

	/* sort eigenvalues into increasing order */
	dsort2((j+1) / 2, j + 1, ritz, bnd);

	/* massage error bounds for very close ritz values */
	error_bound(&ENOUGH, endl, endr, ritz, bnd);

	/* should we stop? */
	if (neig < maxprs) {
	    if (!neig) last = first + 9;
	    else last = first + imax(3, 1 + ((j-5) * (maxprs-neig)) / neig);
	    last = imin(last, lanmax);
	}
	else ENOUGH = TRUE;
	ENOUGH = ENOUGH || first >= lanmax;
	id1++;
    }
    store(n, STORQ, j, wptr[1]);
    return;
}

/***********************************************************************
 *                                                                     *
 *			lanczos_step()                                 *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   Function embodies a single Lanczos step

   Arguments 
   ---------

   (input)
   n        dimension of the eigenproblem for matrix B
   first    start of inbex through loop				      
   last     end of inbex through loop				     
   wptr	    array of pointers pointing to work space		    
   alf	    array to hold diagonal of the tridiagonal matrix T
   eta      orthogonality estimate of Lanczos vectors at step j   
   oldeta   orthogonality estimate of Lanczos vectors at step j-1
   bet      array to hold off-diagonal of T                     
   ll       number of intitial Lanczos vectors in local orthog. 
              (has value of 0, 1 or 2)			
   enough   stop flag			

   Functions used
   --------------

   BLAS		ddot, dscal, daxpy, datx, dcopy
   USER		store
   LAS		purge, ortbnd, startv
   UTILITY	imin, imax

***********************************************************************/

static void lanczos_step(long n, long first, long last, real *wptr[],
		  real *alf, real *eta, real *oldeta,
		  real *bet, long *ll, long *enough)

{
    real t, *mid;
    long i;

    for (j=first; j<last; j++) {
	mid     = wptr[2];
	wptr[2] = wptr[1];
	wptr[1] = mid;
	mid     = wptr[3];
	wptr[3] = wptr[4];
	wptr[4] = mid;

	store(n, STORQ, j-1, wptr[2]);
	if (j-1 < MAXLL) store(n, STORP, j-1, wptr[4]);
	bet[j] = rnm;

	/* restart if invariant subspace is found */
	if (!bet[j]) {
	    rnm = startv(n, wptr);
	    if (ierr) return;
	    if (!rnm) *enough = TRUE;
	}
	if (*enough) {

	    /* bug fix supplied by Doug Rohde (MIT); 
	       email contact: dr+svd@tedlab.mit.edu  (Feb 2004) */

	    mid     = wptr[2];
	    wptr[2] = wptr[1];
	    wptr[1] = mid;
	    break;
	}

	/* take a lanczos step */
	t = 1.0 / rnm;
	datx(n, t, wptr[0], 1, wptr[1], 1);
	dscal(n, t, wptr[3], 1);
	opb(n, wptr[3], wptr[0]);
	daxpy(n, -rnm, wptr[2], 1, wptr[0], 1);
	alf[j] = ddot(n, wptr[0], 1, wptr[3], 1);
	daxpy(n, -alf[j], wptr[1], 1, wptr[0], 1);

	/* orthogonalize against initial lanczos vectors */
	if (j <= MAXLL && (fabs(alf[j-1]) > 4.0 * fabs(alf[j])))
	    *ll = j;  
	for (i=0; i < imin(*ll, j-1); i++) {
	    store(n, RETRP, i, wptr[5]);
	    t = ddot(n, wptr[5], 1, wptr[0], 1);
	    store(n, RETRQ, i, wptr[5]);
	    daxpy(n, -t, wptr[5], 1, wptr[0], 1);
	    eta[i] = eps1;
	    oldeta[i] = eps1;
	}

	/* extended local reorthogonalization */
	t = ddot(n, wptr[0], 1, wptr[4], 1);
	daxpy(n, -t, wptr[2], 1, wptr[0], 1);
	if (bet[j] > 0.0) bet[j] = bet[j] + t;
	t = ddot(n, wptr[0], 1, wptr[3], 1);
	daxpy(n, -t, wptr[1], 1, wptr[0], 1);
	alf[j] = alf[j] + t;
	dcopy(n, wptr[0], 1, wptr[4], 1);
	rnm = sqrt(ddot(n, wptr[0], 1, wptr[4], 1));
	anorm = bet[j] + fabs(alf[j]) + rnm;
	tol = reps * anorm;

	/* update the orthogonality bounds */
	ortbnd(alf, eta, oldeta, bet);

	/* restore the orthogonality state when needed */
	purge(n,*ll,wptr[0],wptr[1],wptr[4],wptr[3],wptr[5],eta,oldeta);
	if (rnm <= tol) rnm = 0.0;
    }
    return;
}

/***********************************************************************
 *                                                                     *
 *                          ortbnd()                                   *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   Funtion updates the eta recurrence

   Arguments 
   ---------

   (input)
   alf      array to hold diagonal of the tridiagonal matrix T         
   eta      orthogonality estimate of Lanczos vectors at step j        
   oldeta   orthogonality estimate of Lanczos vectors at step j-1     
   bet      array to hold off-diagonal of T                          
   n        dimension of the eigenproblem for matrix B		    
   j        dimension of T					  
   rnm	    norm of the next residual vector			 
   eps1	    roundoff estimate for dot product of two unit vectors

   (output)
   eta      orthogonality estimate of Lanczos vectors at step j+1     
   oldeta   orthogonality estimate of Lanczos vectors at step j        


   Functions used
   --------------

   BLAS		dswap

***********************************************************************/

static void ortbnd(real *alf, real *eta, real *oldeta, real *bet)

{
    long i;
    if (j < 1) return;
    if (rnm) {
	if (j > 1) {
	    oldeta[0] = (bet[1] * eta[1] + (alf[0]-alf[j]) * eta[0] -
			 bet[j] * oldeta[0]) / rnm + eps1;
	}
	for (i=1; i<=j-2; i++) 
	    oldeta[i] = (bet[i+1] * eta[i+1] + (alf[i]-alf[j]) * eta[i] +
			 bet[i] * eta[i-1] - bet[j] * oldeta[i])/rnm + eps1;
    }
    oldeta[j-1] = eps1;
    dswap(j, oldeta, 1, eta, 1);  
    eta[j] = eps1;
    return;
}


/***********************************************************************
 *                                                                     *
 *				purge()                                *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   Function examines the state of orthogonality between the new Lanczos
   vector and the previous ones to decide whether re-orthogonalization 
   should be performed


   Arguments 
   ---------

   (input)
   n        dimension of the eigenproblem for matrix B		       
   ll       number of intitial Lanczos vectors in local orthog.       
   r        residual vector to become next Lanczos vector            
   q        current Lanczos vector			           
   ra       previous Lanczos vector
   qa       previous Lanczos vector
   wrk      temporary vector to hold the previous Lanczos vector
   eta      state of orthogonality between r and prev. Lanczos vectors 
   oldeta   state of orthogonality between q and prev. Lanczos vectors
   j        current Lanczos step				     

   (output)
   r	    residual vector orthogonalized against previous Lanczos 
	      vectors
   q        current Lanczos vector orthogonalized against previous ones


   Functions used
   --------------

   BLAS		daxpy,  dcopy,  idamax,  ddot
   USER		store

***********************************************************************/

static void purge(long n, long ll, real *r, real *q, real *ra,  
	   real *qa, real *wrk, real *eta, real *oldeta)

{
    real t, tq, tr, reps1;
    long k, iteration, flag, i;

    if (j < ll+2) return; 

    k = idamax(j - (ll+1), &eta[ll], 1) + ll;
    if (fabs(eta[k]) > reps) {
	reps1 = eps1 / reps;
	iteration = 0;
	flag = TRUE;
	while (iteration < 2 && flag) {
	    if (rnm > tol) {

		/* bring in a lanczos vector t and orthogonalize both 
		 * r and q against it */
		tq = 0.0;
		tr = 0.0;
		for (i = ll; i < j; i++) {
		    store(n,  RETRQ,  i,  wrk);
		    t   = -ddot(n, qa, 1, wrk, 1);
		    tq += fabs(t);
		    daxpy(n,  t,  wrk,  1,  q,  1);
		    t   = -ddot(n, ra, 1, wrk, 1);
		    tr += fabs(t);
		    daxpy(n, t, wrk, 1, r, 1);
		}
		dcopy(n, q, 1, qa, 1);
		t   = -ddot(n, r, 1, qa, 1);
		tr += fabs(t);
		daxpy(n, t, q, 1, r, 1);
		dcopy(n, r, 1, ra, 1);
		rnm = sqrt(ddot(n, ra, 1, r, 1));
		if (tq <= reps1 && tr <= reps1 * rnm) flag = FALSE;
	    }
	    iteration++;
	}
	for (i = ll; i <= j; i++) { 
	    eta[i] = eps1;
	    oldeta[i] = eps1;
	}
    }
    return;
}

/***********************************************************************
 *                                                                     *
 *                         stpone()                                    *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   Function performs the first step of the Lanczos algorithm.  It also
   does a step of extended local re-orthogonalization.

   Arguments 
   ---------

   (input)
   n      dimension of the eigenproblem for matrix B

   (output)
   ierr   error flag
   wptr   array of pointers that point to work space that contains
	    wptr[0]             r[j]
	    wptr[1]             q[j]
	    wptr[2]             q[j-1]
	    wptr[3]             p
	    wptr[4]             p[j-1]
	    wptr[6]             diagonal elements of matrix T 


   Functions used
   --------------

   BLAS		daxpy, datx, dcopy, ddot, dscal
   USER		store, opb
   LAS		startv

***********************************************************************/

static void stpone(long n, real *wrkptr[])

{
    real t, *alf;
    alf = wrkptr[6];

    /* get initial vector; default is random */
    rnm = startv(n, wrkptr);
    if (rnm == 0.0 || ierr != 0) return;

    /* normalize starting vector */
    t = 1.0 / rnm;
    datx(n, t, wrkptr[0], 1, wrkptr[1], 1);
    dscal(n, t, wrkptr[3], 1);

    /* take the first step */
    opb(n, wrkptr[3], wrkptr[0]);
    alf[0] = ddot(n, wrkptr[0], 1, wrkptr[3], 1);
    daxpy(n, -alf[0], wrkptr[1], 1, wrkptr[0], 1);
    t = ddot(n, wrkptr[0], 1, wrkptr[3], 1);
    daxpy(n, -t, wrkptr[1], 1, wrkptr[0], 1);
    alf[0] += t;
    dcopy(n, wrkptr[0], 1, wrkptr[4], 1);
    rnm = sqrt(ddot(n, wrkptr[0], 1, wrkptr[4], 1));
    anorm = rnm + fabs(alf[0]);
    tol = reps * anorm;
    return;
}

/***********************************************************************
 *                                                                     *
 *                         startv()                                    *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   Function delivers a starting vector in r and returns |r|; it returns 
   zero if the range is spanned, and ierr is non-zero if no starting 
   vector within range of operator can be found.

   Parameters 
   ---------

   (input)
   n      dimension of the eigenproblem matrix B
   wptr   array of pointers that point to work space
   j      starting inbex for a Lanczos run
   eps    machine epsilon (relative precision)

   (output)
   wptr   array of pointers that point to work space that contains
	  r[j], q[j], q[j-1], p[j], p[j-1]
   ierr   error flag (nonzero if no starting vector can be found)

   Functions used
   --------------

   BLAS		ddot, dcopy, daxpy
   USER		opb, store
   MISC		random

***********************************************************************/

static real startv(long n, real *wptr[])

{
    real rnm2, *r, t;
    /*long irand; */
    rand_t rstat;
    long id, i;

    /* get initial vector; default is random */
    rnm2 = ddot(n, wptr[0], 1, wptr[0], 1);
    /*irand = 918273 + j; */
    seed_rand(&rstat, 918273);
    r = wptr[0];
    for (id = 0; id < 3; id++) {
	if (id > 0 
	    || j > 0 
	    || rnm2 == 0.0) {
	    for (i = 0; i < n; i++) 
		r[i] = lrand(&rstat);
	}
	dcopy(n, wptr[0], 1, wptr[3], 1);

	/* apply operator to put r in range (essential if m singular) */
	opb(n, wptr[3], wptr[0]);
	dcopy(n, wptr[0], 1, wptr[3], 1);
	rnm2 = ddot(n, wptr[0], 1, wptr[3], 1);
	if (rnm2 > 0.0) break;
    }

    /* fatal error */
    if (rnm2 <= 0.0) {
	ierr = 8192;
	return(-1);
    }
    if (j > 0) {
	for (i = 0; i < j; i++) {
	    store(n, RETRQ, i, wptr[5]);
	    t = -ddot(n, wptr[3], 1, wptr[5], 1);
	    daxpy(n, t, wptr[5], 1, wptr[0], 1);
	}

	/* make sure q[j] is orthogonal to q[j-1] */
	t = ddot(n, wptr[4], 1, wptr[0], 1);
	daxpy(n, -t, wptr[2], 1, wptr[0], 1);
	dcopy(n, wptr[0], 1, wptr[3], 1);
	t = ddot(n, wptr[3], 1, wptr[0], 1);
	if (t <= eps * rnm2) t = 0.0;
	rnm2 = t;
    }
    return(sqrt(rnm2));
}


/************************************************************** 
 *							      *
 * Function finds sqrt(a^2 + b^2) without overflow or         *
 * destructive underflow.				      *
 *							      *
 **************************************************************/ 
/************************************************************** 

   Funtions used
   -------------

   UTILITY	ddmax, ddmin

**************************************************************/ 

static real pythag(real a1, real b)

{
    real p, r, s, t, u, temp;

    p = ddmax(fabs(a1), fabs(b));
    if (p != 0.0) {
	temp = ddmin(fabs(a1), fabs(b)) / p;
	r = temp * temp; 
	t = 4.0 + r;
	while (t != 4.0) {
	    s = r / t;
	    u = 1.0 + 2.0 * s;
	    p *= u;
	    temp = s / u;
	    r *= temp * temp;
	    t = 4.0 + r;
	}
    }
    return(p);
}

/***********************************************************************
 *                                                                     *
 *			error_bound()                                  *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   Function massages error bounds for very close ritz values by placing 
   a gap between them.  The error bounds are then refined to reflect 
   this.


   Arguments 
   ---------

   (input)
   endl     left end of interval containing unwanted eigenvalues
   endr     right end of interval containing unwanted eigenvalues
   ritz     array to store the ritz values
   bnd      array to store the error bounds
   enough   stop flag


   Functions used
   --------------

   BLAS		idamax
   UTILITY	dmin

***********************************************************************/

static void error_bound(long *enough, real endl, real endr, 
		 real *ritz, real *bnd)
{
    long mid, i;
    real gapl, gap;

    /* massage error bounds for very close ritz values */
    mid = idamax(j + 1, bnd, 1);

    for (i=((j+1) + (j-1)) / 2; i >= mid + 1; i -= 1)
	if (fabs(ritz[i-1] - ritz[i]) < eps34 * fabs(ritz[i])) 
	    if (bnd[i] > tol && bnd[i-1] > tol) {
		bnd[i-1] = sqrt(bnd[i] * bnd[i] + bnd[i-1] * bnd[i-1]);
		bnd[i] = 0.0;
	    }
	 

    for (i=((j+1) - (j-1)) / 2; i <= mid - 1; i +=1 ) 
	if (fabs(ritz[i+1] - ritz[i]) < eps34 * fabs(ritz[i])) 
	    if (bnd[i] > tol && bnd[i+1] > tol) {
		bnd[i+1] = sqrt(bnd[i] * bnd[i] + bnd[i+1] * bnd[i+1]);
		bnd[i] = 0.0;
	    }

    /* refine the error bounds */
    neig = 0;
    gapl = ritz[j] - ritz[0];
    for (i = 0; i <= j; i++) {
	gap = gapl;
	if (i < j) gapl = ritz[i+1] - ritz[i];
	gap = ddmin(gap, gapl);
	if (gap > bnd[i]) bnd[i] = bnd[i] * (bnd[i] / gap);
	if (bnd[i] <= 16.0 * eps * fabs(ritz[i])) {
	    neig += 1;
	    if (!*enough) *enough = endl < ritz[i] && ritz[i] < endr;
	}
    }   
    return;
}

/***********************************************************************
 *                                                                     *
 *				imtqlb()			       *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   imtqlb() is a translation of a Fortran version of the Algol
   procedure IMTQL1, Num. Math. 12, 377-383(1968) by Martin and 
   Wilkinson, as modified in Num. Math. 15, 450(1970) by Dubrulle.  
   Handbook for Auto. Comp., vol.II-Linear Algebra, 241-248(1971).  
   See also B. T. Smith et al, Eispack Guide, Lecture Notes in 
   Computer Science, Springer-Verlag, (1976).

   The function finds the eigenvalues of a symmetric tridiagonal
   matrix by the implicit QL method.


   Arguments 
   ---------

   (input)
   n      order of the symmetric tridiagonal matrix                   
   d      contains the diagonal elements of the input matrix           
   e      contains the subdiagonal elements of the input matrix in its
          last n-1 positions.  e[0] is arbitrary	             

   (output)
   d      contains the eigenvalues in ascending order.  if an error
            exit is made, the eigenvalues are correct and ordered for
            indices 0,1,...ierr, but may not be the smallest eigenvalues.
   e      has been destroyed.					    
   ierr   set to zero for normal return, j if the j-th eigenvalue has
            not been determined after 30 iterations.		    

   Functions used
   --------------

   UTILITY	fsign
   MISC		pythag

***********************************************************************/

static void imtqlb(long n, real d[], real e[], real bnd[])

{
    long last, l, m, i, iteration;

    /* various flags */
    long exchange, convergence, underflow;	

    real b, test, g, r, s, c, p, f;

    if (n == 1) return;
    ierr = 0;
    bnd[0] = 1.0;
    last = n - 1;
    for (i = 1; i < n; i++) {
	bnd[i] = 0.0;
	e[i-1] = e[i];
    }
    e[last] = 0.0;
    for (l = 0; l < n; l++) {
	iteration = 0;
	while (iteration <= 30) {
	    for (m = l; m < n; m++) {
		convergence = FALSE;
		if (m == last) break;
		else {
		    test = fabs(d[m]) + fabs(d[m+1]);
		    if (test + fabs(e[m]) == test) convergence = TRUE;
		}
		if (convergence) break;
	    }
	    p = d[l]; 
	    f = bnd[l]; 
	    if (m != l) {
		if (iteration == 30) {
		    ierr = l;
		    return;
		}
		iteration += 1;
		/*........ form shift ........*/
		g = (d[l+1] - p) / (2.0 * e[l]);
		r = pythag(g, 1.0);
		g = d[m] - p + e[l] / (g + fsign(r, g));
		s = 1.0;
		c = 1.0;
		p = 0.0;
		underflow = FALSE;
		i = m - 1;
		while (underflow == FALSE && i >= l) {
		    f = s * e[i];
		    b = c * e[i];
		    r = pythag(f, g);
		    e[i+1] = r;
		    if (r == 0.0) underflow = TRUE;
		    else {
			s = f / r;
			c = g / r;
			g = d[i+1] - p;
			r = (d[i] - g) * s + 2.0 * c * b;
			p = s * r;
			d[i+1] = g + p;
			g = c * r - b;
			f = bnd[i+1];
			bnd[i+1] = s * bnd[i] + c * f;
			bnd[i] = c * bnd[i] - s * f;
			i--;
		    }
		}       /* end while (underflow != FALSE && i >= l) */
		/*........ recover from underflow .........*/
		if (underflow) {
		    d[i+1] -= p;
		    e[m] = 0.0;
		}
		else {
		    d[l] -= p;
		    e[l] = g;
		    e[m] = 0.0;
		}
	    } 		       		   /* end if (m != l) */
	    else {

		/* order the eigenvalues */
		exchange = TRUE;
		if (l != 0) {
		    i = l;
		    while (i >= 1 && exchange == TRUE) {
			if (p < d[i-1]) {
			    d[i] = d[i-1];
			    bnd[i] = bnd[i-1];
			    i--;
			}
			else exchange = FALSE;
		    }
		}
		if (exchange) i = 0;
		d[i] = p;
		bnd[i] = f; 
		iteration = 31;
	    }
	}			       /* end while (iteration <= 30) */
    }				   /* end for (l=0; l<n; l++) */
    return;
}						  /* end main */

/***********************************************************************
 *                                                                     *
 *				imtql2()			       *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   imtql2() is a translation of a Fortran version of the Algol
   procedure IMTQL2, Num. Math. 12, 377-383(1968) by Martin and 
   Wilkinson, as modified in Num. Math. 15, 450(1970) by Dubrulle.  
   Handbook for Auto. Comp., vol.II-Linear Algebra, 241-248(1971).  
   See also B. T. Smith et al, Eispack Guide, Lecture Notes in 
   Computer Science, Springer-Verlag, (1976).

   This function finds the eigenvalues and eigenvectors of a symmetric
   tridiagonal matrix by the implicit QL method.


   Arguments
   ---------

   (input)                                                             
   nm     row dimension of the symmetric tridiagonal matrix           
   n      order of the matrix                                        
   d      contains the diagonal elements of the input matrix        
   e      contains the subdiagonal elements of the input matrix in its
            last n-1 positions.  e[0] is arbitrary	             
   z      contains the identity matrix				    
                                                                   
   (output)                                                       
   d      contains the eigenvalues in ascending order.  if an error
            exit is made, the eigenvalues are correct but unordered for
            for indices 0,1,...,ierr.				   
   e      has been destroyed.					  
   z      contains orthonormal eigenvectors of the symmetric   
            tridiagonal (or full) matrix.  if an error exit is made,
            z contains the eigenvectors associated with the stored 
          eigenvalues.					
   ierr   set to zero for normal return, j if the j-th eigenvalue has
            not been determined after 30 iterations.		    


   Functions used
   --------------
   UTILITY	fsign
   MISC		pythag

***********************************************************************/


static void imtql2(long nm, long n, real d[], real e[], real z[])

{
    long inbex, nnm, j2, last, l, m, i, k, iteration, convergence, underflow;
    real b, test, g, r, s, c, p, f;
    if (n == 1) return;
    ierr = 0;
    last = n - 1;
    for (i = 1; i < n; i++) e[i-1] = e[i];
    e[last] = 0.0;
    nnm = n * nm;
    for (l = 0; l < n; l++) {
	iteration = 0;

	/* look for small sub-diagonal element */
	while (iteration <= 30) {
	    for (m = l; m < n; m++) {
		convergence = FALSE;
		if (m == last) break;
		else {
		    test = fabs(d[m]) + fabs(d[m+1]);
		    if (test + fabs(e[m]) == test) convergence = TRUE;
		}
		if (convergence) break;
	    }
	    if (m != l) {

		/* set error -- no convergence to an eigenvalue after
		 * 30 iterations. */     
		if (iteration == 30) {
		    ierr = l;
		    return;
		}
		p = d[l]; 
		iteration += 1;

		/* form shift */
		g = (d[l+1] - p) / (2.0 * e[l]);
		r = pythag(g, 1.0);
		g = d[m] - p + e[l] / (g + fsign(r, g));
		s = 1.0;
		c = 1.0;
		p = 0.0;
		underflow = FALSE;
		i = m - 1;
		while (underflow == FALSE && i >= l) {
		    f = s * e[i];
		    b = c * e[i];
		    r = pythag(f, g);
		    e[i+1] = r;
		    if (r == 0.0) underflow = TRUE;
		    else {
			s = f / r;
			c = g / r;
			g = d[i+1] - p;
			r = (d[i] - g) * s + 2.0 * c * b;
			p = s * r;
			d[i+1] = g + p;
			g = c * r - b;

			/* form vector */
			for (k = 0; k < nnm; k += n) {
			    inbex = k + i;
			    f = z[inbex+1];
			    z[inbex+1] = s * z[inbex] + c * f;
			    z[inbex] = c * z[inbex] - s * f;
			} 
			i--;
		    }
		}   /* end while (underflow != FALSE && i >= l) */
		/*........ recover from underflow .........*/
		if (underflow) {
		    d[i+1] -= p;
		    e[m] = 0.0;
		}
		else {
		    d[l] -= p;
		    e[l] = g;
		    e[m] = 0.0;
		}
	    }
	    else break;
	}		/*...... end while (iteration <= 30) .........*/
    }		/*...... end for (l=0; l<n; l++) .............*/

    /* order the eigenvalues */
    for (l = 1; l < n; l++) {
	i = l - 1;
	k = i;
	p = d[i];
	for (j2 = l; j2 < n; j2++) {
	    if (d[j2] < p) {
		k = j2;
		p = d[j2];
	    }
	}
	/* ...and corresponding eigenvectors */
	if (k != i) {
	    d[k] = d[i];
	    d[i] = p;
	    for (j2 = 0; j2 < nnm; j2 += n) {
		p = z[j2+i];
		z[j2+i] = z[j2+k];
		z[j2+k] = p;
	    }
	}   
    }
    return;
}		/*...... end main ............................*/


/***********************************************************************
 *                                                                     *
 *				machar()			       *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   This function is a partial translation of a Fortran-77 subroutine 
   written by W. J. Cody of Argonne National Laboratory.
   It dynamically determines the listed machine parameters of the
   floating-point arithmetic.  According to the documentation of
   the Fortran code, "the determination of the first three uses an
   extension of an algorithm due to M. Malcolm, ACM 15 (1972), 
   pp. 949-951, incorporating some, but not all, of the improvements
   suggested by M. Gentleman and S. Marovich, CACM 17 (1974), 
   pp. 276-277."  The complete Fortran version of this translation is
   documented in W. J. Cody, "Machar: a Subroutine to Dynamically 
   Determine Determine Machine Parameters," TOMS 14, December, 1988.


   Parameters reported 
   -------------------

   ibeta     the radix for the floating-point representation       
   it        the number of base ibeta digits in the floating-point
               significand					 
   irnd      0 if floating-point addition chops		      
             1 if floating-point addition rounds, but not in the 
                 ieee style					
             2 if floating-point addition rounds in the ieee style
             3 if floating-point addition chops, and there is    
                 partial underflow				
             4 if floating-point addition rounds, but not in the
                 ieee style, and there is partial underflow    
             5 if floating-point addition rounds in the ieee style,
                 and there is partial underflow                   
   machep    the largest negative integer such that              
                 1.0+float(ibeta)**machep .ne. 1.0, except that 
                 machep is bounded below by  -(it+3)          
   negeps    the largest negative integer such that          
                 1.0-float(ibeta)**negeps .ne. 1.0, except that 
                 negeps is bounded below by  -(it+3)	       

***********************************************************************/

static void machar(long *ibeta, long *it, long *irnd, long *machep, long *negep)

{

    real beta, betain, betah, a2, b, ZERO, ONE, TWO, temp, tempa, temp1;
    long i, itemp;

    ONE = (real) 1;
    TWO = ONE + ONE;
    ZERO = ONE - ONE;

    a2 = ONE;
    temp1 = ONE;
    while (temp1 - ONE == ZERO) {
	a2 = a2 + a2;
	temp = a2 + ONE;
	temp1 = temp - a2;
    }

    b = ONE;
    itemp = 0;
    while (itemp == 0) {
	b = b + b;
	temp = a2 + b;
	itemp = (long)(temp - a2);
    }
    *ibeta = itemp;
    beta = (real) *ibeta;

    *it = 0;
    b = ONE;
    temp1 = ONE;
    while (temp1 - ONE == ZERO) {
	*it = *it + 1;
	b = b * beta;
	temp = b + ONE;
	temp1 = temp - b;
    }
    *irnd = 0; 
    betah = beta / TWO; 
    temp = a2 + betah;
    if (temp - a2 != ZERO) *irnd = 1;
    tempa = a2 + beta;
    temp = tempa + betah;
    if ((*irnd == 0) && (temp - tempa != ZERO)) *irnd = 2;


    *negep = *it + 3;
    betain = ONE / beta;
    a2 = ONE;
    for (i = 0; i < *negep; i++) a2 = a2 * betain;
    b = a2;
    temp = ONE - a2;
    while (temp-ONE == ZERO) {
	a2 = a2 * beta;
	*negep = *negep - 1;
	temp = ONE - a2;
    }
    *negep = -(*negep);

    *machep = -(*it) - 3;
    a2 = b;
    temp = ONE + a2;
    while (temp - ONE == ZERO) {
	a2 = a2 * beta;
	*machep = *machep + 1;
	temp = ONE + a2;
    }
    eps = a2;
    return;
}

/***********************************************************************
 *                                                                     *
 *                     store()                                         *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   store() is a user-supplied function which, based on the input
   operation flag, stores to or retrieves from memory a vector.


   Arguments 
   ---------

   (input)
   n       length of vector to be stored or retrieved
   isw     operation flag:
	     isw = 1 request to store j-th Lanczos vector q(j)
	     isw = 2 request to retrieve j-th Lanczos vector q(j)
	     isw = 3 request to store q(j) for j = 0 or 1
	     isw = 4 request to retrieve q(j) for j = 0 or 1
   s	   contains the vector to be stored for a "store" request 

   (output)
   s	   contains the vector retrieved for a "retrieve" request 

   Functions used
   --------------

   BLAS		dcopy

***********************************************************************/

static void store(long n, long isw, long j2, real *s)

{
    switch(isw) {
    case STORQ:	dcopy(n, s, 1, &a[(j2+MAXLL) * n], 1);
	break;
    case RETRQ:	dcopy(n, &a[(j2+MAXLL) * n], 1, s, 1);
	break;
    case STORP:	if (j2 >= MAXLL) {
	    fprintf(stderr,"store: (STORP) j2 >= MAXLL \n");
	    break;
	}
	dcopy(n, s, 1, &a[j2*n], 1);
	break;
    case RETRP:	if (j2 >= MAXLL) {
	    fprintf(stderr,"store: (RETRP) j2 >= MAXLL \n");
	    break;
	}
	dcopy(n, &a[j2*n], 1, s, 1);
	break;
    }
    return;
}
static real fsign(real a2,real b)
/************************************************************** 
 * returns |a| if b is positive; else fsign returns -|a|      *
 **************************************************************/ 
{

    if ((a2>=0.0 && b>=0.0) || (a2<0.0 && b<0.0)) return(a2);
    if ((a2<0.0 && b>=0.0) || (a2>=0.0 && b<0.0)) return(-a2);
    return a2;
}

static real ddmax(real a2, real b)
/************************************************************** 
 * returns the larger of two real precision numbers         *
 **************************************************************/ 
{

    if (a2 > b) return(a2);
    else return(b);
}

static real ddmin(real a2, real b)
/************************************************************** 
 * returns the smaller of two real precision numbers        *
 **************************************************************/ 
{

    if (a2 < b) return(a2);
    else return(b);
}

static long imin(long a2, long b)
/************************************************************** 
 * returns the smaller of two integers                        *
 **************************************************************/ 
{

    if (a2 < b) return(a2);
    else return(b);
}

static long imax(long a2,long b)
/************************************************************** 
 * returns the larger of two integers                         *
 **************************************************************/ 
{

    if (a2 > b) return(a2);
    else return(b);
}
/************************************************************** 
 * Constant times a vector plus a vector     		      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 
/************************************************************** 
 * Constant times a vector plus a vector     		      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 

static void daxpy (long n,real da,real *dx,long incx,real *dy,long incy)

{
    long i;

    if (n <= 0 || incx == 0 || incy == 0 || da == 0.0) return;
    if (incx == 1 && incy == 1) 
	for (i=0; i < n; i++) {
	    *dy += da * (*dx++);
	    dy++;
	}
    else {
	if (incx < 0) dx += (-n+1) * incx;
	if (incy < 0) dy += (-n+1) * incy;
	for (i=0; i < n; i++) {
	    *dy += da * (*dx);
	    dx += incx;
	    dy += incy;
	}
    }
    return;
}
/************************************************************** 
 * Function forms the dot product of two vectors.      	      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 

static real ddot(long n,real *dx,long incx,real *dy,long incy)

{
    long i;
    real dot_product;

    if (n <= 0 || incx == 0 || incy == 0) return(0.0);
    dot_product = 0.0;
    if (incx == 1 && incy == 1) 
	for (i=0; i < n; i++) dot_product += (*dx++) * (*dy++);
    else {
	if (incx < 0) dx += (-n+1) * incx;
	if (incy < 0) dy += (-n+1) * incy;
	for (i=0; i < n; i++) {
	    dot_product += (*dx) * (*dy);
	    dx += incx;
	    dy += incy;
	}
    }
    return(dot_product);
}
/************************************************************** 
 * function scales a vector by a constant.	     	      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 

static void datx(long n,real da,real *dx,long incx,real *dy,long incy)

{
    long i;

    if (n <= 0 || incx == 0 || incy == 0 || da == 0.0) return;
    if (incx == 1 && incy == 1) 
	for (i=0; i < n; i++) *dy++ = da * (*dx++);

    else {
	if (incx < 0) dx += (-n+1) * incx;
	if (incy < 0) dy += (-n+1) * incy;
	for (i=0; i < n; i++) {
	    *dy = da * (*dx);
	    dx += incx;
	    dy += incy;
	}
    }
    return;
}
/********************************************************************* 
 * Function sorts array1 and array2 into increasing order for array1 *
 *********************************************************************/

static void dsort2(long igap,long n,real *array1,real *array2)

{
    real temp;
    long i, j2, inbex;

    if (!igap) return;
    else {
	for (i = igap; i < n; i++) {
	    j2 = i - igap;
	    inbex = i;
	    while (j2 >= 0 && array1[j2] > array1[inbex]) {
		temp = array1[j2];
		array1[j2] = array1[inbex];
		array1[inbex] = temp;
		temp = array2[j2];
		array2[j2] = array2[inbex];
		array2[inbex] = temp;
	        j2 -= igap;
		inbex = j2 + igap;
	    }
	} 
    }
    dsort2(igap/2,n,array1,array2);
}
/************************************************************** 
 * Function interchanges two vectors		     	      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 

static void dswap(long n,real *dx,long incx,real *dy,long incy)

{
    long i;
    real dtemp;

    if (n <= 0 || incx == 0 || incy == 0) return;
    if (incx == 1 && incy == 1) {
	for (i=0; i < n; i++) {
	    dtemp = *dy;
	    *dy++ = *dx;
	    *dx++ = dtemp;
	}	
    }
    else {
	if (incx < 0) dx += (-n+1) * incx;
	if (incy < 0) dy += (-n+1) * incy;
	for (i=0; i < n; i++) {
	    dtemp = *dy;
	    *dy = *dx;
	    *dx = dtemp;
	    dx += incx;
	    dy += incy;
	}
    }
}

/***************************************************************** 
 * Function finds the inbex of element having max. absolute value*
 * based on FORTRAN 77 routine from Linpack by J. Dongarra       *
 *****************************************************************/ 

static long idamax(long n,real *dx,long incx)

{
    long ix,i,iMax;
    real dtemp, dMax;

    if (n < 1) return(-1);
    if (n == 1) return(0);
    if (incx == 0) return(-1);

    if (incx < 0) ix = (-n+1) * incx;
    else ix = 0;
    iMax = ix;
    dx += ix;
    dMax = fabs(*dx);
    for (i=1; i < n; i++) {
	ix += incx;
	dx += incx;
	dtemp = fabs(*dx);
	if (dtemp > dMax) {
	    dMax = dtemp;
	    iMax = ix;
	}
    }
    return(iMax);
}
/************************************************************** 
 * Function scales a vector by a constant.     		      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 

static void dscal(long n,real da,real *dx,long incx)

{
    long i;

    if (n <= 0 || incx == 0) return;
    if (incx < 0) dx += (-n+1) * incx;
    for (i=0; i < n; i++) {
	*dx *= da;
	dx += incx;
    }
    return;
}
/************************************************************** 
 * Function copies a vector x to a vector y	     	      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 

static void dcopy(long n,real *dx,long incx,real *dy,long incy)

{
    long i;

    if (n <= 0 || incx == 0 || incy == 0) return;
    if (incx == 1 && incy == 1) 
	for (i=0; i < n; i++) *dy++ = *dx++;

    else {
	if (incx < 0) dx += (-n+1) * incx;
	if (incy < 0) dy += (-n+1) * incy;
	for (i=0; i < n; i++) {
	    *dy = *dx;
	    dx += incx;
	    dy += incy;
	}
    }
    return;
}

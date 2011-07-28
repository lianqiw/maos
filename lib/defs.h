/*
  To be included in mat.c, cell.c and matbin.c
*/
#if MAT_VERBOSE == 1
#define matinfo(A...) {fprintf(stderr, A);}
#else
#define matinfo(A...)
#endif
#ifndef MAT_TYPE
#define MAT_TYPE
#ifndef USE_SINGLE
#ifndef USE_COMPLEX
//Double
#define X(A) d##A
#define Y(A) A
#define Z(A) d##A##_
#define T double
#define M_T M_DBL
#define M_TT M_DMAT
#define REAL(A) (A)
#define ABS(A) fabs(A)
#define SQRT(A) sqrt(A)
#define RANDU(A) randu(A)
#define RANDN(A) randn(A)
#define PRINT(A) printf("%10.3e",A);
#define CONJ(x) (x)
#define dot_do dotdbl
#else
//Double Complex
#define X(A) c##A
#define Y(A) c##A
#define Z(A) z##A##_ //blas/lapack convention
#define T dcomplex
#define M_T M_CMP
#define M_TT M_CMAT
#define REAL(A) creal(A)
#define ABS(A) cabs(A)
#define SQRT(A) csqrt(A)
#define RANDU(A) (randu(A)+I*randu(A))
#define RANDN(A) (randn(A)+I*randn(A))
#define PRINT(A) printf("(%10.3e %10.3eI)",creal(A),cimag(A));
#define CONJ(x) conj(x)
#define dot_do dotcmp
#endif
#else //#define USE_SINGLE
//Float
#ifndef USE_COMPLEX
#define X(A) s##A
#define Y(A) s##A
#define Z(A) s##A##_
#define T float
#define M_T M_FLT
#define M_TT M_SMAT
#define REAL(A) (A)
#define ABS(A) fabsf(A)
#define SQRT(A) sqrtf(A)
#define RANDU(A) randuf(A)
#define RANDN(A) randnf(A)
#define PRINT(A) printf("%10.3e",A);
#define CONJ(x) (x)
#define dot_do dotflt
#else
//Single Complex
#define X(A) z##A
#define Y(A) z##A
#define Z(A) c##A##_ //blas/lapack convention
#define T scomplex
#define M_T M_ZMP
#define M_TT M_ZMAT
#define REAL(A) crealf(A)
#define ABS(A) cabsf(A)
#define SQRT(A) csqrtf(A)
#define RANDU(A) (randuf(A)+I*randuf(A))
#define RANDN(A) (randnf(A)+I*randnf(A))
#define PRINT(A) printf("(%10.3e %10.3eI)",creal(A),cimag(A));
#define CONJ(x) conjf(x)
#define dot_do dotzmp
#endif//#define USE_COMPLEX
#endif//#define USE_SINGLE
#endif//#define MATTYPE
#define PMAT(A,pp) T (*restrict pp)[(A)->nx]=(void *)(A)->p
#define PCELL(M,P) X(mat)* (*restrict P)[(M)->nx]=(void*)(M)->p

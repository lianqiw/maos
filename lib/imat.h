#ifndef AOS_LIB_IMAT_H
#define AOS_LIB_IMAT_H
typedef struct imat{
    long nx;
    long ny;
    int *p;
}imat;

typedef struct icell{
    long nx;
    long ny;
    imat **p;
}icell;

imat* inew(long nx, long ny);
icell* icellnew(long nx, long ny);
void ifree(imat *A);
void icellfree(icell *A);
#endif

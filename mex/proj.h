#ifndef AOS_PROJ_H
#define AOS_PROJ_H
typedef struct RECTMAP_T{
    double *p;      /*OPD*/
    long nx,ny;     /*nx*ny, nx changes fastest, specifies x. ny specifies y*/
    double ox, oy;  /*origin*/
    double dx,dy;   /*spacing*/
    double h;       /*height of the grid*/
    double vx, vy;  /*wind speed. useful for atmospheric grid.*/
}RECTMAP_T;

/*map of locs*/
typedef struct LOCMAP_T{
    struct LOC_T *loc;
    long *p; //map is mxn
    double ox;
    double oy;
    int nx,ny;
    int npad;//just for checking. no need in computation.
}LOCMAP_T;

typedef struct LOC_T{
    double *locx;/*points to x coordinates of each point*/
    double *locy;
    long   nloc;/*number of points*/
    double dx; 
    LOCMAP_T *map;/*point to the map used for accphi only.*/
}LOC_T;
#endif

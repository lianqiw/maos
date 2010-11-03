#ifndef AOS_PROJ_H
#define AOS_PROJ_H
typedef struct rectmap_t{
    double *p;      /*OPD*/
    long nx,ny;     /*nx*ny, nx changes fastest, specifies x. ny specifies y*/
    double ox, oy;  /*origin*/
    double dx,dy;   /*spacing*/
    double h;       /*height of the grid*/
    double vx, vy;  /*wind speed. useful for atmospheric grid.*/
}rectmap_t;

/*map of locs*/
typedef struct locmap_t{
    struct loc_t *loc;
    long *p; /*map is mxn*/
    double ox;
    double oy;
    int nx,ny;
    int npad;/*just for checking. no need in computation.*/
}locmap_t;

typedef struct loc_t{
    double *locx;/*points to x coordinates of each point*/
    double *locy;
    long   nloc;/*number of points*/
    double dx; 
    locmap_t *map;/*point to the map used for accphi only.*/
}loc_t;
#endif

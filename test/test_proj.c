#include "../lib/aos.h"

/*
  Rectangular grid. (different X/Y Spacing)
 */

static void test_grid_proj(){
    dmat *junk=dread("M3_p.bin");
    dmat *X=dread("M3_x.bin");
    dmat *Y=dread("M3_y.bin");
    dmat *tmp=dread("M3_theta.bin");
    double bx=tmp->p[0];
    double by=tmp->p[1];
    dfree(tmp);
    rectmap_t *mapin=calloc(1, sizeof(rectmap_t));
    mapin->p=junk->p;
    mapin->ox=X->p[0];
    mapin->oy=Y->p[0];
    mapin->dx=X->p[1]-X->p[0];
    mapin->dy=Y->p[1]-Y->p[0];
    mapin->nx=X->nx*X->ny;
    mapin->ny=Y->nx*Y->ny;
    double d_m3_f=20.;/*from m3 to focus */
    double d_exitpupil_f=46.38661051;
    /*double d_m2_m3=23.59375000; */
    /*double d_m2_f=43.59375000; */
    double d_exitpupil_m3=d_exitpupil_f-d_m3_f;
    double r_exitpupil=1.546220350;
    double r_pupil=15;

    if(X->p[X->nx*X->ny-1]-X->p[0]-mapin->dx*X->nx*X->ny>1.e-10){
	error("X has non even spacing\n");
    }
    if(Y->p[Y->nx*Y->ny-1]-Y->p[0]-mapin->dy*Y->nx*Y->ny>1.e-10){
	error("Y has non even spacing\n");
    }
    free(junk);/*don't dfree */
    dfree(X); dfree(Y);
    /*direction of guide star */
    /*loc_t *loc2=mksqloc2(2000,2000,1./64.); */
    loc_t* loc2=locread("aper_locs.bin");
    
    dmat *amp=dread("aper_amp.bin");
    double *phi2=calloc(1, sizeof(double)*loc2->nloc);
    proj_rect_grid(mapin,M_PI*0.75,M_PI*0.5,
		   loc2,-r_exitpupil/r_pupil,r_exitpupil/r_pupil,
		   amp->p,phi2,-2,d_exitpupil_f,d_exitpupil_m3,-bx,-by);
    /*drawopdamp("test_proj",loc2,phi2,amp,"phi"); */
    
    writedbl(phi2,loc2->nloc,1,"phi");
    /*locwrite(loc2,"loc"); */
}
int main(){
    test_grid_proj();
}

#ifndef AOS_LIB_GENOTF_H
#define AOS_LIB_GENOTF_H
#include "type.h"
#include "loc.h"
void genotf(cmat **otf,    /**<The otf array for output*/
	    LOC_T *loc,    /**<the common aperture grid*/
	    const double *amp,   /**<The amplitude map of all the (sub)apertures*/
	    const double* opdbias,  /**<The static OPD bias. */
	    const double *area, /**<normalized area of the (sub)apertures*/
	    double thres,  /**<The threshold to consider a (sub)aperture as full*/
	    double wvl,    /**<The wavelength. only needef if opdbias is not null*/
	    double dtheta, /**<Sampling of PSF.*/
	    double r0,     /**<Fried parameter*/
	    double l0,     /**<Outer scale*/
	    long ncompx,    /**<Size of OTF*/
	    long ncompy,    /**<Size of OTF*/
	    long nsa,       /**<Number of (sub)apertures*/
	    long pttr,      /**<Remove piston/tip/tilt*/
	    long nthread    /**<Number of threads*/
	    );
#endif

/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#ifndef SKYC_SYPES_H
#define SKYC_SYPES_H
#include "parms.h"
/**
   \file skyc/types.h

   Run time structs
   Renamed _T to _S to prevent from conflict with maos/ folder.
*/
/**
   Detector transfer function */
typedef struct DTF_S{ 
    cmat *nominal;     /**<The FFT of the pixel functions*/
    dsp *si;           /**<The pixel selection*/
    cmat *U;           /**<Special frequency vector along x*/
}DTF_S;
/**
   Struct for POWFS*/
typedef struct POWFS_S{
    int     ipowfs;    /**<Which ipowfs this is*/
    DTF_S  *dtf;       /**<array of dtf for each wvl*/
    loc_t  *loc;       /**<explicit pts in a regular grid. */
    dmat   *amp;       /**<amplitude map defined on loc*/
    double *locxamp;   /**<dot(loc->locx,amp);*/
    double *locyamp;   /**<dot(loc->locy,amp);*/
    loc_t  *saloc;     /**<subaperture location*/
    double dxwvf;      /**<sampling for the wvf.*/
    int    nxwvf;      /**<number of points each subaps*/
}POWFS_S;
/**
   Struct for pixel intensity statistics*/
typedef struct PISTAT_S{
    dcell *psf;        /**<short exposure PSF*/
    dcell *neaspec;    /**<noise equivalent angle due to spec noise*/
    dcell *i0;         /**<normalized pixel intensity, doesn't contain siglev*/
    dcell *gx;         /**<gradient of i0*/
    dcell *gy;         /**<gradient of i0*/
    dcell *i0s;        /**<i0 multiplied with signal level*/
    dcell *gxs;        /**<gradient of i0s*/
    dcell *gys;        /**<gradient of i0s*/
    dcell **mtche;     /**<matched filter estimator*/
    dcell *sanea;      /**<noise equivalent angle, including all effect*/
    dcell *sanea0;     /**<noise equivalent angle due to photon noise alone*/
    dmat  *scale;      /**<extra scaling factor due to bicubic spline on strehl*/
    dmat  *snr;        /**<signal to noise ratio*/
}PISTAT_S;
/**
   Struct for WFS*/
typedef struct WFS_S{
    int ipowfs;        /**<type of powfs for each wfs.*/
    int istar;         /**<index to star;*/
    double thetax;     /**<location of wfs*/
    double thetay;     /**<location of wfs*/
    dmat *mags;        /**<magnitude in each band*/
    dmat *siglev;    /**<signal level at maos.dt sampling for each wvl*/
    double siglevtot;  /**<sum of siglev*/
    double bkgrnd;     /**<background level per pixel at maos.dt*/
    ccell **wvfout;    /**<complex wavefront output from maos run*/
    dmat *ztiltout;   /**<ztilt out from maos run*/
    dmat *goff;        /**<gradient offset for ncpa calibration*/
    PISTAT_S *pistat;  /**<information about pixel intensities. first got from
			  star_t, then updated.*/
}WFS_S;

/**
   Data for each available star.
*/
typedef struct STAR_S{
    double thetax;     /**<star location*/
    double thetay;     /**<star location*/
    double scbcs;      /**<scale factor to match strehl for bi-cubic spline*/
    dmat  *mags;       /**<star magnitude at each wavelength*/
    /*The following is for each ipowfs */
    dcell  *siglev;    /**<signal level of size npowfs, nwvl*/
    dmat  *siglevtot;  /**<total signal level of size npowfs*1*/
    dmat  *bkgrnd;     /**<pixel background of size npowfs*1*/
    int   *use;        /**<whether this star is used in phy simulations, of size npowfs*1*/
    PISTAT_S *pistat;  /**<pixel intensity statstics npowf*1*/
    dcell *g;          /**<gradient operator of size npowfs*1*/
    ccell ***wvfout;   /**<complex wavefront output from maos run, of size npowfs*1*/
    dcell *ztiltout;  /**<ztilt output from maos run, of size npowfs*1*/
    dcell *goff;       /**<gradient offset for NCPA calibration.*/
    int   nstep;       /**<number of time steps available.*/
    dmat* idtrat;       /**<prefered dtrat for minimum snr*/
}STAR_S;
/**
   asterism dependent data.*/
typedef struct ASTER_S{
    int iaster;        /**<index of aster*/
    int nstep;         /**<number of time steps available*/
    int tsa;           /**<total number of subapertures.*/
    int nwfs;          /**<number of WFS.*/
    int use;           /**<use this asterism for further physical optics simulation.*/
    WFS_S *wfs;        /**<list of wfs for this asterism*/
    dcell *g;          /**<NGS mode to grad operator*/
    dmat *gm;          /**<matrix version of g.*/

    /*The following are for each dtrat */
    dcell *pgm;        /**<mode reconstructor */
    dcell *gain;       /**<type II gain vector*/
    dcell **neam;        /**<measurement error covariance matrix (full matrix with values in diagonal)*/
    dcell *sigman;     /**<NGS, TT noise propagated from WFS measurement noise.*/
    dmat *res_ws;      /**<residual windshake after servo rejection.*/
    dmat *res_ngs;     /**<residual ngs mode error after servo. */
    int mdtrat;        /**<dtrat of minimum rms in OL estimation.*/
    double mresol;     /**<miminum rms on OL restimation.*/
    rand_t rand;  /**<random stream*/
    int idtratmin;     /**<minimum index of dtrat allowed*/
    int idtratmax;     /**<maximum index of dtrat allowed*/
    kalman_t**kalman;
    dmat *idtrats;
    dmat *dtrats;
    long *ngs;         /**<number of gradients for each wfs*/
}ASTER_S;
/**
   A few simulation parameters.*/
typedef struct SIM_S{
    const PARMS_S *parms; /**<input parameters*/
    POWFS_S *powfs;    /**<structs about POWFS*/
    dcell *stars;      /**<randomly generated star fields.*/
    dmat *mideal;      /**<ideal NGS modes*/
    dmat *mideal_oa;   /**<ideal NGS modes on axis (dot product only)*/
    double rmsol;       /**<open loop error*/
    STATUS_T *status;  /**<to report status to the scheduler*/
    STAR_S *star;      /**<STAR_S*/
    int iseed;         /**<Current seed index*/
    int seed_maos;     /**<Current MAOS seed to read in PSF*/
    int nstep;         /**<Number of steps*/
    rand_t rand;       /**<Random stream*/
    dmat *res;         /**<residual error. 5*nsky. 
			  - Total
			  - NGS mode
			  - TT mode
			  - Residual wind shake
			  - Residual focus error
		       */
    dmat *res_oa;      /**<residual error on axis 5*nsky. in the same format as res*/
    dmat *res_geom;      /**<residual error estimated from servo analysis. in the same foramt as res*/
    dcell *mres;       /**<residual NGS modes. 5*nsky*/
    dcell *sel;        /**<selected stars for each star field.*/
    dmat *fss;         /**<optimal sampling freq for each star field*/
    dcell *gain;       /**<optimal gains for each star field.*/
    dmat *psd_ngs;     /**<PSD of NGS(tip/tilt + plate scale) modes*/
    dmat *psd_tt;      /**<PSD of Tip/tilt modes*/
    dmat *psd_ps;      /**<PSD of plate scale modes*/
    dmat *psd_focus;   /**<PSD of focus error due to atmosphere and sodium layer*/
    dmat *gain_x;      /**<The sampled sigma in m2*/
    dcell **gain_ngs;  /**<TypeII gain and result of different sigman and different dtrat*/
    dcell **gain_ps;   /**<TypeII gain and result of different sigman and different dtrat*/
    dcell **gain_tt;   /**<TypeII gain and result of different sigman and different dtrat*/
    dcell **gain_focus;/**<TypeII gain and result of different sigman and different dtrat*/
    dcell ***bspstrehl;/**<Coeffiecients of bicubic spline fit of Strehl on the grid.*/
    dmat *bspstrehlxy; /**<Coordinate of the grid.*/
    int isky;          /**<current star field being evaluated*/
    int isky_start;    /**<first star field to evaluate*/
    int isky_end;      /**<last star field to evalaute (exclusive)*/
    double tk_0;       /**<initial star time*/
    pthread_mutex_t mutex_status;/**<mutex for status reporting*/
    dmat *sdecoeff;    /**<sde coefficient*/
    dcell *psdi;        /**<PSD of each mode computed from time series*/
    dcell **nonlin;     /**<nonlinearity*/
    dmat *neaspec_dtrats;
}SIM_S;
#endif

/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
typedef struct dtf_s{ 
    cmat *nominal;     /**<The FFT of the pixel functions*/
    dsp *si;           /**<The pixel selection*/
    cmat *U;           /**<Special frequency vector along x*/
}dtf_s;
/**
   Struct for POWFS*/
typedef struct powfs_s{
    int     ipowfs;    /**<Which ipowfs this is*/
    dtf_s  *dtf;       /**<array of dtf for each wvl*/
    loc_t  *loc;       /**<explicit pts in a regular grid. */
    dmat   *amp;       /**<amplitude map defined on loc*/
	dmat   *sacent;	   /**<Amplituded weighted average of locx and locy */
    loc_t  *saloc;     /**<subaperture location*/
    real 	dxwvf;      /**<sampling for the wvf.*/
    int   	nxwvf;      /**<number of points each subaps*/
}powfs_s;
/**
   Struct for pixel intensity statistics*/
typedef struct pistat_s{
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
    dmat  *scale;      /**<extra scaling factor for pistat and psf due to bicubic spline on strehl*/
    dmat  *snr;        /**<signal to noise ratio*/
}pistat_s;
/**
   Struct for WFS*/
typedef struct wfs_s{
    int ipowfs;        /**<type of powfs for each wfs.*/
    int istar;         /**<index to star;*/
    real thetax;     /**<location of wfs*/
    real thetay;     /**<location of wfs*/
    dmat *mags;        /**<magnitude in each band*/
    dmat *siglev;    /**<signal level at maos.dt sampling for each wvl*/
    real siglevtot;  /**<sum of siglev*/
    real bkgrnd;     /**<background level per pixel at maos.dt*/
    ccell **wvfout;    /**<complex wavefront output from maos run*/
    dmat *ztiltout;   /**<ztilt out from maos run*/
    dmat *goff;        /**<gradient offset for ncpa calibration*/
    pistat_s *pistat;  /**<information about pixel intensities. first got from
			  star_t, then updated.*/
}wfs_s;

/**
   Data for each available star.
*/
typedef struct star_s{
    real thetax;     /**<star location*/
    real thetay;     /**<star location*/
    dmat  *mags;       /**<star magnitude at each wavelength*/
    /*The following is for each ipowfs */
    dcell  *siglev;    /**<signal level of size npowfs, nwvl*/
    dmat  *siglevtot;  /**<total signal level of size npowfs*1*/
    dmat  *bkgrnd;     /**<pixel background of size npowfs*1*/
    int   *use;        /**<whether this star is used in phy simulations, of size npowfs*1*/
    pistat_s *pistat;  /**<pixel intensity statstics npowf*1*/
    dcell *g;          /**<gradient operator of size npowfs*1*/
    ccell ***wvfout;   /**<complex wavefront output from maos run, of size npowfs*1*/
    dcell *ztiltout;  /**<ztilt output from maos run, of size npowfs*1*/
    dcell *goff;       /**<gradient offset for NCPA calibration.*/
    int   nstep;       /**<number of time steps available.*/
    dmat* minidtrat;   /**<prefered dtrat for minimum snr*/
}star_s;
/**
   asterism dependent data.*/
typedef struct aster_s{
    int iaster;        /**<index of aster for this sky*/
	int iaster_all;	   /**<index of aster for all sky */
    int nstep;         /**<number of time steps available*/
    int tsa;           /**<total number of subapertures.*/
    int nwfs;          /**<number of WFS.*/
    int use;           /**<use this asterism for further physical optics simulation.*/
    wfs_s *wfs;        /**<list of wfs for this asterism*/
    dcell *g;          /**<NGS mode to grad operator*/
    dmat *gm;          /**<matrix version of g.*/
    lmat *mdirect;     /**<in multirate mode, directly output such modes for slower mode*/
    /*The following are for each dtrat */
    dcell *pgm;        /**<mode reconstructor */
    dcell *gain;       /**<type II gain vector*/
    //dccell *neam;        /**<measurement error covariance matrix (full matrix with values in diagonal)*/
    dcell *sigman;     /**<NGS, TT noise propagated from WFS measurement noise.*/
    dmat *res_ws;      /**<residual windshake after servo rejection.*/
    dmat *res_ngs;     /**<residual ngs mode error after servo. */
    real minest;    /**<miminum rms on servo restimation.*/
    rand_t rand;       /**<random stream*/
    kalman_t **kalman;	/**<Kalman filter for LQG controller. */
    lmat *idtrats;	   /**<idtrat of each wfs*/
    lmat *dtrats;	   /**<dtrat of each wfs */
    long *ngs;         /**<number of gradients for each wfs*/
    dcell *phyres;      /**<Wavefront variance result from physical optics simulations.*/
    dcell *phymres;     /**<Residual modes from physical optics simulation*/
	real minphyres;	/**<Best PO performance (at idtratphy) */
	int idtratphy;		/**<idtrat for the best PO performance */
	int idtratmin;     /**<minimum index of dtrat allowed*/
    int idtratmax;     /**<maximum index of dtrat allowed*/
	int idtratest;     /**<dtrat of minimum rms in OL estimation.*/
}aster_s;
/**
   A few simulation parameters.*/
typedef struct sim_s{
    const parms_s *parms; /**<input parameters*/
    powfs_s *powfs;    /**<structs about POWFS*/
    dcell *stars;      /**<randomly generated star fields.*/
    dmat *mideal;      /**<ideal NGS modes*/
    dmat *mideal_oa;   /**<ideal NGS modes on axis (dot product only)*/
    real varol;      /**<open loop error*/
    status_t *status;  /**<to report status to the scheduler*/

    rand_t rand;       /**<Random stream used to generate stars and time series from psds.*/
    dmat *res;         /**<final residual error for saving. 5*nsky. Total WFE, NGS mode, TT mode, wind shake (if separately estimated), estimated. */
    //dmat *res_est;     /**<residual error estimated from servo analysis. in the same foramt as res*/
    dmat *res_aster;   /**<Performance and parameter of all asterisms evaluated.*/
    
    dcell *mres;       /**<residual NGS modes. 5*nsky*/
    dcell *gain;       /**<optimal gains for each star field.*/
    dcell *psds;        /**<PSD of All (gsplit=0) or Each mode (gsplit=1)*/
    dmat *gain_x;      /**<The sampled sigma in m2*/
    dcccell *gain_pre; /**<TypeII gain and result of different sigman and different dtrat*/
    dcell ***bspstrehl;/**<Coeffiecients of bicubic spline fit of Strehl on the grid.*/
    dmat *bspstrehlxy; /**<Coordinate of the grid.*/

    real tk_0;       /**<initial star time*/
    pthread_mutex_t mutex_status;/**<mutex for status reporting*/
    dmat *sdecoeff;    /**<sde coefficient*/
    dccell *nonlin;     /**<nonlinearity*/
    dmat *neaspec_dtrats;
	unsigned int res_iaster;   /**<counter for res_aster*/
	unsigned int isky; /**<current star field being evaluated*/
    unsigned int isky_start;    /**<first star field to evaluate*/
    unsigned int isky_end;      /**<last star field to evalaute (exclusive)*/
    unsigned int isky_print;    /**<printed isky*/
    unsigned int iseed;         /**<Current seed index*/
    int seed_maos;     /**<Current MAOS seed to read in PSF*/
    unsigned int nstep;         /**<Number of steps*/
}sim_s;
#endif

/*
  Copyright 2009-2018 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#ifndef SKYC_PARMS_H
#define SKYC_PARMS_H
#include <photon.h>
/**
   \file skyc/parms.h
   Parameters for simulation.
*/
/** Parameters for skycoverage simulation relates to high order
   pre-processing. Some of the parameters are to be exported by MAOS. The rest
   are to be supplied by the end user.  */

typedef struct MAOS_S{
    double r0z;    /**<Fried parameter, used to read in PSD information.*/
    double dt;     /**<simulation period time. 1/800*/
    double zadeg;  /**<Zenith angle in degree*/
    double za;     /**<Derived: zenith angle*/
    double hc;     /**<Height conjugation of upper DM*/
    double hs;     /**<Height of LGS*/
    double D;      /**<Diameter of telescope*/
    double *wvl;   /**<Wavelength vector*/
    int nwvl;      /**<Number of wavelength*/
    int npowfs;    /**<Number of types of wfs.*/
    int nmod;      /**<Number of modes controlled by LOW NGS WFS;*/
    int *seeds;    /**<List of seed for maos PSF outputs.*/
    int nseed;     /**<Number of seeds*/
    double *wddeg; /**<wind direction list for each layer. just to propagate to
		      maos for integrated simulation.*/
    int nwddeg;    /**<number of wddeg*/

    /*The following are vectors for each powfs. */
    int *nsa;      /**<total number of subapertures*/
    int *msa;      /**<order of correction in 1-d*/
    int *ncomp;    /**<size of saved psf/wvf*/
    int *embfac;   /**<embeding factor of saved psf. determines the PSF sampling.*/
    double *dxsa;  /**<subaperture size*/
    char **fnwfsloc; /**<file name for WFS grid*/
    char **fnwfsamp; /**<file name for WFS amptitude defined on grid*/
    char **fnsaloc;  /**<file name for subaperture location coordinate*/
    char *fnmideal;  /**<file name for ideal NGS modes*/
    char *fnmidealp; /**<file name for ideal NGS mode computation for each direction*/
    int evlindoa;    /**<index of the On axis evaluation point.*/
    dmat *mcc;       /**<NGS modal cross coupling matrix.*/
    dmat *mcc_oa;    /**<NGS modal cross coupling matrix for on axis only.*/
    dmat *mcc_tt;    /**<NGS tip/tilt modal cross coupling matrix*/
    dmat *mcc_oa_tt; /**<NGS tip/tilt modal cross coupling matrix for on axis only, 2x2*/
    dmat *mcc_oa_tt2;/**<modal cross coupling matrix between NGS and tip/tilt for on axis only, 2x5*/
    double ngsgrid;  /**<spacing of NGS grid in as*/
    int nstep;       /**<Number of time steps (sim.nsim) in MAOS simulation*/
    int ahstfocus;   /**<1: The focal plate scale mode does not include global focus*/
    int mffocus;     /**<1: maos already does sodium tracking*/
    int indfocus;    /**<Index of focus mode. 0: disable*/
    int indps;       /**<Index of plate scale mode. 0: disable*/
    int indastig;    /**<Index of astigmatism mode. 0: disable*/
    char *fnrange;   /**<sodium range time series. forward to maos.*/
}MAOS_S;

/**
   Contains parameters for skycoverage.
*/
typedef struct SKYC_S{
    ZB_S zb;         /**<Sky background and zero magnitude flux*/
    int verbose;     /**<be verbose at output.*/
    //int reest;       /**<reestimate the error after gain estimation.*/
    int dbg;         /**<save intermediate information for debugging.*/
    int dbgsky;      /**<only run this sky frame if not -1*/
    int dbgaster;    /**<only run this asterism if not -1*/
    int keeporder;   /**<1: Keep order of stars as input is skyc.star.*/
    int interpg;     /**<Interpolate gain based in dtrat and signan*/
    int save;        /**<save information for latter running in MAOS.*/
    int maxstar;     /**<maximum number of (brightest) stars for each star field.*/
    int maxaster;    /**<maximum number of best asterism for final processing.*/
    int maxdtrat;    /**<maximum number of dtrats to try for each asterism.*/
    int nthread;     /**<number of threads*/
    int start;       /**<start number of sky field*/
    int nsky;        /**<number of sky field. 500.*/
    int seed;        /**<seed for generating asterism and random numbers.*/
    int navg;        /**<Number of simulations to average*/
    int noisefull;   /**<use full noise insteaded servo filtered noise to regularize. =0*/
    int psd_scale;   /**<scale the PSD to equal to open loop error*/
    int noisy;       /**<noise simulations.*/
    int ttfbrightest;/**<make ttf the brightest always.*/
    int bspstrehl;   /**<Do bicubic spline interpolation on strehl*/
    int npowfs;      /**<number of powfs, has to match MAOS_S.npowfs*/
    double lat;      /**<Galactic latitude*/
    double lon;      /**<Galactic longitude.*/
    double catscl;   /**<Scale the catlog star count*/
    double patfov;   /**<Patrol FoV in arcsec (diameter)*/
    double minrad;   /**<minimum radius of the stars (keep out of center science region)*/
    double imperrnm; /**<Implementation error in nm in the center.*/
    double imperrnmb;/**<Implementation error slopt in nm:
			imperrnm(theta)=sqrt(imperrnm^2+imperrnmb^2*theta^2). The
			location in fov, theta, is normalized by patfov/2*/
    /** The following are vectors for each powfs; */
    int *nwfsmax;    /**<maximum # of wfs for each powfs*/
    int  nwfstot;    /**<Sum of nwfsmax*/
    double *pixtheta;/**<size of WFS ccd pixel.*/
    double *pixblur; /**<blurring of pixels.*/
    int *pixpsa;     /**<number of detector pixels in each direction per sa*/
    int *pixguard;   /**<guard window size in each direction.*/
    double *pixoffx; /**<pixel offset along x in units of pixel*/
    double *pixoffy; /**<pixel offset along y in units of pixel*/
    double keepout;  /**<NGS probe keep out range in arcsec.*/
    double rne;      /**<detector read out noise in electron. -1 to use the formula.*/
    double excess;   /**<detector excess photon noise factor*/
    dmat *rnefs;     /**<derived, actual read out noise, may be frame rate dependent.*/
    double *telthruput;/**<Telescope throughput at each wvl*/
    double *qe;      /**<quantum efficiency at each wvl.*/
    double wvlmean;  /**<Mean wavelength*/
    dmat *wvlwt;     /**<Weight for each wavelength*/
    int ngsalign;    /**<align NGS to the grid.*/
    int limitnstep;  /**<Limit the number of steps in each simulation*/
    int evlstart;    /**<time step to start evaluate performance*/
    int phystart;    /**<time step to start physical optics*/
    int mtchcr;      /**<constraint in matched filter*/
    int mtch;        /**<matched filter*/
    int neaaniso;     /**<use additional measurement error caused by
			 anisoplanatism in regularization.*/
    int neanonlin;   /**<use additional measurement error caused by WFS nonlinearity*/
    int ndtrat;      /**<number of dtrat*/
    dmat *dtrats;     /**<ratio between NGS and LGS WFS sampling period*/
    dmat *dtrats_mr; /**<For multirate*/
    double *fss;     /**<sampling frequency at each dtrat*/
    int servo;       /**<servo type of NGS LOOP. 2: type II*/
    int ngain;       /**<Number of parameters for gain*/
    int gsplit;      /**<use separate gains for tip/tilt and plate scale modes*/
    dmat *psd_ngs;   /**<PSD of NGS(tip/tilt + plate scale) modes*/
    dmat *psd_tt;    /**<PSD of Tip/tilt modes*/
    dmat *psd_ps;    /**<PSD of plate scale modes*/
    dmat *psd_ws;    /**<PSD of windshake*/
    dmat *psd_focus; /**<PSD of residual focus due to sodium layer variation, propagated by LGS WFS*/
    double zc_f;     /**<focus zoom corrector frequency*/
    double zc_zeta;  /**<focus zoom corrector dampling */
    double na_alpha; /**<sodium PSD parameter. PSD is beta*f^alpha*/
    double na_beta;  /**<sodium PSD parameter. PSD is beta*f^alpha*/
    char *fnrange;   /**<Change of sodium height for LGS. in meter. Replaes na_alpha, na_beta*/
    char *stars;     /**<file name of not NULL to load stars from*/
    int addws;       /**<add wind shake time series to simulation*/
    double pmargin;  /**<phase margin of type II*/
    int psdcalc;     /**<Calculate PSD from time series*/
    char **fnpsf1;   /**<file name for additional otf to be interpolated and
			multiplied to dtfq. 2 columns. first column is coordinate
			of otf, and second column of value.*/
    double sdetmax;  /**<tmax for SDE fitting*/
    int multirate;   /**<Each OIWFS can run at different dtrat*/
    dmat* snrmin;   /**<Minimum SNR to determine minimum dtrat. SNR computed as pixtheta/nea*/
    int usephygrad;  /**<1: Use physical optics grad instead of ztilt*/
    int estimate;    /**<1: Estiamte performance only, without time domain simulation*/
}SKYC_S;
/**
   Parameters for skycoverage.
 */
typedef struct PARMS_S{
    MAOS_S maos;     /**<parameters exported by maos*/
    SKYC_S skyc;     /**<parameters supplied by user*/
    int *fdlock;     /**<Records the fd of the seed lock file. if -1 will skip the seed*/
}PARMS_S;
/**
   ARG_S is used for command line parsing.
*/
typedef struct ARG_S{
    int detach;      /**<Detach from the command line and run in background*/
    int override;    /**<Override result even if Res?_?.done exist*/
    int force;       /**<For start, bypassing scheduler*/
    int nthread;     /**<Number of threads*/
    char *dirout;    /**<Result output directory*/
    char *conf;      /**master .conf file. nfiraos.conf by default. -c to change*/
    char *confcmd;
}ARG_S;

PARMS_S *setup_parms(const ARG_S *arg);
void free_parms(PARMS_S *parms);
#endif

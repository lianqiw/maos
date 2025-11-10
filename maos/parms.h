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
#ifndef __AOS_PARMS_H__
#define __AOS_PARMS_H__
#include "../lib/aos.h"
/**
   \file maos/parms.h

   Configuration parameters that remain constant during simulation.
*/
#define MAX_AMPS 5
extern real TOMOSCALE;
/*********************************************/
/*             Structs for parms             */
/*********************************************/
/**
   contains input parameters for the atmospheric turbulence.
*/
typedef struct atm_cfg_t{
    real r0z;   /**<r0 at zenith*/
    real r0;    /**<derived from r0z for zenith angle za*/
    real dx;    /**<sampling of turbulence screens*/
    real hmax;  /**<maximum in ht*/
    dmat *L0;     /**<outer scale. One number or one per layer*/
    dmat *r0logpsds; /**<[alpha beta]: temporal PSD of log(r0) is beta*f^alpha. f is in hz.*/
    dmat *r0logpsdt; /**<[alpha beta]: spatial  PSD of log(r0) is beta*f^alpha. f is in m.*/
    dmat *ht;   /**<height of each layer*/
    dmat *wt;   /**<weight of each layer (relative strength of \f$C_n^2\f$)*/
    dmat *ws;   /**<wind speed of each layer*/
    dmat *wddeg;/**<wind direction of each layer*/
    dmat *size; /**<size of atm in meter, [0 0]: automatic*/
    lmat *overx;  /**<maximum pixel distance in x direction the beam can be without wrapping*/
    lmat *overy;  /**<maximum pixel distance in y direction the beam can be without wrapping*/
    int nps;      /**<number of phase screens*/
    int wdrand;   /**<randomize wind direction*/
    int iground;  /**<index into the ground layer*/
    lmat *ipsr;    /**<corresponding reconstruction layer*/
    int nx;       /**<turbulence screen size along x*/
    int ny;       /**<turbulence screen size along y*/
    lmat *nxn;    /**<minimum turbulence screen size along to cover meta pupil*/
    int nxnmax;   /**<max of nxn*/
    int method;   /**<0: FFT Von Karman. 1: FFT Biharmonic. 2: Fractal method.*/
    int frozenflow;  /**<frozen flow. automatic if closeloop=1*/
    int ninit;    /**<Initial size of the screen in fractal method. >=2*/
    int share;    /**<0: disable sharing of atmosphere using file backend*/
    int r0evolve; /**<Evolve r0 according to r0logpsd*/
	int dtrat;    /**<Only used if atm are loaded with frames that need to playback in time. */
	int interp;   /**<Interpolation method when atm.dtrat>0. 0: stepwise(no interpolation), 1:linear, 2:sin^2.*/
}atm_cfg_t;
/**
   contains input parameters for the atmospheric reconstruction.  */
typedef struct atmr_cfg_t{
    real r0z;   /**<r0 at zenith*/
    real r0;    /**<derived from r0z for zenith angle za*/
    real L0;    /**<outer scale*/
    real hs;    /**<height of the high order guide star. derived*/
    real hmax;  /**<maximum of ht*/
    dmat *ht;   /**<height of each layer*/
    dmat *wt;   /**<weight of each layer (relative strength of \f$C_n^2\f$)*/
    real dx;    /**<baseline sampling (when os=1). matches to high order wfs.*/
    lmat *indps;   /**<Mapping atmr.ps to atm.ps*/
    lmat *os;      /**<over sampling factor of xloc over actuator spacing */
    int nps;      /**<number of phase screens*/
}atmr_cfg_t;
/**
   contains input parameters about the aperture, like the diameter,
   amplitude map, etc */
typedef struct aper_cfg_t{
    real d;     /**<Telescope aperture diameter*/
    real din;   /**<Telescope inner blocking diameter*/
    real rot;	/**<pupil rotation*/
	map_t *amp; /**<Aperture amplitude map from aper.fnamp config. */
    dmat *misreg;  /**<Calibrated misregistration of the telescope pupil. 2x1*/
    dmat *misregu; /**<Uncalibrated misregistration of the telescope pupil. */
    char *pupmask;/**<The pupil cold stop*/
}aper_cfg_t;
/**
   contains input parameters for laser launch telescope
*/
typedef struct llt_cfg_t{
    real d;      /**<LLT clear aperture diameter*/
    real widthp; /**<Gaussian beam width percentage of d*/
    real focus;  /**<RMS focus error in nm of LLT.*/
    real ttrat;  /**<Ratio of uplink jitter to science jitter due to M2 windshake.*/
    real fcfsm;  /**<corner frequency for offloading FSM to a common path pointing mirror in LLT. 0: disabled*/
	real dhs;    /**<Spacing of sublayers to simulate for LGS*/
    char *ttpsd;   /**<PSD of uplink beam jitter*/
    char *fnrange; /**<File contains range to sodium layer*/
    char *fnprof;  /**<File contains sodium profile*/
    char *fnprep;/**<File contains sodium profiled used for computing i0. if NULL: equal to fnprof*/
    char *fnsurf;  /**<Pupil Surface OPD error*/
    char *fnamp;   /**<Pupil amplitude map. overrides widthp*/
    dmat *ox;    /**<location x of LLT center wrt telescope aperture center*/
    dmat *oy;    /**<see ox.*/
    dmat *misreg; /*beam to pupil misregistration (centering error)*/

    int ttfr;      /**<Remove piston/tip/tilt and focus (if = 2) from ncpa*/
    int colprep;   /**<starting column to use in fn for ETF in preparation of matched filter*/
    int colsim;    /**<starting column to use in fn for ETF in simulation*/
    int coldtrat;  /**<change to next sodium profile during simulation every coldtrat time step*/
    int nhs;       /**<Number of sublayer to simulate for LGS*/
	//Moved from dbg_cfg_t
	int na_smooth;   /**<1: smooth sodium profile to coarser grid before computing etf*/
	int na_interp;   /**<1: Interpolate sodium profile and use FFT to build etf. 0: direct sum, slow*/
	int na_fit_maxit;/**<Number of iterations. 0: auto, 1 for CMF, 3 for COG. see wfsgrad.c*/
	real na_fit_svdthres;/**<threshold for SVD inverse in sodium fitting.*/
	real na_fit_dh;  /**<sampling in height in sodium fitting*/
	real na_thres;   /**<altitude error threshold to move trombone, in unit of meter.*/

    //Computed
    lmat *i;        /**<Index into llt for this iwfs.*/
    int nllt;       /**<number of launch telescopes in this powfs*/
    //real lpfsm;    /**<Low pass filter for LLT FSM offloading*/
    real epfsm;    /**<Integrator gain for LLT FSM offloading*/
} llt_cfg_t;
/**
* parameters for dithering
* */
typedef struct dither_cfg_t{
	int mode;
	real amp; /**<Dither amplitude.*/
	real gpll;/**<Gain of phase locked loop*/
	real gog; /**<Gain for updating optical gain for cog*/
	real gdrift;/**<Gain for drift control*/
	real glpf; /**<LPF gain for i0,gx,gy update (matched filter)*/
	int npoint;/**<Number of points in each dither peroid (4)*/
	int pllskip;/**<Skip WFS frames for uplink loop to stable*/
	int pllrat; /**<Number of WFS frames for updating PLL.*/
	int ogskip; /**<Number of WFS frames to skip before computing averaged images*/
	int ograt;  /**<Number of WFS frames to update pixel processing algorithm (MF/CoG)*/
	int ogsingle;/**<*Force using single gain update (when dither==1 for SHWFS)*/
} dither_cfg_t;

/**
   contains input parameters for each type of wfs (powfs).
*/
typedef struct powfs_cfg_t{
    dmat *wvl;     /**<list of wavelength in ascending order.*/
    dmat *wvlwts;  /**<weights for each wavelength. can be overriden by wfs.wvlwts.*/
	real wvlmean;  /**<Average wavelength*/
    char *saloc;   /**<saloc override file*/
	real misregx;  /**<misregistration wrt telescope pupil. shift along x */
	real misregy;  /**<misregistration wrt telescope pupil. shift along y */
	real misregc;  /**<misregistration wrt telescope pupil. camera is rotated CCW.*/
	real misregrmax;  /**<Maximum misregistration for all wfs in this powfs. */
    char *amp;     /**<amplitude override file*/
    char *piinfile;/**<input averaged pixel intensities for matched filter. NULL to disable*/
    char *sninfile;/**<Speckle noisy input file. NULL to disable. not used*/
    real hs;       /**<height of guide star*/
    real hc;       /**<conjugation height of WFS pupil*/
    real saat;     /**<subaperture area (normalized) threshold to drop subaperture.*/
    real safill2d; /**<subaperture lenslet throughgput. value is used  to alter amplitude map*/
    real saspherical;/**<Subaperture spherical aberration in nm RMS at best focus.*/
    real safocuspv;   /**<Subaperture focus error in nm p/v*/
    char  *neareconfile;/**<File contains noise equivalent angle in radian^2. Contains cell array of nwfsx1.*/
    real nearecon; /**<NEA used in reconstruction in milli-arcsecond, sim.dtref integration time. Will be scaled by powfs.dtrat and subaperture area before use.*/
    real neasim;   /**<NEA used in simulation. -1 to use nearecon*/
    char*  neasimfile;/**<read NEA used in simulation from file. Defined at sim.dt sampling rate, in radian. neasim must be -1*/
    real neaextra; /**<Extra NEA to add in quadrature to the NEA determined by matched filter or CoG*/
    real neamin;   /**<Minimum NEA to limit the NEA determined by matched filter or CoG*/
    real bkgrnd;   /**<background in electron per pixel per LGS frame*/
    real bkgrndc;  /**<How much of the background in bkgrnd can be calibrated out. depends on variability.*/
    char *bkgrndfn;/**<file contains sky background/rayleigh scatter input for each subaperture in each wfs. */
    char *bkgrndfnc;/**<How much of the background in bkgrndfn can be calibrated out. depends on variability.*/
    dmat *qe;      /**<File containing matrix of pixpsax*pixpsay specifying QE of each pixel. To simulate PCCD non uniform response*/
    real rne;      /**<read out noise in electron per pixel per frame*/
    real pixblur;  /**<pixel bluring due to leakage. relative to pixel size.*/
    real dsa;      /**<Size of subaperture in 1 dimension*/
    real dx;       /**<sampling of opd points in each subaperture. usually matches atmosphere sampling for LGS. may be coraser for NGS.*/
    real pixtheta; /**<size of pixel pitch along x/y or azimuthal if radial ccd. Converted to radian from user input*/
    real radpixtheta; /**<size of pixel pitch along radial direction. -1 for square pixel*/
    real fieldstop;/**<size of field stop in arcsec.*/
	real astscale; /**<Scale wfs.thetax and wfs.thetay.*/
    real pixoffx;  /**<offset of image center from center of detector*/
    real pixoffy;  /**<see pixoffx*/
    real sigscale; /**<scale the signal level for simulation.*/
    real siglev;   /**<signal level at dtrat=1. will be override by wfs.siglev is specified.*/
    dmat* siglevs;  /**<in array format. 1x1 or nwfsx1, scaled from wfs.siglev by dtrat. */
    real sigrecon; /**<signal level for NEA computation*/
    struct llt_cfg_t *llt;/**<configuration for LLT*/
    char* fnllt;   /**<filename of LLT configuration. empty means no llt.*/
    pywfs_cfg_t *pycfg; /**<Set only for Pyramid WFS.*/
    char *pywfs; /**<Pyramid WFS configuration*/
    int type;      /**<WFS type: 0: SHWFS, 1:Pyramid WFS*/
    int step;      /**<frame to start using WFS*/
    int trs;       /**<tip/tilt removal flag. True for LGS, False for NGS*/
    int frs;       /**Global focus removal flag. Optional for LGS, False for NGS*/
    int lo;        /**<whether this is a low order wfs. False for LGS, True for NGS*/
    int skip;      /**<skip in high order tomography, for split tomo (derived parameter)*/
    int psol;      /**<Compute pseudo open loop gradients (derived parameter)*/
    lmat *wfs;     /**<array of wfs belongs to this powfs*/
    lmat *wfsr;    /**<array of reconstruction wfs belongs to this powfs*/
    lmat *wfsind;  /**<wfsind[iwfs] gives the index of the wfs in this powfs group*/
    int nwfs;      /**<number of wfs belonging to this powfs*/
    int nwfsr;     /**<number of wfs for reconstruction belonging to this powfs*/
    int neaphy;    /**<use nea from physical optical precomputation in geometric simulations.*/
    int phyusenea; /**<force using supplied noise equivalent angle in physical optics simulations*/
    int order;     /**<order of wavefront sensing along one dimension.*/
    int pixpsa;    /**<number of detector pixels along x/y (or azimuthal if radial CCD).*/
    int radpix;    /**<number of detector pixels along radial direction if radial CCD*/
    int radgx;     /**<1: gx/gy is along R/A coordinate. Only valid if radpix is set*/
    int notf;      /**<PSF is extended to this size before FFT into OTF.  0 for automatic*/
    int embfac;    /**<Embed subaperture atm OPD before fft. set to 2.*/
    int nwvl;      /**<Number of wavelength. 1 for LGS, 2 for NGS J+H sensing.*/
    int gtype_sim; /**<wfs type if not using physical optics in reconstruction.
		       - 0: geometric
		       - 1: ztilt.*/
    int gtype_recon;/**<wfs type if not using physical optics in simulation.
		       - 0: geometric
		       - 1: ztilt.*/
    int phytype_recon;/**<physical optics type for reconstruction. 1: mtch, 2: tcog, 3: MAP*/
    int phytype_sim;  /**<physical optics type for simulation. -1 to follow phytype_recon*/
    int phytype_sim1; /**<Save phytype_sim initial value*/
    int phytype_sim2; /**<physical optics type after dithering update. -1 to follow phytype_sim*/
    int phystep;    /**<frame to start using physical optics.
		       -  0: means from frame 0.
		       - >0: need to compute GS0 to calculate geometric optics
		       - -1: never, doesn't need to compute DTF
		    */
    int usephy;     /**<whether physical optics is used at all during simulation.(derived parameter)*/
    real r0;      /**<Fried parameter  for matched filter generation. Uses atm.r0, atm.L0 is not set*/
    real L0;      /**<Outerscale for matched filter generation. Uses atm.r0, atm.L0 is not set*/
    int mtchcr;     /**<use constrained matched filter (0: disable, 1: both axis. 2: radial/x only, 3: az/y only)*/
    int mtchcpl;    /**<use coupling between r/a measure error. useful for LGS with x-y ccd.*/
    int mtchstc;    /**<shift peak in the time averaged short exposure PSF to center using fft.*/
    int sigmatch;   /**<scale subaperture image to have the same intensity as i0. Keep false.*/
    int mtchadp;    /**<Using adaptive matched filter. When the number of pixels
		       in the image brighter than half maximum is more than this
		       value, use constraint. introduced on 2011-02-21.*/
    int mtchfft;    /**<Compute gx, gy using i0 with FFT derivative instead of PSF.*/
    real cogthres; /**<CoG threshold, relative to max(im)*/
    real cogoff;   /**<CoG offset to remove, relative to max(im). */
	real cogmask;  /**<Using this threshold to setup mask for CoG based on i0*/
	int cogcalib;   /**<Calibrate CoG gradient offset using i0. */
	dmat* ncpa;     /**<Description of NCPA; 2xn; first row is rms in meter, second row is zernike mode or negative for power law.*/
    int needGS0;    /**<need to compute GS0 (derived parameter)*/
    int noisy;      /**<noisy or not during *simulation* */
    int ncpa_method;/**<Method to correct ncpa.
		       - 0: do nothing.
		       - 1: derive gradient offset from OPD
		       - 2: derive gradient offset(CoG) or matched filter from i0 with bias OPD*/
    int pistatout;  /**<output time averaged short exposure image. 1: shift to center, 2: do not shift*/
    int pistatstart;/**<time step to compute pistatout*/
    int pistatstc;  /**<1: shift to center using fft method. 0: use geometric gradients.*/
    int psfout;     /**<output time history of low order wfs PSF. never do this for LGS.*/
    int dtrat;      /**<ratio of sample period over sim.dt. Note that sim.dt is assumed to be the fast, high-order loop (e.g. LGS) and the slow or low-order TT WFS cannot run faster. powfs.dtrat must be an integer value.*/
    int idtrat;     /**<Index of dtrat into parms->sim.dtrats*/
    int i0scale;    /**<scale i0 to matched subaperture area.*/
    int moao;       /**<index into MOAO struct. -1: no moao*/

    int i0save;     /**<Save time averaged subaperture images.*/
    char *i0load;   /**<load i0,gx,gy from this folder.*/
    real gradscale; /**<Scale CL gradients. For testing*/
	//char* fndither; /**<Configuration for dither*/
    int dither;     /**<Turn on/off dithering to update centroid gain or matched filter*/
    int dither_mmd;/**<Enable multi-mode dithering. Only effective when recon.modal is set*/
	real dither_amp; /**<Dither amplitude.*/
	real dither_gpll;/**<Gain of phase locked loop*/
	real dither_gog; /**<Gain for updating optical gain for cog*/
	real dither_gdrift;/**<Gain for drift control*/
	real dither_glpf; /**<LPF gain for i0,gx,gy update (matched filter)*/
	int dither_npoint;/**<Number of points in each dither peroid (4)*/
	int dither_pllskip;/**<Skip WFS frames for uplink loop to stable*/
	int dither_pllrat; /**<Number of WFS frames for updating PLL.*/
	int dither_ogskip; /**<Number of WFS frames to skip before computing averaged images*/
	int dither_ograt;  /**<Number of WFS frames to update pixel processing algorithm (MF/CoG)*/
	int dither_ogsingle;/**<*Force using single gain update (when dither==1 for SHWFS)*/
    //options for zoom corrector
    int zoomdtrat;   /**<Separated from llt.coldtrat to work with statc sodium profile*/
    int zoomshare;   /**<1: All LGS share the same trombone*/
    real zoomgain; /**<gain of the trombone controller*/
    real zoomgain_drift; /**<gain for the trombone controller with i0 drift input*/
    int ng;         /**<number of gradients per subaperture. 2 for SHWFS. >2 for raw PWFS*/

    //options for FSM
    real apfsm;     /**<servo coefficient for for LGS uplink pointing loop.*/
    real epfsm;     /**<error gain for uplink pointing*/
    real alfsm;      /**<Additional latency (*sim.dt) of the uplink loop*/
    real zetafsm;    /**<Damping of FSM modeled as second harmonic oscillater (SHO).*/
    real f0fsm;      /**<Resonance frequency of FSM (SHO). 0: infinite.*/
    int idealfsm;    /**<ideal compensation for uplink pointing*/
    int commonfsm;   /**<Make FSM common for each powfs (LLT). Keep at 0. */
}powfs_cfg_t;
/**
   contains input parmaeters for each wfsr for reconstruction
*/
typedef struct wfsr_cfg_t{
	real thetax;  /**<x direction*/
	real thetay;  /**<y direction*/
	real misregx; /**<misregistration wrt telescope pupil: shift along x */
	real misregy; /**<misregistration wrt telescope pupil. shift along y */
	real misregc; /**<misregistration wrt telescope pupil. clocking error.*/
	//real hc;    /**<conjugation height of WFS pupil is wfs.hc=powfs.hc+wfs.delta_hc (input)*/
	real hs;      /**height of star is wfs.hs=powfs.hs+wfs.delta_hs (input)*/
	int powfs;    /**<powfs type*/
}wfsr_cfg_t;
/**
   contains input parmaeters for each wfs
*/
typedef struct wfs_cfg_t{
    dmat *wvlwts; /**<Weights of signal value for each wavelength. if not specified in config, will use powfs.wvlwts*/
    dmat *sabad;  /**<coordinate of bad subaperture due to bad detector or lenslet array.*/
    real thetax;  /**<x direction*/
    real thetay;  /**<y direction*/
	real misregx; /**<misregistration wrt telescope pupil: shift along x */
	real misregy; /**<misregistration wrt telescope pupil. shift along y */
	real misregc; /**<misregistration wrt telescope pupil. clocking error.*/
    real hc;      /**<conjugation height of WFS pupil is wfs.hc=powfs.hc+wfs.delta_hc (input)*/
    real hs;      /**height of star is wfs.hs=powfs.hs+wfs.delta_hs (input)*/
    real siglev; /**<Total signal value for all wavelength. if not specified in config, will use powfs.siglev*/
    real sigsim; /**<Signal value used for simulation. (derived parameter)*/
    real fitwt;  /**<Include wfs in fitting directions if corresponding wfs[iwfs].fitwt is greater than 0*/
    int powfs;   /**<powfs type*/
}wfs_cfg_t;
/**
   contains input parameters for each deformable mirror.
*/
typedef struct dm_cfg_t{
    real guard;   /**<extra DM actuator rings outside of aper.d*/
    dmat *stroke; /**<Stroke of DM (surface). OPD goes to \f$\pm\f$ stroke. nactx2 array: min and max per actuator$*/
    real iastroke;/**<Inter actuator stroke (surface)*/
    dcell *strokescale;   /**< describes polynomials that convert opd to voltage (first cell), and voltage to opd
			    	* (second cell). The two operations has to be strict inverse of each other*/
	real dratio;     /**<telescope diameter to DM diameter ratio (beam angle magnification factor)*/
    real ht;      /**<height conjugation range*/
    real dx;      /**<actuator separation along x (derived from order)*/
    real ar;      /**<[in] aspect ratio: dy/dx*/
    real dy;      /**<actuator separation along y (derived from dx and ar*/
    real offset;  /**<Center-most actuator offset from origin
		       - =0 means there is a act on center.
		       - 1/2 means no one in the center.*/
    real iac;     /**<If !=0: use cubic influence function with this Inter-Actuator Coupling coefficient.*/
    real histbin; /**<The bin width for histogram.*/
	real nmod;    /**<Maximum number of modes to control in modal controller*/

    int histn;      /**<Number of bins in histogram.*/
    int hist;       /**<Compute histogram of commands of each actuator*/
    int order;      /**<Nominal order of the DM within telescope clear subaperture*/
    int isground;   /**<Is this DM the ground DM (derived)*/
    dmat *actfloat; /**<floating actuators. nx2 coordinate*/
    dmat *actstuck; /**<stuck actuators. nx2 coordinate.*/

    real hyst;       /**<The hysteresis amount (ratio)*/
    real hyst_alpha; /**<The DM hysteresis model alpha parameter*/
    real hyst_stroke;/**<The surface stroke that the hysteresis is measured at*/
}dm_cfg_t;
/**
   contarins input parameters all evaluation directions.  */
typedef struct evl_cfg_t{
    dmat *thetax;   /**<x Coordinate of evaluation directions*/
    dmat *thetay;   /**<y Coordinate of evaluation directions*/
    dmat *wt;       /**<weight of each direction*/
    dmat *wvl;      /**<wavelength for PSF and strehl computation*/
    dmat *hs;       /**<height of each science object*/
	const char **wvlname; /**<Common name of each wavelength. */
    real dx;        /**<sampling of aperture for evaluation*/
    int nwvl;       /**<Number of wavelength*/
    lmat *psf;      /**<1: participate in psf evaluation.*/
    lmat *psfr;     /**<1: participate in psf reconstruction telemetry*/
    int npsf;       /**<how many directions we compute psf for*/
    int rmax;       /**<max radial mode for performance evaluation.
                    - 0: piston only
                    - 1: piston/tip/tilt.*/
    int nmod;       /**<Number of modes. derived from rmax. (nmax+1)*(nmax+2)/2*/

    int psfol;      /**<compute Open loop PSF.
                    - 1: on axis only.
                    - 2: all directions and average them.*/
    int psfhist;    /**<output history of the psf (a lot of storage)*/
    int psfmean;    /**<output time averaged psf*/
    int cov;        /**<save covairance of science OPD ,every this time step, for directions where evl.psf is 1*/
    int opdmean;    /**<save science OPD time average every `evlopdmean` time steps*/
    lmat *pttr;     /**<remove p/t/t from psf. 1 number for each evl.*/

    int psfisim;    /**<time step to start psfmean.*/
    lmat *psfsize;  /**<save this number of pixels of the center of the psf. 1
			number for each wvl.*/
    lmat *psfgridsize;/**<grid size for FFT to generate PSF. Becareful about FFT
			speed and enough padding. Determines the sampling of the
			generated PSF. 0 or negative for automatic. 1 number for
			each wvl.*/
    int nevl;       /**<Number of evaluation directions. (derived)*/
    int tomo;       /**<evaluate tomography performance.*/
    int indoa;      /**<index of the on axis evluation point.*/
    int moao;       /**<index into MOAO struct. -1: no MOAO*/
	int split;	    /**<evaluate split tomography low order*/
}evl_cfg_t;

/**
   contains input parameters for wavefront tomography.
*/
typedef struct tomo_cfg_t{
    real tikcr;    /**<tikhonov regularization.*/
    real minwt;    /**<minimum layer weight allowed. if less than this will force to this.*/
    real iac;      /**<!=0: use cubic influence function with this Inter-Actuator Coupling coefficient.*/
    real cxxscale; /**<scale the Cxx^-1 term.*/
    real svdthres; /**<Threshold in SVD inversion*/
    int square;      /**<use square/rectangular grid instead of tighter irregular grid*/
    int cone;        /**<use cone coordinate in xloc: keep true*/
    int cxxalg;      /**<method to compute Cxx^-1. 0: bihormonic approx. 1: inverse psd. 2: fractal*/
    int guard;       /**<guard rings of reconstruction grid ploc and xloc*/
    int pos;         /**<over sampling factor of ploc over actuator spacing*/
    int nxbase;      /**<Each layer xloc grid size is tomo.os*tomo.nxbase is not zero. same for ploc.*/
    int piston_cr;   /**<single point piston constraint. */

    int ahst_wt;     /**<Weight used to compute low order model removal in AHST
			1: remove effect on NGS WFS (not good if WFS is outside of science FoV)
			2: remove effect on Science
			3: Identity weighting (bad)
		     */
    int ahst_idealngs;/**<ideal correction on NGS modes. For skycoverage preprocessing.*/
    int ahst_focus;   /**<1: Make magnification mode free of focus in science (only effective when sim.mffocus=1*/
	int ahst_keepfocus;/**<keep LGS focus in ngs mode removal*/
    int alg;         /**<Tomography algorithm to solve the linear equation.\todo implement BGS, MG
			0: Cholesky direct solve for the large matrix.  (CBS)
			1: CG or PCG.
			2: SVD or EVD: Eigen value decompsition
		     */
    int bgs;         /**<1: use BGS, block Gaussia Seidel then use alg to solve each block.*/
    int precond;     /**<Tomography preconditioner.
			0: No preconditioner.             (CG)
			1: Fourier Domain Preconditioner. (FDPCG)
		     */
    int maxit;       /**<max iterations. Usually 30 for CG, 3 for FDPCG in
			closed loop warm restart. x10 in open loop*/
    int cgwarm;      /**<Warm restart in CG. */
    int assemble;    /**<force assemble tomography matrix in CG*/
    int predict;     /**<test predictive control.*/
    int ninit;       /**<like atm.ninit, the initial screen to generate from covariance directly*/
    int splitlrt;    /**<1: low rank terms also in LHS.*/
}tomo_cfg_t;
/**
   contains input parameters for deformable mirror fitting.
*/
typedef struct fit_cfg_t{
    dmat *thetax;  /**<x Coordinate of DM fitting directions. */
    dmat *thetay;  /**<y Coordinate of DM fitting directions. */
    dmat *wt;      /**<weight of each direction*/
    dmat *hs;      /**<height of target in each direction*/
    real tikcr;    /**<tikhonov regularization*/
    real svdthres; /**<Threshold in SVD inversion*/
    real actthres; /**<Threshold for slaving value of weakly coupled actuators*/
    real actthres2;/**<Threshold for reducing jump across weakly coupled actuators*/
    int actslave;    /**<Enable slaving for non-active actuators. Useful in CBS method*/
    int actextrap;   /**<extrapolate actuator results to non-active actuators.*/
    int nfit;        /**<Number of DM fit directions */
    int lrt_piston;  /**<Piston constraint low rank term in fit coefficient matrix*/
    int lrt_tt;      /**<differential tip/tilt constraint on two DMs or tt on upper dms.*/
    int alg;         /**<Fitting algorithm to solve the linear equation.
			0: Cholesky direct solve for the large matrix.  (CBS)
			1: CG or PCG.
			2: SVD or EVD: Eigen value decompsition
		     */
    int bgs;         /**<1: use BGS, block Gaussia Seidel then use alg to solve each block.*/
    int precond;     /**<Preconditioner. Not available.*/
    int maxit;       /**<max iterations. Usually 4 for CG*/
	int guard;       /**<guard rings of reconstruction grid ploc*/
    int square;      /**<using square grid on DM and ploc.*/
    int assemble;    /**<force assemble fit matrix in CG*/
    int pos;         /**<over sampling of floc over aloc. for fitting. normally equal to tomo.pos*/
    int indoa;       /**<Index of on axis point.*/
    int cachedm;     /**<Cache DM command in intermediate plane*/
    int cachex;      /**<Cache X (xloc) in intermediate plane*/
    int cgwarm;      /**<Warm restart in CG. */
}fit_cfg_t;
/**
   contains input parameters for the least square reconstructor.
*/
typedef struct lsr_cfg_t{
    real tikcr;    /**<tikhonov regularization*/
    real svdthres; /**<Threshold in SVD inversion*/
    real actthres; /**<Threshold for slaving value of weakly coupled actuators*/
    real actthres2;/**<Threshold for reducing jump across weakly coupled actuators*/
    char  *fnreg;    /**<File containing a regularization term to add to LL.M*/
    int actextrap;   /**<extrapolate actuator results to non-active actuators .*/
    int actslave;    /**<Enable slaving for non-active actuators. Useful in CBS method*/
    int splitlrt;    /**<1: low rank terms also in LHS.*/
    int bgs;         /**<1: use BGS, block Gaussia Seidel then use alg to solve each block.*/
    int alg;         /**<algorithm to solve the linear equation.
			0: Cholesky direct solve for the large matrix.  (CBS)
			1: CG or PCG.
			2: SVD or EVD: Eigen value decompsition
		     */
    int maxit;       /**<max iterations. Usually 30 for CG*/
    int cgwarm;      /**<Warm restart in CG. */

}lsr_cfg_t;
/**
   contains input parameters for wavefront reconstruction.
*/
typedef struct recon_cfg_t{
    real psdservo_gain; /**<Gain used to update servo parameter*/
    real poke;       /**<How much WFE (meter) to apply to OPD for computing experimental interaction matrix*/

    int alg;         /**<algorithm for reconstruction. 0: MVR. 1: LSR. moved from sim.recon*/
    int glao;        /**<whether we are in GLAO mode where all WFS in each powfs are averaged*/
    int split;       /**<split reconstruction/tomography type.
						- 0: integrated tomography
						- 1: adhoc split tomography
						- 2: minimum variance split tomography (only valid if recon.alg=0)*/
    int modal;       /**-2: emulate zonal, -1: zernike, 0: zonal, 1: KL modes*/
    int psol;        /**<Use pseudo open loop gradients for wavefront reconstruction*/
    int mvm;         /**<Use the various algorithms recon.alg to assemble a control
						matrix to multiply to gradients to get DM commands. If
						the algorithm needs PSOL gradient, we will have an
						auxillary matrix to multiply to the DM actuators and
						subtract from the result.*/

    int psd;         /**<Flag: compute PSDs of DM error signal averaged over aperture and field points (m^2/Hz).*/
    int psddtrat_hi; /**<how many time step to sample for PSD computation.*/
    int psddtrat_lo; /**<how many time step to sample for low order PSD computation.*/

    int psdnseg;     /**<how many overlapping partitions of the time history to compute PSD.*/
    int twfs_rmin; 	 /**<minimum zernike order (inclusive)*/
    int twfs_rmax;	 /**<maximum zernike order (inclusive)*/
    int twfs_radonly;/**<1: radial only, 0: all modes*/
	int petal;  	 /**1: enable petal mode control*/
	int petaldtrat;	 /**<how many time steps to average for petaling mode control.*/
	int petalstep; 	 /**<simulation step to enable petal mode control*/
	int petalnpsf;   /**<number of pixels (each dimension) used for petal reconstruction*/
	int petaltt;    /**<whether include t/t modes in the reconstruction*/
    char **distortion_dm2wfs; /**<Distortion from DM to each WFS model used in reconstruction. Affects GA*/
    char **distortion_dm2sci; /**<Distortion from DM to each science model used in reconstruction. Affects HA*/
    char **distortion_tel2wfs;/**<Distortion from Telescope to each WFS model used in reconstruction. Affects HXW*/
}recon_cfg_t;
/**
   contains input parameters for simulation, like loop gain, seeds, etc.
*/
typedef struct sim_cfg_t{
    real dt;         /**<sampling period (s) for simulation.*/
    real dtref;      /**<sampling period (s) for setting siglev or nearecon.*/
    real za;         /**<zenith angle in radian*/
    real htel;       /**<Height of telescope. Used to adjust sodium profile range*/
    int start;       /**<time step to start simulation. 0*/
    int end;         /**<time step to stop simulation. exclusive*/
    int pause;       /**<Pause at the end of every time step*/
    lmat *seeds;      /**<simulation seeds*/
    int nseed;       /**<How many simulation seed*/
    int closeloop;   /**<closed loop or open loop*/
    dmat *wspsd;     /**<Telescope wind shake PSD input. Nx2. First column is
						freq in Hz, Second column is PSD in rad^2/Hz.*/
    int wsseq;       /**<sequence of wind shake time series.*/
    /*control */
    dmat *aphi;      /**<servo coefficient for high order dm.  A is command. e is
						error signal. at time step n, the command is updated by
						A(n)=A(n-1)*apdm(0)+A(n-2)*ap(1)+...+e(n-2)*ep
						*/
    dmat *ephi;      /**<error gain for DM commands (high order)*/

    real f0dm;      /**<Natural frequency of the DMs.*/
    real zetadm;    /**<Damping of the DMs.*/
    //real f0tt;      /**<Natural frequency of tip/tilt mirror. We do not implement SHO for tip/tilt as it offloads from DM which does high speed tip/tilt */
    //real zetatt;    /**<Damping frequency of tip/tilt mirror.*/
    dmat *aplo;      /**<servo coefficient for ngs modes.*/
    dmat *eplo;      /**<error gain for NGS modes (low order)*/

    real alhi;     /**<Additional latency (*sim.dt) of the high order loop besides 2 cycle delay.*/
    real allo;        /**<Additional latnecy (*sim.dt) of the low order loop*/

    real aptwfs;   /**<Twfs reference vector servo coefficient.*/
    real eptwfs;   /**<Twfs reference vector servo gain.*/
    real eptsph;   /**<Twfs reference vector servo gain for spherical mode*/
    real fcttm;    /**<cross over frequency of tip/tilt split. 0 to disable ttm.*/
    real fcfocus;  /**<cross-over frequency of the focus LPF.*/
    real fov;      /**<User specified fov diameter*/
	real foveff;   /**<The effective fov diameter */
    int focus2tel;   /**<Offload focus to telescope*/
    real epfocus2tel;/*Gain for telescope focus control*/
    int mffocus;     /**<method for focus blending between LGS and LO NGS
			- 0: no focus blending.
			- 1: Focus blending using CL gradients, for each LGS independently.
			- 2: Focus blending using CL gradinets, for common LGS focus only (not preferred).
		     */
    int cachedm;     /**<cache dm shape on fine sampled grid matched WFS or Science grid*/
    int fuseint;     /**<fuse the high and low order integrators in split tomography */
    int skysim;      /**<1: we are doing skycoverage preprocessing*/
    int evlol;       /**<evaluate open loop error only*/
    int noatm;       /**<disable atmosphere*/
    int idealtomo;   /**<Use downsampled turbulence directly as tomography outout (to evaluate fitting error without tomography effect)*/
    int psfr;        /**<do PSF reconstruction telemetry*/
    int ecnn;        /**<Calculate WF covariance due to WFS noise cov Cnn.*/
    int wfsalias;    /**<Study the wfs aliasing effect by projecting turbulence
			onto the NULL space of DM.*/
    int idealwfs;    /**<Generates ideal WFS by sensing turbulence with DM range.*/
    int idealevl;    /**<Evaluate performance within DM range.*/

    /* A few derived parameters*/
    real dtlo;     /**<low order wfs sampling period*/
    real dthi;     /**<high order wfs sampling period*/
    int dtrat_hi;    /**<ratio of sampling period over clock of high order wfs*/
    int dtrat_lo;    /**<highest dtrat of the lower order loop.*/
    int dtrat_lo2;   /**<lowest dtrat of the lower order loop.*/
    int dtrat_lof;   /**<lowest dtrat of the lower order focus loop.*/
    int dtrat_skip;  /**<dtrat (over sim.dt) for frame drop. Be careful when powfs.dtrat is not one.*/
    int noisy_hi;    /**<whether high order WFS is noisy*/
    int noisy_lo;    /**<whether low order WFS is noisy*/
    real lpfocushi;	/**<derived: lpfocus=2*pi*fc*sim.dthi*/
    real lpfocuslo;	/**<derived: lpfocus=2*pi*fc*sim.dtlo*/
    real lpttm;    	/**<los path filter for ttm. derived: lpttm=2*pi*fcttm*sim.dt*/

    int dmclip;      /**<derived: Need to clip actuator stroke*/
    int dmclipia;    /**<derived: Need to clip inter-actuator stroke*/
    int dmproj;      /**<derived: Need to projection atmosphere onto DMspace. */
    int mvmport;     /**<Non zero: specify which port does the MVM server run on and connect to it for MVM reconstruction.*/
    char *mvmhost;   /**<Which host does the MVM server run*/
    int mvmsize;     /**<number of gradients to send each time. 0 is all.*/
    int mvmngpu;     /**<number of GPUs to use in server*/
	char *dmadd;     /**<Containing dm vector to simulate turbulence (added to integrator output).
			 		 	 It should be cell array (time steps) of cell arry (DMs) of vectors. Can be empty*/
}sim_cfg_t;
typedef struct ncpa_cfg_t{
    dmat *thetax; 	/**<Coordinate for NCPA calibration (arcsec)*/
    dmat *thetay; 	/**<Coordinate for NCPA calibration (arcsec)*/
    dmat *wt;   	/**<Weight for each point.*/
    dmat *hs;   	/**<Height of star.*/
	int calib;  	/**<calibrate NCPA. 1: with all DMs. 2: with only ground DM.*/
	int ttr;    	/**<Remove average t/t from NCPA for WFS. Equivalent as repositioning WFS. default 1.*/
	int rmsci;  	/**<1: do not include calibration residual in science path.*/
	int preload; 	/**<preload integrator with DM sys flat*/
    int ndir;   	/**<Number of points for NCPA calibration*/
	char **surf;    /**<OPD surfaces*/
	int nsurf;      /**<Number of OPD surfaces*/
	char **tsurf;   /**<Tilted surfaces, surface, not OPD*/
	int ntsurf;     /**<Number of tilted surfaces*/
}ncpa_cfg_t;
/**
   Parameters for Cn square estimation.
*/
typedef struct cn2est_cfg_t{
    dmat *pair;      /**<If non empty, paris of WFS to use for cn2
			 estimation. Empty: disable cn2 estimation*/
    int step;        /**<do cn2 estimation every this time step*/
    int reset;       /**<reset the accumulated cn2 after every cn2step.*/
    int tomo;        /**<update tomography parameters if non zero*/
    int verbose;     /**<1:Print out estimated r0, cn2 during simulation.*/
    int keepht;      /**<>0: use the layer ht specified by atmr.ht. 2: also do slodar
			directly on these layers.*/
    int nhtomo;      /**<number of layers to feed into reconstructor. only
			effective if keepht=0*/
    int moveht;      /**<1: move the ht used for reconstructor to near strongest
			layers. only effective if keepht=0.*/
    int psol;        /**<Use pseudo open loop gradients. 0 to probe residual*/
    real hmax;     /**<maximum height to estimat*/
    real saat;     /**<subaperture area threashold to use in cn2 estimation*/
}cn2est_cfg_t;
/**
   contains input parameters for plotting during simulation. For debug purpose
*/
typedef struct plot_cfg_t{
    int setup;       /**<Plot various information in setup process*/
    int atm;         /**<Plot the generated atmosphere*/
    int run;         /**<Plot information during simulation*/
    int opdx;        /**<Plot turbulence projected onto xloc.*/
    int psf;         /**<Plot PSF in log (1) or linear (2) mode*/
    int grad2opd;    /**<Plot gradients as reconstructed opd*/
    int all;         /**<Enables setup, atm, run*/
	real opdmax;     /**<Set zlim for OPD drawing*/
	real gmax;       /**<Set zlim for gradient drawing*/
    real psfmin;     /**<Set zlim for psf drawing*/
}plot_cfg_t;
/**
   contains input parameters for debugging.
*/
typedef struct dbg_cfg_t{
    int mvstlimit;   /**<Limit number of modes controled on MVST*/
    int annular_W;   /**<Define the W0/W1 on annular aperture instead of circular*/
    lmat *tomo_maxit; /**<if not empty, will study these maxits in open loop*/
    int tomo_hxw;    /**<1: Force use hxw always instead of ray tracing from xloc to ploc.*/
    int ecovxx;      /**<save the xx used to calculate ecov in psfr.*/
    int useopdr;     /**<use opdr in psf reconstruction*/
    int force;       /**<Force run even if Res_${seed}.done exists*/
    int cmpgpu;      /**<1: cpu code follows GPU implementation.*/
    int pupmask;     /**<Testing pupil mask for NGS WFS to be within LGS volume.*/
    int wfslinearity;/**<Study the linearity of this wfs*/
    int nocgwarm;    /**<Disable warm restart in CG*/
    int test;        /**<Temporary any testing purpose*/
    int dmfullfov;   /**<let DM cover full FoV (sim.fov)*/
    int tomo;        /**<Comparing tomography in GPU and CPU*/
    int fit;         /**<Comparing DM fitting in GPU and CPU*/

    int gp_noamp;    /**<Use annular instead of aper.amp for GP*/
    dmat *atm;       /**<test special atmosphere. <0: fourier mode with spatial frequency 1/dbg.atm m^-1. >0: zernike mode*/
    real gradoff_scale;/**<Scale the reference vector*/
	int gradoff_reset;/**<reset gradoff after creating matched filter with dithering*/

    dcell *dmoff;    /**<DM offset for simulating turbulence on the DM. dimension: ndm*nstep*/
    dcell *gradoff;  /**<Introduced additional gradient offset. dimension: nwfs*nstep*/
    int twfsflag;    /**<use TWFS to control 0: all modes, 1: radial only*/
	int twfsrmax;    /**<TWFS maximum zernike radial order.*/

    int wfs_iac;     /**<Cubic spline coupling factor for turbulence fitting onto wfs grid.*/
    int fullatm;     /**<Always copy full atm to GPU.*/
    int lo_blend;    /**<Low order multi-rate control blending scheme.*/
    real eploscale;/**<Scale of eplo*/
    int recon_stuck; /**<Whether to handle stuck actuator in reconstruction.*/
}dbg_cfg_t;
/**
   Configure GPU usage for different parts.
*/
typedef struct gpu_cfg_t{
    int wfs;         /**<Use GPU for wavefront sensor*/
    int evl;         /**<Use GPU for performance evaluation*/
    int tomo;        /**<Use GPU for tomography*/
    int fit;         /**<Use GPU for DM fitting. Options: 1) assembled matrix 2) matrix free for RHS.*/
    int lsr;         /**<Use GPU for least square reconstruction*/
    int recon;       /**<Use GPU for reconstructor any of(tomo, fit, lsr)*/
    int psf;         /**<Use GPU for accumulating PSF. */
    int moao;        /**<Use GPU for moao.*/
}gpu_cfg_t;
/**
   contains input parameters for each MOAO type.
*/
typedef struct moao_cfg_t{
    real dx;       /**<Spacing of MOAO DM act*/
    int order;       /**<Nominal order of this MOAO*/
    int used;        /**<This moao is used*/
    int actslave;    /**<Do we do actuator slaving*/
    int lrt_ptt;     /**<Piston/tip/tilt constraint*/
    real iac;      /**<If !=0: use cubic influence function with this Inter-Actuator Coupling coefficient.*/
    real stroke;   /**<Stroke of the MOAO DM*/
    real gdm;      /**<The gain of type I controller. a[n]=a[n-1]+e*g where g=o[n]-a[n-1]*/
    real ar;       /**<Aspect ratio dy/dx*/
    real guard;
    dmat *actfloat;  /**<file containing floating actuators. nx2 coordinate*/
    dmat *actstuck;  /**<file containing stuck actuators. nx2 coordinate.*/
}moao_cfg_t;
/**
   contains input parameters for reusing of saved variables.
*/
typedef struct load_cfg_t{
    char *atm;       /**<load atmosphere from. Contains cell array of square matrix*/
    char *locs;      /**<load aper_locs from*/
    char *aloc;      /**<load DM aloc from*/
    char *xloc;      /**<load xloc for recon from*/
    char *ploc;      /**<load ploc for recon from*/
    char *floc;      /**<load floc for recon from*/
    char *cxx;       /**<load laplacian from to do Cxx^-1 in tomo.*/
    //char *HXF;       /**<load HXF from.*/
    char *HXW;       /**<load HXW from.*/
    //char *HA;        /**<load HA from.*/
    char *GP;        /**<load GP from.*/
    char *GA;        /**<load GA from.*/
    char *mvm;       /**<load mvm from.*/
    char *mvmi;      /**<load mvmi from.*/
    char *mvmf;      /**<load mvmf from.*/
    int mvst;        /**<load MVST mvst_U and mvst_FU. see recon.c*/
    int GS0;         /**<if 1, load GS0 from powfs%d_GS0.bin*/
    int tomo;        /**<if 1, load tomo matrix*/
    int fit;         /**<if 1, load fit matrix*/
    int W;           /**<if 1, load W0, W1*/
    char *ncpa;      /**<Load ncpa from this path. saved by save.ncpa*/
	char *saneai;    /**<load subaperture NEA inverse from file for wavefront reconstruction */
}load_cfg_t;
/**
   contains input parameters for saving variables.
*/
typedef struct save_cfg_t{
    int extra;       /**<Save extra results, namely clep, olep, cleNGSmp, etc*/

    int all;         /**<save absolutely everything. mainly for debugging*/
    int setup;       /**<save preparation matrices*/
    int recon;       /**<save reconstructor information. large*/
    int mvst;        /**<MVST computation intermediate matrices*/
    int ncpa;        /**<save NCPA surface OPD on aper and powfs*/
    int fdpcg;       /**<save FDPCG matrices*/
    /*run time */
    int atm;         /**<save atmosphere*/
    int run;         /**<save run time informaton for each time step*/
    int opdr;        /**<save reconstructed OPD on XLOC for each time step*/
    int opdx;        /**<save ATM propagated to XLOC for each time step*/
    int dm;          /**<save computed DM actuator commands for each time step*/
    int evlopd;      /**<save science OPD every `evlopd` time steps*/
    int dither;      /**<save estimated matched filter from dithering*/
    int gradoff;     /**<save gradient reference vector*/
    /*for WFS. Converted from scalar or vector input.
      Scalar input: 1: both, 2: high order only, 3: lo only
      Vector input: Equal to the number of WFS.
    */
    lmat *wfsopd;    /**<save WFS OPD:*/
    lmat *ints;      /**<save WFS subaperture image*/
    lmat *grad;      /**<save WFS gradients*/
    lmat *gradnf;    /**<save WFS noise free gradients*/
    lmat *gradpsol;  /**<save WFS PSOL gradients*/
    lmat *gradgeom;  /**<save WFS geometric gradient during physical optics simu*/
    /*The following are derived from above. */
    int wfsopdhi;    /**<save high order WFS OPD(derived)*/
    int wfsopdlo;    /**<save low order WFS OPD(derived)*/
    int intshi;      /**<save high orrder WFS integration(derived)*/
    int intslo;      /**<save low orrder WFS integration(derived)*/
    int gradhi;      /**<save WFS gradients for high order wfs (derived)*/
    int gradlo;      /**<save WFS gradients for low order wfs (derived)*/
    int gradgeomhi;  /**<save WFS geometric gradient during physical optics simulations.(derived)*/
    int gradgeomlo;  /**<save WFS geometric gradient during physical optics simulations.(derived)*/

    /*Others*/
    int gcovp;       /**<output cumulative gradient covariance average every gcovp step*/
    int ngcov;       /**<number of pairs of gradient covariance to compute*/
    lmat *gcov;      /**<size of 2*ngcov, specifying wfs for each pair*/
    int ecov;        /**<save covariance of DM error vector*/
    int mvmi;        /**<save TomoL output of mvm control matrix assembly for warm restart.*/
    int mvmf;        /**<save FitR output  of mvm control matrix assembly*/
    int mvm;         /**<save computed mvm control matrix*/
}save_cfg_t;
typedef struct dist_cfg_t{
    char **tel2wfs;  /**<Distortion from telescope pupil to each WFS*/
    char **dm2wfs;   /**<Distortion from DM to each WFS. Displacement due to altitude should not be included here*/
    char **dm2sci;   /**<Distortion from DM to science. Not specified for individual science*/
}dist_cfg_t;
/**
   is a wrapper of all _CFG_T data types.
*/
typedef struct parms_t{
    atm_cfg_t    atm;   /**<atmospheric parameters*/
    atmr_cfg_t   atmr;  /**<information about reconstructed atm*/
    aper_cfg_t   aper;  /**<aperture parameters*/
    tomo_cfg_t   tomo;  /**<tomography parameters*/
    fit_cfg_t    fit;   /**<DM fit parameters*/
    lsr_cfg_t    lsr;   /**<LSR parameters*/
    recon_cfg_t  recon; /**<general reconstruction parameters*/
    evl_cfg_t    evl;   /**<Performance evaluation parameters*/

    /*the following are pointers because there may be several*/
    powfs_cfg_t *powfs; /**<Array of wfs type*/
    wfs_cfg_t   *wfs;   /**<Array of wfs*/
    wfsr_cfg_t   *wfsr;  /**<Array of wfs used in reconstruction. Has only 1 wfs
			   per powfs in glao mode, otherwise same as wfs.*/
    //pywfs_cfg_t *pywfs; /**<Array of Pyramid WFS pywfs*/
    dm_cfg_t    *dm;    /**<Array of DM*/
    moao_cfg_t  *moao;  /**<Array of MOAO*/

    sim_cfg_t    sim;   /**<Simulation information*/
	ncpa_cfg_t	 ncpa;  /**<Surface and NCPA caibration parameters.*/
    cn2est_cfg_t cn2;   /**<Parameters for Cn2 estimation*/
    plot_cfg_t   plot;  /**<Specify what to plot during simulation.*/
    dbg_cfg_t    dbg;   /**<Specify debugging parameters*/
    gpu_cfg_t    gpu;   /**<Specify GPU options.*/
    load_cfg_t   load;  /**<Specify what matrices to load for debugging*/
    save_cfg_t   save;  /**<Specify what to save to file for debugging*/
    dist_cfg_t   distortion;/**<Field distortion. Not misregistration*/
    int npowfs;      /**<Number of wfs types*/
    int nwfs;        /**<Number of wfs*/
    int nwfsr;       /**<Number of wfs used in reconstruction. =npowfs in glao, =nwfs otherwise*/
    //int npywfs;      /**<Number of Pyramid WFS*/
    int ndm;         /**<Number of DMs*/
    int nmoao;       /**<Number of different MOAO type*/

    int *fdlock;    /**<Records the fd of the seed lock file. if -1 will skip the seed*/
	char **fnlock;	 /**<Records the filename of the seed lock file. */
    int nlopowfs;    /**<Number of low order wfs types*/
    lmat *lopowfs;   /**<List of low order powfs*/
    int nhipowfs;    /**<Number of high order wfs types*/
    lmat *hipowfs;   /**<List of high order powfs*/
    int ntrpowfs;    /**<Number of tip/tilt removed wfs types*/
    int ntipowfs;    /**<Number of tip/tilt include wfs types*/
    int nphypowfs;   /**<Number of powfs with local/uplink tip/tilt loop*/
    int nlowfs;      /**<Number of low order wfs.*/
    int nhiwfs;      /**<Number of high order wfs*/
    dmat *dirs;      /**<Collect for beam directions*/
    int dither;      /**<Some WFS is doing dithering*/
    int ilgspowfs;   /**<Index of LGS WFS*/
    int nlgspowfs;   /**<Number of LGS WFS*/
	int ittfpowfs;   /**<Index of ttf lo powfs. -1 if none*/
	int ittpowfs;   /**<Index of tt lo powfs. -1 if none*/
    int itpowfs;     /**<Index of twfs*/
    int idmground;   /**<Index of ground dm. default to 0*/
    int step_lo;     /**<Enabling step for low order wfs*/
    int step_hi;     /**<Enabling step for high order wfs*/
    real hipowfs_hsmin; /**<high order wfs minimum height*/
	real hipowfs_hsmax; /**<high order wfs maximum height*/
    int itwfssph;    /**<index of TWFS spherical mode*/
}parms_t;
/**
   arg_t is used for command line parsing.
*/
typedef struct arg_t{
    int detach;      /**<Detach from the command line and run in background*/
    int over;        /**<Run simulation even if Res_${seed}.done exists*/
    int force;       /**<For start, bypassing scheduler*/
    int *gpus;       /**<Index of GPU to use. -1 to disable*/
    int ngpu;        /**<Number of entries in gpus*/
    int ngpu2;       /**<Number of GPUs to use. Ignore of gpus is set.*/
    int server;      /**<MAOS acting as server*/
    char *dirout;    /**<Result output directory*/
    char *confmain;  /**<master .conf file. nfiraos.conf by default. -c to change*/
    char *confcmd;   /**<Additional configuration options supplied in command line.*/
    char *host;      /**<Run in another host*/
    char *execmd;    /**<concatenation of argv*/
}arg_t;
parms_t* setup_parms(const char *mainconf, const char *extracmd, int override);
void setup_parms_gpu(parms_t *parms, int *gpus, int ngpu);
void free_parms(parms_t *parms);
/*The following are here so that we don't have to include type.h or utils.h */
/*convenient constants. used in utils.c */
typedef enum T_TYPE{
    T_PLOC=0,
    T_ALOC,
    T_XLOC,
    T_ATM,
}T_TYPE;

enum WFS_TYPE{ //WFS type
    WFS_SH=0,//shack-hartmann WFS
    WFS_PY=1, //pyramid WFS
    WFS_TOT
};
enum GTYPE{//geometric WFS gradient algoirthm
    GTYPE_G=0,//G tilt: average gradient over aperture (high order WFS)
    GTYPE_Z=1,//Z tilt: zernike best fit tilt (low order WFS)
    GTYPE_TOT
};
enum PTYPE{//physical optics pixel processing algorithm
    PTYPE_MF=1, //matched filter
    PTYPE_COG=2,//COG
    PTYPE_MAP=3,//maximum apriori (deprecated)
    PTYPE_CORR=4,//Correlation peak first
    PTYPE_CORRS=5, //Correlation sum first (deprecated)
    PTYPE_TOT
};

enum NCPA_METHOD{//method to handle NCPA
    NCPA_G=1,  //with gradient offset
    NCPA_I0=2, //with i0 offset (opd)
    NCPA_TOT
};

enum RECON_ALG{
    RECON_MVR=0,//least squares reconstructor
    RECON_LSR=1,//minimum reconstructor
    RECON_TOT
};
enum ALG_TOMO_FIT{
    ALG_CBS=0, //cholesky backsolve
    ALG_CG=1,  //conjugate gradients
    ALG_SVD=2, //svd inverse
    ALG_BGS=3, //block gauss seidel
    ALG_TOT
};

enum PCG_ALG{
    PCG_NONE=0,  //no PCG
    PCG_FD=1,    //fourier domain PCG
    PCG_TOT
};

#endif

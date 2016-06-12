/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#define MAX_AMPS 5
extern double TOMOSCALE;
/*********************************************/
/*             Structs for parms             */
/*********************************************/
/**
   \file maos/parms.h
   Store the input parameters.
*/
/**
   contains input parameters for the atmospheric turbulence.
*/
typedef struct ATM_CFG_T{
    double r0z;   /**<r0 at zenith*/
    double r0;    /**<derived from r0z for zenith angle za*/
    double L0;    /**<outer scale*/
    double dx;    /**<sampling of turbulence screens*/
    double hmax;  /**<maximum in ht*/
    dmat *r0logpsds; /**<[alpha beta]: temporal PSD of log(r0) is beta*f^alpha. f is in hz.*/
    dmat *r0logpsdt; /**<[alpha beta]: spatial  PSD of log(r0) is beta*f^alpha. f is in m.*/
    dmat *ht;   /**<height of each layer*/
    dmat *wt;   /**<weight of each layer (relative strength of \f$C_n^2\f$)*/
    dmat *ws;   /**<wind speed of each layer*/
    dmat *wddeg;/**<wind direction of each layer*/
    dmat *size; /**<size of atm in meter, [0 0]: automatic*/
    lmat *overx;  /**<maximum pixel distance in x direction the beam can be without wrapping*/
    lmat *overy;  /**<maximum pixel distance in y direction the beam can be without wrapping*/
    mapcell *(*fun)(GENATM_T*); /**<Points to the function used to generated atmosphere. (derived)*/
    int nps;      /**<number of phase screens*/
    int wdrand;   /**<randomize wind direction*/
    int iground;  /**<index into the ground layer*/
    lmat *ipsr;    /**<corresponding reconstruction layer*/
    int nx;       /**<turbulence screen size along x*/
    int ny;       /**<turbulence screen size along y*/
    int nxm;      /**<minimum turbulence screen size along x for fft*/
    int nym;      /**<minimum turbulence screen size along y for fft*/
    int nxn;      /**<minimum turbulence screen size along to cover meta pupil*/
    int nyn;      /**<minimum turbulence screen size along to cover meta pupil*/
    int fractal;  /**<1: Use fractal method to generate atmosphere screen.*/
    int frozenflow;  /**<frozen flow. automatic if closeloop=1*/
    int ninit;    /**<Initial size of the screen in fractal method. >=2*/
    int share;    /**<0: disable sharing of atmosphere using file backend*/
    int r0evolve; /**<Evolve r0 according to r0logpsd*/
}ATM_CFG_T;
/**
   contains input parameters for the atmospheric reconstruction.  */
typedef struct ATMR_CFG_T{
    double r0z;   /**<r0 at zenith*/
    double r0;    /**<derived from r0z for zenith angle za*/
    double L0;    /**<outer scale*/
    double hs;    /**<height of the high order guide star. derived*/
    double hmax;  /**<maximum of ht*/
    dmat *ht;   /**<height of each layer*/
    dmat *wt;   /**<weight of each layer (relative strength of \f$C_n^2\f$)*/
    double dx;    /**<baseline sampling (when os=1). matches to high order wfs.*/
    lmat *indps;   /**<Mapping atmr.ps to atm.ps*/
    lmat *os;      /**<over sampling factor of xloc over actuator spacing */
    int nps;      /**<number of phase screens*/
}ATMR_CFG_T;
/**
   contains input parameters about the aperture, like the diameter,
   amplitude map, etc */
typedef struct APER_CFG_T{
    double d;     /**<Telescope aperture diameter*/
    double din;   /**<Telescope inner blocking diameter*/
    double rotdeg;/**<pupil rotation in degree*/
    char *fnamp;  /**amplitude maps. expected to be square or rectangular mxn, with 0 at
		     [m/2,n/2] (count from 0)*/
    int fnampuser;/**<User provided amplitude map (not default)*/
    char *pupmask;/**<The pupil cold stop*/
}APER_CFG_T;
/**
   contains input parameters for laser launch telescope
*/
typedef struct LLT_CFG_T{
    double d;      /**<LLT clear aperture diameter*/
    double widthp; /**<Gaussian beam width percentage of d*/
    char *ttpsd;   /**<PSD of uplink beam jitter*/
    char *fnrange; /**<File contains range to sodium layer*/
    char *fnprof;  /**<File contains sodium profile*/
    char *fnsurf;  /**<Pupil Surface OPD error*/
    char *fnamp;   /**<Pupil amplitude map. overrides widthp*/
    dmat *ox;    /**<location x of LLT center wrt telescope aperture center*/
    dmat *oy;    /**<see ox.*/
    dmat *misreg;
    lmat *i;        /**<Index into llt for this iwfs.*/
    int ttfr;      /**<Remove piston/tip/tilt/focus from ncpa*/
    int n;         /**<number of launch telescopes in this powfs*/
    int colprep;   /**<starting column to use in fn for ETF in preparation of
		      matched filter*/
    int colsim;    /**<starting column to use in fn for ETF in simulation*/
    int colsimdtrat;/**<change to next sodium profile during simulation every
		       colsimdtrat time step*/
    double ttrat;  /**<Ratio of uplink jitter to downlink jitter due to M2 windshake.*/
} LLT_CFG_T;
/**
   contains input parameters for each type of wfs (powfs).
*/
typedef struct POWFS_CFG_T{
    dmat *wvl;   /**<list of wavelength*/
    dmat *wvlwts; /**<weights for each wavelength. can be overriden by wfs.wvlwts.*/
    char *saloc;   /**<saloc override file*/
    char *piinfile;/**<input averaged pixel intensities for matched filter. NULL
		      to disable*/
    char *sninfile;/**<Speckle noisy input file. NULL to disable. not used*/
    double hs;     /**<height of guide star*/
    double saat;   /**<subaperture area (normalized) threshold to drop subaperture.*/
    double safill2d;/**<subaperture lenslet throughgput. value is used  to alter amplitude map*/
    double saspherical;/**<Subaperture spherical aberration in nm RMS at best focus.*/
    double safocuspv;   /**<Subaperture focus error in nm p/v*/
    char  *neareconfile;/**<prefix of file contains noise equivalent angle in
			   radian. _wfs# is added when reading file.*/
    double nearecon;/**<NEA used in reconstruction*/
    double neasim;  /**<NEA used in simulation. -1 to use nearecon*/
    char*  neasimfile;/**<read NEA used in simulation from file. Defined at
			 sim.dt sampling rate, in radian. neasim must be -1*/
    double neaspeckle;/**<NEA caused by speckle noise. Added to matched filter
			 estimation of NEA due to photon and detector noise in
			 physical optics mode for reconstructor*/
    double bkgrnd;  /**<background in electron per pixel per LGS frame*/
    double bkgrndc;/**<How much of the background in bkgrnd can be calibrated
		      out. depends on variability.*/
    char *bkgrndfn; /**<file contains sky background/rayleigh scatter input for
		       each subaperture in each wfs. */
    char *bkgrndfnc;/**<How much of the background in bkgrndfn can be
		       calibrated out. depends on variability.*/
    double rne;     /**<read out noise in electron per pixel per frame*/
    double pixblur; /**<pixel bluring due to leakage. relative to pixel size.*/
    double dsa;     /**<Size of subaperture in 1 dimension*/
    double dx;      /**<sampling of opd points in each subaperture. usually
		       matches atmosphere sampling for LGS. may be coraser for NGS.*/
    double pixtheta;/**<size of pixel pitch along x/y or azimuthal if radial
		       ccd. Converted to radian from user input*/
    double radpixtheta; /**<size of pixel pitch along radial direction. -1 for square pixel*/
    double radgx;   /*Create and/or use gx, gy along radial/azimuthal direction.*/
    double fieldstop;/**<size of field stop in arcsec.*/
    double pixoffx; /**<offset of image center from center of detector*/
    double pixoffy; /**<see pixoffx*/
    double sigscale;/**<scale the signal level for simulation.*/
    double siglev;  /**<signal level. will be override by wfs.siglev is specified.*/
    struct LLT_CFG_T *llt;/**<configuration for LLT*/
    char* fnllt;    /**<filename of LLT configuration. empty means no llt.*/
    int type;       /**<WFS type: 0: SHWFS, 1:Pyramid WFS*/
    int step;       /**<frame to start using WFS*/
    int trs;        /**<tip/tilt removal flag. True for LGS, False for NGS*/
    int dfrs;       /**<differential focus removal flag. True for LGS, False for NGS*/
    int lo;         /**<whether this is a low order wfs. False for LGS, True for NGS*/
    int skip;       /**<skip in high order tomography, for split tomo (derived parameter)*/
    int psol;       /**<Compute pseudo open loop gradients (derived parameter)*/
    lmat *wfs;       /**<array of wfs belongs to this powfs*/
    lmat *wfsind;    /**<wfsind[iwfs] gives the index of the wfs in this powfs group*/
    int nwfs;       /**<number of wfs belonging to this powfs*/
    int nwfsr;      /**<number of wfs for reconstruction belonging to this powfs*/
    int neaphy;     /**<use nea from physical optical precomputation in geometric simulations.*/
    int phyusenea;  /**<force using supplied noise equivalent angle in physical
		       optics simulations*/
    int order;      /**<order of wavefront sensing along one dimension.*/
    int pixpsa;     /**<number of detector pixels along x/y or azimuthal if radial CCD.*/
    int radpix;     /**<number of detector pixels along radial direction if radial CCD*/
    int radrot;     /**<For radial format CCD, rotate OTF into coordinate plane. uses less memory*/
    int ncomp;      /**<number of PSF points before forming detector image. 0 for automatic*/
    int embfac;     /**<Embed subaperture atm OPD before fft. set to 2.*/
    int nwvl;       /**<Number of wavelength. 1 for LGS, 2 for NGS J+H sensing.*/
    int gtype_sim;  /**<wfs type if not using physical optics in reconstruction. 
		       - 0: geometric
		       - 1: ztilt.*/
    int gtype_recon;/**<wfs type if not using physical optics in simulation. 
		       - 0: geometric
		       - 1: ztilt.*/
    int phytype;    /**<physical optics type for reconstruction. 1: mtch, 2: tcog, 3: MAP*/
    int phytypesim; /**<physical optics type for simulation. -1 to follow phytype*/
    int phytypesim2;/**<physical optics type after dithering update. -1 to follow phytypesim*/
    int phystep;    /**<frame to start using physical optics. 
		       -  0: means from frame 0.
		       - >0: need to compute GS0 to calculate geometric optics
		       - -1: never, doesn't need to compute DTF
		    */
    int usephy;     /**<whether physical optics is used at all during
		       simulation.(derived parameter)*/
    double r0;  /**<Fried parameter  for matched filter generation. Uses atm.r0, atm.L0 is not set*/
    double L0;  /**<Outerscale for matched filter generation. Uses atm.r0, atm.L0 is not set*/
    double mtchcr;  /**<if >0 use constrained mtch for this amount of pixels*/
    double mtchcra; /**<if >0 use constrained mtch for azimuthal for this amount of pixels*/
    int mtchcpl;    /**<use coupling between r/a measure error. useful for LGS with x-y ccd.*/
    int mtchstc;    /**<shift peak in the time averaged short exposure PSF to center using fft.*/
    int mtchscl;    /**<scale subaperture image to have the same intensity as i0. Keep false.*/
    int mtchadp;    /**<Using adaptive matched filter. When the number of pixels
		       in the image brighter than half maximum is more than this
		       value, use constraint. introduced on 2011-02-21.*/
    double cogthres;/**<CoG threshold, relative to max(im)*/
    double cogoff;  /**<CoG offset to remove, relative to max(im). */
    int needGS0;    /**<need to compute GS0 (derived parameter)*/
    int noisy;      /**<noisy or not during *simulation* */
    int ncpa_method;/**<Method to correct ncpa.
		       - 0: do nothing.
		       - 1: apply gradient electronic offset. 
		       - 2: apply ncpa to average pixel intensity i0, better than 1 */
    int pistatout;  /**<output time averaged short exposure image. 1: shift to center, 2: do not shift*/
    int pistatstart;/**<time step to compute pistatout*/
    int pistatstc;  /**<1: shift to center using fft method. 0: use geometric gradients.*/
    int psfout;     /**<output time history of low order wfs PSF. never do this for LGS.*/
    int dtrat;      /**<ratio of sample period over fast loop (LGS)*/
    int idtrat;     /**<Index of dtrat into parms->sim.dtrats*/
    int i0scale;    /**<scale i0 to matched subaperture area.*/
    int moao;       /**<index into MOAO struct. -1: no moao*/
    int dither;     /**<Turn on/off dithering to update centroid gain or matched filter*/
    double gradscale;/**<Scale CL gradients. For testing*/
    double dither_amp; /**<Dither amplitude in arcsec for tip/tilt mode*/
    double dither_gpll;/**<Gain of phase locked loop*/
    double dither_gog; /**<Gain for updating optical gain for cog*/
    int dither_npoint;/**<Number of points in each dither peroid (4)*/
    int dither_pllskip;/**<Skip WFS frames for uplink loop to stable*/
    int dither_pllrat; /**<Number of WFS frames for updating PLL.*/
    int dither_ogskip; /**<Number of WFS frames to skip before computing averaged images*/
    int dither_ograt;  /**<Number of WFS frames to update pixel processing algorithm (MF/CoG)*/
    //options for zoom corrector
    int zoomdtrat;   /**<dtrat of the trombone averager*/
    int zoomshare;   /**<1: All LGS share the same trombone*/
    double zoomgain; /**<gain of the trombone controller*/
    int zoomset;     /**<Set zoom position from the beginning*/
    /*Options for Pywfs*/
    double modulate;  /**<Pyramid modulation diamter in arcsec*/
    int    modulpos;  /**<Number of positions per modulation cycle*/
}POWFS_CFG_T;
/**
   contains input parmaeters for each wfs
*/
typedef struct WFS_CFG_T{
    dmat *wvlwts; /**<Weights of signal value for each wavelength. if not
		       specified in config, will use powfs.wvlwts*/
    char *sabad;    /**<coordinate of bad subaperture due to bad detector or lenslet array.*/
    double thetax;  /**<x direction*/
    double thetay;  /**<y direction*/
    double siglev;  /**<Total signal value for all wavelength. if not specified
		       in config, will use powfs.siglev*/
    double siglevsim;/**<Signal value used for simulation. (derived parameter)*/
    double hs;      /**height of star. Derived from powfs.hs or from input*/
    double fitwt;   /**<Include wfs in fitting directions if corresponding wfs[iwfs].fitwt is greater than 0*/
    int powfs;      /**<powfs type*/
}WFS_CFG_T;
/**
   contains input parameters for each deformable mirror.
*/
typedef struct DM_CFG_T{
    double guard;   /**<extra DM actuator rings outside of aper.d*/
    dmat *stroke;  /**<Stroke of DM (surface). OPD goes to \f$\pm\f$ stroke. array: per actuator$*/
    double iastroke;/**<Inter actuator stroke (surface)*/
    char *iastrokefn;   /**< describes polynomials that convert
			  * opd to voltage (first cell), and voltage to opd
			  * (second cell). The two operations has to be strict
			  * inverse of each other*/
    dcell *iastrokescale;/**<Input from iastrokefn*/
    double vmisreg; /**<vertical misregistration*/
    double ht;      /**<height conjugation range*/
    double dx;      /**<actuator separation along x (derived from order)*/
    double ar;      /**<[in] aspect ratio: dy/dx*/
    double dy;      /**<actuator separation along y (derived from dx and ar*/
    double offset;  /**<Center-most actuator offset from origin
		       - =0 means there is a act on center. 
		       - 1/2 means no one in the center.*/
    double iac;     /**<If !=0: use cubic influence function with this Inter-Actuator Coupling coefficient.*/
    double histbin; /**<The bin width for histogram.*/
    int histn;      /**<Number of bins in histogram.*/
    int hist;       /**<Compute histogram of commands of each actuator*/
    double order;   /**<Order of the DM within telescope clear subaperture*/
    int isground;   /**<Is this DM the ground DM (derived)*/
    char *actfloat; /**<file containing floating actuators. nx2 coordinate*/
    char *actstuck; /**<file containing stuck actuators. nx2 coordinate.*/

    char *hyst;     /**<File containing a matrix that describes the
		       hysterisis. The matrix should have 3 rows, and multiple
		       columns. Each column describe the parameters for a mode,
		       and the rows are: 1) weight 2) alpha, 3) beta. See
		       "Hysteresis Modeling and Simulation" by Luc Gilles*/
}DM_CFG_T;
/**
   contarins input parameters all evaluation directions.  */
typedef struct EVL_CFG_T{
    dmat *thetax; /**<x Coordinate of evaluation directions*/
    dmat *thetay; /**<y Coordinate of evaluation directions*/
    dmat *wt;     /**<weight of each direction*/
    dmat *wvl;    /**<wavelength for PSF and strehl computation*/
    dmat *hs;     /**<height of each science object*/
    double dx;     /**<sampling of aperture for evaluation*/
    int nwvl;       /**<Number of wavelength*/
    lmat *psf;       /**<1: participate in psf evaluation.*/
    lmat *psfr;      /**<1: participate in psf reconstruction telemetry*/
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
    int cov;        /**<save covairance of science OPD ,every this time step,
		       for directions where evl.psf is 1*/
    lmat *pttr;      /**<remove p/t/t from psf. 1 number for each evl.*/
    lmat *psfngsr;   /**<remove ngs modes from psf.*/


    int psfisim;    /**<time step to start psfmean.*/
    lmat *psfsize;    /**<save this number of pixels of the center of the psf. 1
			number for each wvl.*/
    lmat *psfgridsize;/**<grid size for FFT to generate PSF. Becareful about FFT
			speed and enough padding. Determines the sampling of the
			generated PSF. 0 or negative for automatic. 1 number for
			each wvl.*/
    int nevl;       /**<Number of evaluation directions. (derived)*/
    int tomo;       /**<evaluate tomography performance.*/
    int indoa;      /**<index of the on axis evluation point.*/
    int moao;       /**<index into MOAO struct. -1: no MOAO*/
}EVL_CFG_T;

/**
   contains input parameters for wavefront tomography.
*/
typedef struct TOMO_CFG_T{
    double tikcr;    /**<tikhonov regularization.*/
    double minwt;    /**<minimum layer weight allowed. if less than this will force to this.*/
    double iac;      /**<!=0: use cubic influence function with this Inter-Actuator Coupling coefficient.*/
    double cxxscale; /**<scale the Cxx^-1 term.*/
    double svdthres; /**<Threshold in SVD inversion*/
    double cgthres;  /**<Repeat cg if residual is not reached*/
    int square;      /**<use square/rectangular grid instead of tighter irregular grid*/
    int cone;        /**<use cone coordinate in xloc: keep true*/
    int cxx;         /**<method to compute Cxx^-1. 0: bihormonic approx. 1: inverse psd. 2: fractal*/
    int guard;       /**<guard rings of reconstruction grid xloc*/
    int pos;         /**<over sampling factor of ploc over actuator spacing*/
    int nxbase;      /**<Each layer xloc grid size is tomo.os*tomo.nxbase is not zero. same for ploc.*/
    int piston_cr;   /**<single point piston constraint. */
 
    int ahst_wt;     /**<0: use Wg, 1: using Wa*/
    int ahst_idealngs;/**<ideal correction on NGS modes. For skycoverage preprocessing.*/
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
    int assemble;    /**<force assemble tomography matrix in CG*/
    int predict;     /**<test predictive control.*/
    int ninit;       /**<like atm.ninit, the initial screen to generate from covariance directly*/
    int splitlrt;    /**<1: use LGS low rank terms in split tomography.*/
}TOMO_CFG_T;
/**
   contains input parameters for deformable mirror fitting.
*/
typedef struct FIT_CFG_T{
    dmat *thetax;  /**<x Coordinate of DM fitting directions. */
    dmat *thetay;  /**<y Coordinate of DM fitting directions. */
    dmat *wt;      /**<weight of each direction*/
    dmat *hs;      /**<height of target in each direction*/
    double tikcr;    /**<tikhonov regularization*/
    double svdthres; /**<Threshold in SVD inversion*/
    double actthres; /**<When actuator coupling coefficient drops below this, start slaving.*/
    int actslave;    /**<slaving constraint for non-active actuators. Useful in CBS method*/
    int actinterp;   /**<interpolate actuator results to non-active actuators.*/
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
    int square;      /**<using square grid on DM and ploc.*/
    int assemble;    /**<force assemble fit matrix in CG*/
    int pos;         /**<over sampling of floc over aloc. for fitting. normally equal to tomo.pos*/
    int indoa;       /**<Index of on axis point.*/
    int cachedm;     /**<Cache DM command in intermediate plane*/
    int cachex;      /**<Cache X (xloc) in intermediate plane*/
}FIT_CFG_T;
/**
   contains input parameters for the least square reconstructor.
*/
typedef struct LSR_CFG_T{
    double tikcr;    /**<tikhonov regularization*/
    double svdthres; /**<Threshold in SVD inversion*/
    double actthres; /**<When actuator coupling coefficient drops below this, start slaving.*/
    char  *fnreg;    /**<File containing a regularization term to add to LL.M*/
    int actinterp;   /**<interpolate actuator results to non-active actuators .*/
    int actslave;    /**<slaving constraint for non-active actuators. Useful in CBS method*/
    int bgs;         /**<1: use BGS, block Gaussia Seidel then use alg to solve each block.*/
    int alg;         /**<algorithm to solve the linear equation.
			0: Cholesky direct solve for the large matrix.  (CBS)
			1: CG or PCG.
			2: SVD or EVD: Eigen value decompsition
		     */
    int maxit;       /**<max iterations. Usually 30 for CG*/
}LSR_CFG_T;
/**
   contains input parameters for wavefront reconstruction.
*/
typedef struct RECON_CFG_T{
    int alg;        /**<algorithm for reconstruction. 0: MVR. 1: LSR. moved from sim.recon*/
    int glao;       /**<whether we are in GLAO mode where all WFS in each powfs are averaged*/
    int split;      /**<split reconstruction/tomography type.
		       - 0: integrated tomography
		       - 1: adhoc split tomography
		       - 2: minimum variance split tomography (only valid if recon.alg=0)*/
    int modal;       /**-2: emulate zonal, -1: zernike, 0: zonal, 1: KL modes*/
    int nmod;        /**<Maximum number of modes to control in modal controller*/
    int warm_restart; /**<Warm restart in CG*/
    int psol;        /**<Use pseudo open loop gradients for wavefront reconstruction*/
    int mvm;        /**<Use the various algorithms recon.alg to assemble a control
		       matrix to multiply to gradients to get DM commands. If
		       the algorithm needs PSOL gradient, we will have an
		       auxillary matrix to multiply to the DM actuators and
		       subtract from the result.*/
    char **misreg_dm2wfs; /**<Distortion from DM to each WFS model used in reconstruction. Affects GA*/
    char **misreg_dm2sci; /**<Distortion from DM to each science model used in reconstruction. Affects HA*/
    char **misreg_tel2wfs;/**<Distortion from Telescope to each WFS model used in reconstruction. Affects HXW*/
    double poke;    /**<How much WFE (meter) to apply to OPD for computing experimental interaction matrix*/
    int psd;        /**<Flag: compute PSDs of DM error signal averaged over aperture and field points.*/
    int psddtrat;   /**<how many time step to sample for PSD computation.*/
    int psddtrat_lo;   /**<how many time step to sample for low order PSD computation.*/
    int psddtrat_twfs; /**<how many time step to sample for TWFS PSD computation.*/
    int psdnseg;    /**<#how many overlapping partitions of the time history to compute PSD.*/
    char *fnsphpsd; /**<PSD of spherical aberration due to profile evolution.*/
}RECON_CFG_T;
/**
   contains input parameters for simulation, like loop gain, seeds, etc.
*/
typedef struct SIM_CFG_T{
    double dt;       /**<sampling period. 1/800*/
    double dtlo;       /**<low order wfs sampling period*/
    double dthi;       /**<high order wfs sampling period*/
    double za;       /**<zenith angle in radian*/
    int dtrat_hi;      /**<ratio of sampling period over clock of high order wfs*/
    int dtrat_lo;      /**<dtrat of the lower order loop.*/
    int dtrat_skip;   /**<dtrat (over sim.dt) for frame drop. Be careful when powfs.dtrat is not one.*/
    int start;       /**<time step to start simulation. 0*/
    int end;         /**<time step to stop simulation. exclusive*/
    int pause;       /**<Pause at the end of every time step*/
    lmat *seeds;      /**<simulation seeds*/
    int nseed;       /**<How many simulation seed*/
    int closeloop;   /**<closed loop or open loop*/
    char *wspsd;     /**<Telescope wind shake PSD input. Nx2. First column is
			freq in Hz, Second column is PSD in rad^2/Hz.*/
    int wsseq;       /**<sequence of wind shake time series.*/
    /*control */
    dmat *apdm;      /**<servo coefficient for high order dm.  A is command. e is
			error signal. at time step n, the command is updated by
			A(n)=A(n-1)*apdm(0)+A(n-2)*ap(1)+...+e(n-2)*ep
		     */
    dmat *epdm;      /**<error gain for DM commands (high order)*/
    dmat *aplo;      /**<servo coefficient for ngs modes.*/
    dmat *eplo;      /**<error gain for NGS modes (low order)*/
    dmat *apfsm;     /**<servo coefficient for for LGS uplink pointing loop.*/
    dmat *epfsm;     /**<error gain for uplink pointing*/
    int aldm;        /**<Additional latency (*sim.dt) of the high order loop*/
    int allo;        /**<Additional latnecy (*sim.dt) of the low order loop*/
    int alfsm;       /**<Additional latency (*sim.dt) of the uplink loop*/
    int commonfsm;   /**<Make FSM common for each powfs (LLT)*/
    double zetafsm;  /**<Damping of FSM modeled as second harmonic oscillater (SHO).*/
    double f0fsm;    /**<Resonance frequency of FSM (SHO). 0: infinite.*/
    double aptwfs;   /**<Twfs reference vector servo coefficient.*/
    double eptwfs;   /**<Twfs reference vector servo gain.*/
    double fcttm;    /**<cross over frequency of tip/tilt split*/
    double lpttm;    /**<los path filter for ttm. derived: lpttm=2*pi*fcttm*sim.dt*/
    double fcfocus;  /**<cross-over frequency of the focus LPF.*/
    double lpfocushi;  /**<derived: lpfocus=2*pi*fc*sim.dthi*/
    double lpfocuslo;  /**<derived: lpfocus=2*pi*fc*sim.dtlo*/
    double fov;      /**<The diameter of total fov in arcsec*/
    int focus2tel;   /**<Offload focus to telescope*/
    double epfocus2tel;/*Gain for telescope focus control*/
    int mffocus;     /**<method for focus tracing.
			- 0: no focus tracking.
			- 1: Focus tracking using CL gradients, for each LGS independently.
			- 2: Focus tracking using CL gradinets, for common LGS focus only.
		     */
    int idealfsm;    /**<ideal compensation for uplink pointing*/
    int servotype_hi;/**<servo type for high order loop. 1: simple integrator*/
    int servotype_lo;/**<servo type for low order loop. 1: simple integrator. 2: type II*/
    int cachedm;     /**<cache dm shape on fine sampled grid matched WFS or Science grid*/
    int fuseint;     /**<fuse the high and low order integrators in split tomography */
    int skysim;      /**<1: we are doing skycoverage preprocessing*/
    int evlol;       /**<evaluate open loop error only*/
    int noatm;       /**<disable atmosphere*/
    int idealfit;    /**<do ideal DM fitting from atmosphere directly.*/
    int idealtomo;   /**<ideal tomography without wfs (directly propagate from
			turbulence). conflicts with idealfit. combine with
			evl.tomo to evaluate its performance.*/
    int psfr;        /**<do PSF reconstruction telemetry*/
    int ecnn;        /**<Calculate WF covariance due to WFS noise cov Cnn.*/
    int wfsalias;    /**<Study the wfs aliasing effect by projecting turbulence
			onto the NULL space of DM.*/
    int idealwfs;    /**<Generates ideal WFS by sensing turbulence with DM range.*/
    int idealevl;    /**<Evaluate performance within DM range.*/
    /* A few derived parameters*/
    int dmclip;      /**<derived: Need to clip actuator stroke*/
    int dmclipia;    /**<derived: Need to clip inter-actuator stroke*/
    int dmproj;      /**<derived: Need to projection atmosphere onto DMspace. */
    int ahstfocus;   /**<New new mode split in ahst + focus tracking*/

    int mvmport;     /**<Non zero: specify which port does the MVM server run on and connect to it for MVM reconstruction.*/
    char *mvmhost;   /**<Which host does the MVM server run*/
    int mvmsize;     /**<number of gradients to send each time. 0 is all.*/
    int mvmngpu;     /**<number of GPUs to use in server*/
    int ncpa_calib;  /**<calibrate NCPA. 1: with all DMs. 2: with only ground DM.*/
    int ncpa_ttr;    /**<Remove average t/t from NCPA for WFS. Equivalent as repositioning WFS. default 1.*/
    dmat *ncpa_thetax; /**<Coordinate for NCPA calibration (arcsec)*/
    dmat *ncpa_thetay; /**<Coordinate for NCPA calibration (arcsec)*/
    dmat *ncpa_wt;     /**<Weight for each point.*/
    dmat *ncpa_hs;     /**<Height of star.*/
    int ncpa_ndir;       /**<Number of points for NCPA calibration*/
    char *dmadd;      /**<Containing dm vector to simulate turbulence (added to integrator output). 
			 It should be cell array (time steps) of cell arry (DMs) of vectors. Can be empty*/
}SIM_CFG_T;
/**
   Parameters for Cn square estimation.
*/
typedef struct CN2EST_CFG_T{
    dmat *pair;       /**<If non empty, paris of WFS to use for cn2
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
    double hmax;     /**<maximum height to estimat*/
    double saat;     /**<subaperture area threashold to use in cn2 estimation*/
}CN2EST_CFG_T;
/**
   contains input parameters for plotting during simulation. For debug purpose
*/
typedef struct PLOT_CFG_T{
    int setup;       /**<Plot various information in setup process*/
    int atm;         /**<Plot the generated atmosphere*/
    int run;         /**<Plot information during simulation*/
    int opdx;        /**<Plot turbulence projected onto xloc.*/
    int all;         /**<Enables setup, atm, run*/
}PLOT_CFG_T;
/**
   contains input parameters for debugging.
*/
typedef struct DBG_CFG_T{
    int wamethod;    /**<method to compute wa for ngsmod removal.*/
    int atm;         /**<test special atmosphere*/
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
    int na_smooth;   /**<1: smooth sodium profile to coarser grid before computing etf*/
    int na_interp;   /**<1: Interpolate sodium profile and use FFT to build etf. 0: direct sum, slow*/
    int ncpa_preload;/**<preload integrator with DM sys flat*/
    int ncpa_nouncorr;/**<1: do not include uncorrelatable error in science path.*/
    int i0drift;     /**<Control drift of i0 by driving it toward gradncpa*/
    double gradoff_scale;/**<Scale the reference vector*/
    dmat *pwfs_psx;  /**<pyramid WFS pupil shift along x (in pixel). pupil ordering: -x+y, +x+y, -x-y, +x-y.*/
    dmat *pwfs_psy;  /**<pyramid WFS pupil shift along y (in pixel).*/
    double pwfs_flate;/**<pyramid flat edge angular width */
    double pwfs_flatv;/**<pyramid flat vertex angular width*/
    double pwfs_pupelong;/**<pyramid pupil (detector) elongation ratio (long axis / short axis).*/
}DBG_CFG_T;
/**
   Configure GPU usage for different parts.
*/
typedef struct GPU_CFG_T{
    int wfs;         /**<Use GPU for wavefront sensor*/
    int evl;         /**<Use GPU for performance evaluation*/
    int tomo;        /**<Use GPU for tomography*/
    int fit;         /**<Use GPU for DM fitting*/
    int lsr;         /**<Use GPU for least square reconstruction*/
    int psf;         /**<Use GPU for accumulating PSF. */
    int moao;        /**<Use GPU for moao.*/
}GPU_CFG_T;
/**
   contains input parameters for each MOAO type.
*/
typedef struct MOAO_CFG_T{
    double dx;       /**<Spacing of MOAO DM act*/
    double order;    /**<Order of this MOAO*/
    int used;        /**<This moao is used*/
    int actslave;    /**<Do we do actuator slaving*/
    int lrt_ptt;     /**<Piston/tip/tilt constraint*/
    double iac;      /**<If !=0: use cubic influence function with this Inter-Actuator Coupling coefficient.*/
    double stroke;   /**<Stroke of the MOAO DM*/
    double gdm;      /**<The gain of type I controller. a[n]=a[n-1]+e*g where g=o[n]-a[n-1]*/
    double ar;       /**<Aspect ratio dy/dx*/
    double guard;
    char *actfloat;  /**<file containing floating actuators. nx2 coordinate*/
    char *actstuck;  /**<file containing stuck actuators. nx2 coordinate.*/
}MOAO_CFG_T;
/**
   contains input parameters for reusing of saved variables.
*/
typedef struct LOAD_CFG_T{
    char *atm;       /**<load atmosphere from. Contains cell array of square matrix*/
    char *locs;      /**<load aper_locs from*/
    char *aloc;      /**<load DM aloc from*/
    char *xloc;      /**<load xloc for recon from*/
    char *ploc;      /**<load ploc for recon from*/
    char *floc;      /**<load floc for recon from*/
    char *cxx;       /**<load laplacian from to do Cxx^-1 in tomo.*/
    char *HXF;       /**<load HXF from.*/
    char *HXW;       /**<load HXW from.*/
    char *HA;        /**<load HA from.*/
    char *GP;        /**<load GP from.*/
    char *GA;        /**<load GA from.*/
    char *mvm;       /**<load mvm from.*/
    char *mvmi;      /**<load mvmi from.*/
    char *mvmf;      /**<load mvmf from.*/
    char *i0;        /**<load i0 for powfs from.*/
    int mvst;        /**<load MVST mvst_U and mvst_FU. see recon.c*/
    int GS0;         /**<if 1, load GS0 from powfs%d_GS0.bin*/
    int tomo;        /**<if 1, load tomo matrix*/
    int fit;         /**<if 1, load fit matrix*/
    int W;           /**<if 1, load W0, W1*/
    char *ncpa;      /**<Load ncpa from this path. saved by save.ncpa*/
}LOAD_CFG_T;
/**
   contains input parameters for saving variables.
*/
typedef struct SAVE_CFG_T{
    int extra;       /**<Save extra results, namely clep, olep, cleNGSmp, etc*/
    
    int all;         /**<save absolutely everything. mainly for debugging*/
    int setup;       /**<save preparation matrices*/
    int recon;       /**<save reconstructor information. large*/
    int mvst;        /**<MVST computation intermediate matrices*/
    int ncpa;        /**<Save NCPA surface OPD on aper and powfs*/
    /*run time */
    int atm;         /**<save atmosphere*/
    int run;         /**<save run time informaton for each time step*/
    int opdr;        /**<save reconstructed OPD on XLOC for each time step*/
    int opdx;        /**<save ATM propagated to XLOC for each time step*/
    int dm;          /**<save computed DM actuator commands for each time step*/
    int evlopd;      /**<save science OPD for each time step*/
    int dither;      /**<save estimated matched filter from dithering*/
    /*for WFS. Converted from scalar or vector input.
      Scalar input: 1: both, 2: high order only, 3: lo only
      Vector input: Equal to the number of WFS.
    */
    lmat *wfsopd;      /**<save WFS OPD:*/
    lmat *ints;        /**<save WFS subaperture image*/
    lmat *grad;        /**<save WFS gradients*/
    lmat *gradgeom;    /**<save WFS geometric gradient during physical optics simu*/
    /*The following are derived from above. */
    int wfsopdhi;    /**<save high order WFS OPD(derived)*/
    int wfsopdlo;    /**<save low order WFS OPD(derived)*/
    int intshi;      /**<save high orrder WFS integration(derived)*/
    int intslo;      /**<save low orrder WFS integration(derived)*/
    int gradhi;      /**<save WFS gradients for high order wfs (derived)*/
    int gradlo;      /**<save WFS gradients for low order wfs (derived)*/
    int gradgeomhi;  /**<save WFS geometric gradient during physical optics simulations.(derived)*/
    int gradgeomlo;  /**<save WFS geometric gradient during physical optics simulations.(derived)*/

    int gcovp;       /**<output cumulative gradient covariance average every gcovp step*/
    int ngcov;       /**<number of pairs of gradient covariance to compute*/
    lmat *gcov;       /**<size of 2*ngcov, specifying wfs for each pair*/

    int mvmi;        /**<save TomoL output of mvm control matrix assembly for warm restart.*/
    int mvmf;        /**<save FitR output  of mvm control matrix assembly*/
    int mvm;         /**<save computed mvm control matrix*/
}SAVE_CFG_T;
typedef struct MISREG_CFG_T{
    char **tel2wfs;  /**<Distortion from telescope pupil to each WFS*/
    char **dm2wfs;   /**<Distortion from DM to each WFS. Displacement due to altitude should not be included here*/
    char **dm2sci;   /**<Distortion from DM to science. Not specified for individual science*/
    dmat   *pupil;   /**<Misregistration of the telescope pupil*/
}MISREG_CFG_T;
/**
   is a wrapper of all _CFG_T data types.
*/
typedef struct PARMS_T{
    ATM_CFG_T    atm;   /**<atmospheric parameters*/
    ATMR_CFG_T   atmr;  /**<information about reconstructed atm*/
    APER_CFG_T   aper;  /**<aperture parameters*/
    TOMO_CFG_T   tomo;  /**<tomography parameters*/
    FIT_CFG_T    fit;   /**<DM fit parameters*/
    LSR_CFG_T    lsr;   /**<LSR parameters*/
    RECON_CFG_T  recon; /**<general reconstruction parameters*/
    EVL_CFG_T    evl;   /**<Performance evaluation parameters*/
    
    /*the following are pointers because there may be several*/
    POWFS_CFG_T *powfs; /**<Array of wfs type*/
    WFS_CFG_T   *wfs;   /**<Array of wfs*/
    WFS_CFG_T   *wfsr;  /**<Array of wfs used in reconstruction. Has only 1 wfs
			   per powfs in glao mode, otherwise same as wfs.*/
    DM_CFG_T    *dm;    /**<Array of DM*/
    MOAO_CFG_T  *moao;  /**<Array of MOAO*/

    SIM_CFG_T    sim;   /**<Simulation information*/
    CN2EST_CFG_T cn2;   /**<Parameters for Cn2 estimation*/
    PLOT_CFG_T   plot;  /**<Specify what to plot during simulation.*/
    DBG_CFG_T    dbg;   /**<Specify debugging parameters*/
    GPU_CFG_T    gpu;   /**<Specify GPU options.*/
    LOAD_CFG_T   load;  /**<Specify what matrices to load for debugging*/
    SAVE_CFG_T   save;  /**<Specify what to save to file for debugging*/
    MISREG_CFG_T misreg;
    int npowfs;      /**<Number of wfs types*/
    int nwfs;        /**<Number of wfs*/
    int nwfsr;       /**<Number of wfs used in reconstruction. =npowfs in glao, =nwfs otherwise*/
    int ndm;         /**<Number of DMs*/
    int nmoao;       /**<Number of different MOAO type*/
    char **surf;     /**<OPD surfaces*/
    int nsurf;       /**<Number of OPD surfaces*/
    char **tsurf;    /**<Tilted surfaces, surface, not OPD*/
    int ntsurf;      /**<Number of tilted surfaces*/
    lmat *fdlock;     /**<Records the fd of the seed lock file. if -1 will skip the seed*/
    int nlopowfs;    /**<Number of low order wfs types*/
    lmat *lopowfs;    /**<List of low order powfs*/
    int nhipowfs;    /**<Number of high order wfs types*/
    lmat *hipowfs;    /**<List of high order powfs*/
    int ntrpowfs;    /**<Number of tip/tilt removed wfs types*/
    int ntipowfs;    /**<Number of tip/tilt include wfs types*/
    int nphypowfs;   /**<Number of powfs with local/uplink tip/tilt loop*/
    int nlowfs;      /**<Number of low order wfs.*/
    int nhiwfs;      /**<Number of high order wfs*/
    dmat *dirs;      /**<Collect for beam directions*/
    int dither;      /**<Some WFS is doing dithering*/
    int ilgspowfs;     /**<Index of LGS WFS*/
    int nlgspowfs;   /**<Number of LGS WFS*/
    int itpowfs;       /**<Index of twfs*/
    int idmground;   /**<Index of ground dm. default to 0*/
    int step_lo;     /**<Enabling step for low order wfs*/
    int step_hi;     /**<Enabling step for high order wfs*/
}PARMS_T;
/**
   ARG_T is used for command line parsing.
*/
typedef struct ARG_T{
    int detach;      /**<Detach from the command line and run in background*/
    int override;    /**<Run simulation even if Res_${seed}.done exists*/
    int force;       /**<For start, bypassing scheduler*/
    int *gpus;       /**<Index of GPU to use. -1 to disable*/
    int ngpu;        /**<Number of entries in gpus*/
    int ngpu2;       /**<Number of GPUs to use. Ignore of gpus is set.*/
    int server;      /**<MAOS acting as server*/
    char *dirout;    /**<Result output directory*/
    char *conf;      /**<master .conf file. nfiraos.conf by default. -c to change*/
    char *confcmd;   /**<Additional configuration options supplied in command line.*/
}ARG_T;
PARMS_T* setup_parms(char *main, char *extra, int override);
void setup_parms_gpu(PARMS_T *parms, int *gpus, int ngpu);
void free_parms(PARMS_T *parms);
/*The following are here so that we don't have to include type.h or utils.h */
/*convenient constants. used in utils.c */
typedef enum T_TYPE{
    T_PLOC=0,
    T_ALOC,
    T_XLOC,
    T_ATM,
}T_TYPE;
void plotdir(char *fig, const PARMS_T *parms, double totfov, char *format,...);
#endif

/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
    double l0;    /**<outer scale*/
    double dx;    /**<sampling of turbulence screens*/
    double *ht;   /**<height of each layer*/
    double *wt;   /**<weight of each layer (relative strength of \f$C_n^2\f$)*/
    double *ws;   /**<wind speed of each layer*/
    double *wddeg;/**<wind direction of each layer*/
    double *size; /**<size of atm in meter, [0 0]: automatic*/
    long *overx;  /**<maximum pixel distance in x direction the beam can be without wrapping*/
    long *overy;  /**<maximum pixel distance in y direction the beam can be without wrapping*/
    map_t **(*fun)(GENSCREEN_T*); /**<Points to the function used to generated atmosphere. (derived)*/
    int nps;      /**<number of phase screens*/
    int wdrand;   /**<randomize wind direction*/
    int iground;  /**<index into the ground layer*/
    int *ipsr;    /**<corresponding reconstruction layer*/
    int nx;       /**<turbulence screen size alog x*/
    int ny;       /**<turbulence screen size along y*/
    int fractal;  /**<1: Use fractal method to generate atmosphere screen.*/
    int evolve;   /**<evolve the atm in additional to frozen flow. developed for
		     fractal since it does not wrap.*/
    int frozenflow;  /**<frozen flow. automatic if closeloop=1*/
    int ninit;    /**<Initial size of the screen in fractal method. >=2*/
    int share;    /**<0: disable sharing of atmosphere using file backend*/
}ATM_CFG_T;
/**
   contains input parameters for the atmospheric reconstruction.  */
typedef struct ATMR_CFG_T{
    double r0z;   /**<r0 at zenith*/
    double r0;    /**<derived from r0z for zenith angle za*/
    double l0;    /**<outer scale*/
    double hs;    /**<height of the high order guide star. derived*/
    double *ht;   /**<height of each layer*/
    double *wt;   /**<weight of each layer (relative strength of \f$C_n^2\f$)*/
    double dx;    /**<baseline sampling (when os=1). matches to high order wfs.*/
    int *os;      /**<over sampling factor of xloc over actuator spacing */
    int nps;      /**<number of phase screens*/
}ATMR_CFG_T;
/**
   contains input parameters about the aperture, like the diameter,
   amplitude map, etc */
typedef struct APER_CFG_T{
    double d;     /**<Telescope aperture diameter*/
    double din;   /**<Telescope inner blocking diameter*/
    double dx;    /**<sampling of aperture evaluation grid aper_locs*/
    double rotdeg;/**<pupil rotation in degree*/
    char *fnamp;  /**amplitude maps. expected to be square or rectangular mxn, with 0 at
		     [m/2,n/2] (count from 0)*/
    double *misreg;/**<misregistration of the pupil along x and y. Shifts the
		      amplitude map and atmosphere.*/
    int ismisreg; /**<true if misreg contains nonzero numbers*/
    char *pupmask;/**<The pupil cold stop*/
}APER_CFG_T;
/**
   contains input parameters for laser launch telescope
*/
typedef struct LLT_CFG_T{
    double d;      /**<LLT clear aperture diameter*/
    double widthp; /**<Gaussian beam width percentage of d*/
    char *fnrange; /**<File contains range to sodium layer*/
    char *fn;      /**<File contains sodium profile*/
    char *fnsurf;  /**<Surface OPD error*/
    double *ox;    /**<location x of LLT center wrt telescope aperture center*/
    double *oy;    /**<see ox.*/
    double *misreg;
    int *i;        /**<Index into llt for this iwfs.*/
    int smooth;    /**<smooth the sodium profile or not*/
    int n;         /**<number of launch telescopes in this powfs*/
    int colprep;   /**<starting column to use in fn for ETF in preparation of
		      matched filter*/
    int colsim;    /**<starting column to use in fn for ETF in simulation*/
    int colsimdtrat;/**<change to next sodium profile during simulation every
		       colsimdtrat time step*/
} LLT_CFG_T;
/**
   contains input parameters for each type of wfs (powfs).
*/
typedef struct POWFS_CFG_T{
    double *wvl;   /**<list of wavelength*/
    char *piinfile;/**<input averaged pixel intensities for matched filter. NULL
		      to disable*/
    char *sninfile;/**<Speckle noisy input file. NULL to disable. not used*/
    double hs;     /**<height of guide star*/
    double saat;   /**<subaperture area (normalized) threshold to drop subaperture.*/
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
    double dx;      /**<sampling of opd points in each subaperture. usually
		       matches atmosphere sampling for LGS. may be coraser for NGS.*/
    double pixtheta;/**<size of pixel pitch along x or y in radian.*/
    double pixoffx; /**<offset of image center from center of detector*/
    double pixoffy; /**<see pixoffx*/
    double sigscale;/**<scale the signal level for simulation.*/
    double siglev;  /**<signal level. will be override by wfs.siglev is specified.*/
    double *wvlwts; /**<weights for each wavelength. can be overriden by wfs.wvlwts.*/
    struct LLT_CFG_T *llt;/**<configuration for LLT*/
    char* fnllt;    /**<filename of LLT configuration. empty means no llt.*/
    int trs;        /**<tip/tilt removal flag. True for LGS, False for NGS*/
    int dfrs;       /**<differential focus removal flag. True for LGS, False for NGS*/
    int lo;         /**<whether this is a low order wfs. False for LGS, True for NGS*/
    int skip;       /**<skip in high order tomography, for split tomo (derived parameter)*/
    int psol;       /**<Compute pseudo open loop gradients (derived parameter)*/
    int *wfs;       /**<array of wfs belongs to this powfs*/
    int *wfsind;    /**<wfsind[iwfs] gives the index of the wfs in this powfs group*/
    int nwfs;       /**<number of wfs belonging to this powfs*/
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
    int phytype;    /**<physical optics type. 1: mtch, 2: tcog*/
    int phytypesim; /**<physical optics type for simulation. -1 to follow phytype*/
    int phystep;    /**<frames to start using physical optics. 
		       -  0: means from frame 0.
		       - >0: need to compute GS0 to calculate geometric optics
		       - -1: never, doesn't need to compute DTF
		    */
    int usephy;     /**<whether physical optics is used at all during
		       simulation.(derived parameter)*/
    double mtchcrx; /**<if >0 use constrained mtch for x for this amount of pixels*/
    double mtchcry; /**<if >0 use constrained mtch for y for this amount of pixels*/
    int mtchcpl;    /**<use coupling between r/a measure error. useful for LGS with x-y ccd.*/
    int mtchstc;    /**<shift peak in the time averaged short exposure PSF to center using fft.*/
    int mtchscl;    /**<scale subaperture image to have the same intensity as i0. Keep false.*/
    int mtchadp;    /**<Using adaptive matched filter. When the number of pixels
		       in the image brighter than half maximum is more than this
		       value, use constraint. introduced on 2011-02-21.*/
    int hasGS0;     /**<need to compute GS0 (derived parameter)*/
    int noisy;      /**<noisy or not during *simulation* */
    char* misreg;   /**misregistration of the WFS, described by file containing
		       a cell array A of 2xp, where p is #1 or number of wfs
		       belonging to this powfs. The misregistered coordinate is
		       computed as

		       xm(ip)=\sum_{ic}(A{1,ip}(1,ic)*pow(x,A{1,ip}(2,ic))*pow(y,A{1,ip}(3,ic)))
		       ym(ip)=\sum_{ic}(A{2,ip}(1,ic)*pow(x,A{1,ip}(2,ic))*pow(y,A{1,ip}(3,ic)))
		       where ip is ipowfs. ic is index of column in entries of
		       A.*/
    char *ncpa;     /**<none common path aberrations for powfs. Used in
		       perparation of matched filter *and* in simulation. Each
		       file contains 2xn cell array, where n can be 1 or number
		       of wfs belonging to this powfs. The format of each 2x1
		       cell is the same as for surf.*/
    int ncpa_method;/**<Method to correct ncpa.
		       - 0: do nothing.
		       - 1: apply gradient electronic offset. 
		       - 2: apply ncpa to average pixel intensity i0, better than 1 */
    int pistatout;  /**<output time averaged short exposure image.*/
    int pistatstart;/**<time step to compute pistatout*/
    int pistatstc;  /**<1: shift to center using fft method. 0: use geometric gradients.*/
    int psfout;     /**<output time history of low order wfs PSF. never do this for LGS.*/
    int dtrat;      /**<ratio of sample period over fast loop (LGS)*/
    int i0scale;    /**<scale i0 to matched subaperture area.*/
    int *scalegroup;/**<scale group for dm propergation cache.(derived parameter)*/
    int moao;       /**<index into MOAO struct. -1: no moao*/
    int dl;         /**<is diffraction limited. derived from comparing pixtheta
		       with diffraction limited image size.*/

}POWFS_CFG_T;
/**
   contains input parmaeters for each wfs
*/
typedef struct WFS_CFG_T{
    double *wvlwts; /**<Weights of signal value for each wavelength. if not
		       specified in config, will use powfs.wvlwts*/
    double thetax;  /**<x direction*/
    double thetay;  /**<y direction*/
    double siglev;  /**<Total signal value for all wavelength. if not specified
		       in config, will use powfs.siglev*/
    double siglevsim;/**<Signal value used for simulation. (derived parameter)*/
    int powfs;      /**<powfs type*/
}WFS_CFG_T;
/**
   contains input parameters for each deformable mirror.
*/
typedef struct DM_CFG_T{
    double guard;   /**<extra DM actuator rings outside of aper.d*/
    double stroke;  /**<Stroke of DM. OPD goes to \f$\pm\f$ stroke$*/
    double vmisreg; /**<vertical misregistration*/
    double ht;      /**<height conjugation range*/
    double dx;      /**<actuator separation*/
    double offset;  /**<Center-most actuator offset from origin
		       - =0 means there is a act on center. 
		       - 1/2 means no one in the center.*/
    double iac;     /**<Inter-Actuator Coupling coefficient.*/
    double histbin; /**<The bin width for histogram.*/
    int histn;      /**<Number of bins in histogram.*/
    int hist;       /**<Compute histogram of commands of each actuator*/
    int cubic;      /**<use cubic spline. better than linear. need to specify iac*/ 
    int order;      /**<Order of the DM within telescope clear subaperture*/
    int isground;   /**<Is this DM the ground DM (derived)*/
    char *actfloat; /**<file containing floating actuators. nx2 coordinate*/
    char *actstuck; /**<file containing stuck actuators. nx2 coordinate.*/

    int ncache;
    double *dxcache;/**<the sampling of plane to cache dm for each scale
		       group. (derived)*/
    char* misreg;   /**<misregistration of the DM with respect to the pupil,
		       described by file containing a cell array A of 2x1.  The
		       misregistered coordinate is computed as

		       xm(ip)=\sum_{ic}(A{1,ip}(1,ic)*pow(x,A{1,ip}(2,ic))*pow(y,A{1,ip}(3,ic)))
		       ym(ip)=\sum_{ic}(A{2,ip}(1,ic)*pow(x,A{1,ip}(2,ic))*pow(y,A{1,ip}(3,ic)))
		       where ip is ipowfs. ic is index of column in entries of
		       A.*/
    char *hyst;     /**<File containing a matrix that describes the
		       hysterisis. The matrix should have 3 rows, and multiple
		       columns. Each column describe the parameters for a mode,
		       and the rows are: 1) weight 2) alpha, 3) beta. See
		       "Hysteresis Modeling and Simulation" by Luc Gilles*/
}DM_CFG_T;
/**
   contarins input parameters all evaluation directions.  */
typedef struct EVL_CFG_T{
    double *thetax; /**<x Coordinate of evaluation directions*/
    double *thetay; /**<y Coordinate of evaluation directions*/
    double *wt;     /**<weight of each direction*/
    double *wvl;    /**<wavelength for PSF and strehl computation*/
    double *hs;     /**<height of each science object*/
    double *misreg; /**<Misregistration wrt to nominal pupil.*/
    int ismisreg;   /**<Science evl is misregistered*/
    int nwvl;       /**<Number of wavelength*/
    int *psf;       /**<1: participate in psf evaluation.*/
    int *psfr;      /**<1: participate in psf reconstruction telemetry*/
    int npsf;       /**<how many directions we compute psf for*/
    int psfol;      /**<compute Open loop PSF.
		       - 1: on axis only.
		       - 2: all directions and average them.*/
    int rmax;       /**<max radial mode for performance evaluation. 
		       - 0: piston only
		       - 1: piston/tip/tilt.*/
    int nmod;       /**<Number of modes. derived from rmax. (nmax+1)*(nmax+2)/2*/

    int psfhist;    /**<output history of the psf (a lot of storage)*/
    int *psfpttr;   /**<remove p/t/t from psf. 1 number for each evl.*/
    int psfmean;    /**<output time averaged psf*/
    int opdcov;     /**<save covairance of science OPD ,every this time step,
		       for directions where evl.psf is 1*/

    int psfisim;    /**<time step to start psfmean.*/
    int *psfsize;    /**<save this number of pixels of the center of the psf. 1
			number for each wvl.*/
    int *psfgridsize;/**<grid size for FFT to generate PSF. Becareful about FFT
			speed and enough padding. Determines the sampling of the
			generated PSF. 0 or negative for automatic. 1 number for
			each wvl.*/
    int nevl;       /**<Number of evaluation directions. (derived)*/
    int tomo;       /**<evaluate tomography performance.*/
    int indoa;      /**<index of the on axis evluation point.*/
    int *scalegroup;/**<scale group for dm cache. havenumber of nevl*dm
		       elements(derived parameter)*/
    int moao;       /**<index into MOAO struct. -1: no MOAO*/
    int nthread;    /**<number of threads for evaluation parallelization*/
}EVL_CFG_T;

/**
   contains input parameters for wavefront tomography.
*/
typedef struct TOMO_CFG_T{
    double tikcr;    /**<tikhonov regularization.*/
    double minwt;    /**<minimum layer weight allowed. if less than this will force to this.*/
    double iac;      /**<#inter-actuator-coupling in cubic influence function (testing)*/
    double cxxscale; /**<scale the Cxx^-1 term.*/
    double svdthres; /**<Threshold in SVD inversion*/
    int square;      /**<use square/rectangular grid instead of tighter irregular grid*/
    int cone;        /**<use cone coordinate in xloc: keep true*/
    int cxx;         /**<method to compute Cxx^-1. 0: bihormonic approx. 1: inverse psd. 2: fractal*/
    int guard;       /**<guard rings of reconstruction grid xloc*/
    int pos;         /**<over sampling factor of ploc over actuator spacing*/
    int piston_cr;   /**<single point piston constraint. */
    int split;       /**<split tomography type.
			- 0: integrated tomography
			- 1: adhoc split tomography
			- 2: minimum variance split tomography*/
    int ahst_wt;     /**<0: use Wg, 1: using Wa*/
    int ahst_idealngs;/**<ideal correction on NGS modes. For skycoverage preprocessing.*/
    int ahst_rtt;    /**<remote tip/tilt in high order DM fit output in split mode*/
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
    int windest;     /**<estimate wind. \todo finish implement it.*/
    int windshift;   /**<shift opdr using wind velocity (from input if windest=0)*/
    int cubic;       /**<cubic influence function in tomography (testing)*/
    int ninit;       /**<like atm.ninit, the initial screen to generate from covariance directly*/
}TOMO_CFG_T;
/**
   contains input parameters for deformable mirror fitting.
*/
typedef struct FIT_CFG_T{
    double *thetax;  /**<x Coordinate of DM fitting directions. */
    double *thetay;  /**<y Coordinate of DM fitting directions. */
    double *wt;      /**<weight of each direction*/
    double *ht;      /**<height of each direction*/
    double tikcr;    /**<tikhonov regularization*/
    double svdthres; /**<Threshold in SVD inversion*/
    int actslave;    /**<slaving constraint for non-active actuators. Useful in CBS method*/
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
}FIT_CFG_T;

/**
   contains input parameters for simulation, like loop gain, seeds, etc.
*/
typedef struct SIM_CFG_T{
    double dt;       /**<sampling period. 1/800*/
    double za;       /**<zenith angle in radian*/
    int start;       /**<time step to start simulation. 0*/
    int end;         /**<time step to stop simulation. exclusive*/
    int *seeds;      /**<simulation seeds*/
    int nseed;       /**<How many simulation seed*/
    int nthread;     /**<Number of threads to run the simulation*/
    int closeloop;   /**<closed loop or open loop*/
    char *gtypeII_lo;/**<contains 3x1 or 3xnmod type II gains.*/
    char *wspsd;     /**<Telescope wind shake PSD input. Nx2. First column is
			freq in Hz, Second column is PSD in rad^2/Hz.*/
    //control
    double *apdm;    /**<servo coefficient for high order dm.  A is command. e is
			error signal. at time step n, the command is updated by
			A(n)=A(n-1)*apdm(0)+A(n-2)*ap(1)+...+e(n-2)*ep
		     */
    double *apngs;   /**<servo coefficient for ngs modes.*/
    double *apupt;   /**<servo coefficient for for LGS uplink pointing loop.*/
    double epdm;     /**<error gain for DM commands (high order)*/
    double epngs;    /**<error gain for NGS modes (low order)*/
    double epupt;    /**<error gain for uplink pointing*/
    double dpupt;    /**<derivative tracking for uplink pointer. keep 0 to disable*/
    double epfocus;  /**<error gain for LGS focus tracking with zoom optics*/
    double lpfocus;  /**<parameter for low pass filter of LGS focus tracking with offset*/
    double fov;      /**<The diameter of total fov in arcsec*/
    int napdm;       /**<number of entries in apdm*/
    int napngs;      /**<number of entries in apngs*/
    int napupt;      /**<number of entries in apupt*/
    int mffocus;     /**<method for focus tracing.
			- 0: no focus tracking.
			- use CL grads + DM grads - Xhat grad for LGS and NGS.
			- use CL grads + DM grads - Xhat grad for LGS; */
    int uptideal;    /**<ideal compensation for uplink pointing*/
    int servotype_hi;/**<servo type for high order loop. 1: simple integrator*/
    int servotype_lo;/**<servo type for low order loop. 1: simple integrator. 2: type II*/
    int cachedm;     /**<cache dm shape on fine sampled grid matched WFS or Science grid*/
    int cachesurf;   /**<cache surface on fine sampled grid matched WFS or Science grid*/
    int fuseint;     /**<fuse the high and low order integrators in split tomography */
    int skysim;      /**<1: we are doing skycoverage preprocessing*/
    int recon;       /**<reconstruction method. 0: minimum variance, 1: least square*/
    int glao;        /**<reconstruction in GLAO mode: average WFS measurements*/
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
    int dmttcast;    /**<derived: cast tip/tilt from DM commands to study saturation or
			histogram and then add back*/
    int dmclip;      /**<drived: Need to clip actuators*/
}SIM_CFG_T;
/**
   Parameters for Cn square estimation.
*/
typedef struct CN2EST_CFG_T{
    int *pair;       /**<If non empty, paris of WFS to use for cn2
			estimation. Empty: disable cn2 estimation*/
    int npair;       /**<Derived: number of entries in cn2pair*/
    int step;        /**<do cn2 estimation every this time step*/
    int reset;       /**<reset the accumulated cn2 after every cn2step.*/
    int tomo;        /**<update tomography parameters if non zero*/
    int keepht;      /**<>0: use the layer ht specified by atmr.ht. 2: also do slodar
			directly on these layers.*/
    int nhtomo;      /**<number of layers to feed into reconstructor. only
			effective if keepht=0*/
    int moveht;      /**<1: move the ht used for reconstructor to near strongest
			layers. only effective if keepht=0.*/
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
    int psol;        /**<test add dm command offseted by 1 frame in the future to psol grad*/
    int wamethod;    /**<method to compute wa for ngsmod removal.*/
    int atm;         /**<test special atmosphere*/
    int keepshm;     /**<keep the atmospehre in the shared memory.*/
    int mvstlimit;   /**<Limit number of modes controled on MVST*/
    int annular_W;   /**<Define the W0/W1 on annular aperture instead of circular*/
    int *tomo_maxit; /**<if not empty, will study these maxits in open loop*/
    int ntomo_maxit; /**<Number of elements in tomo_maxit*/
    int tomo_hxw;    /**<1: Force use hxw always instead of ray tracing from xloc to ploc.*/
    int parallel;    /**<The parallel scheme. 1: fully parallel. 0: do not parallel the big loop (sim, wfsgra,d perfevl)*/
    int splitlrt;    /**<1: use LGS low rank terms in split tomography.*/
    int ecovxx;      /**<save the xx used to calculate ecov in psfr.*/
    int useopdr;     /**<use opdr in psf reconstruction*/
    int force;       /**<Force run even if Res_${seed}.done exists*/
    int usegwr;      /**<GA/GX method: 0: GP, 1: GS0*/
    int gpu_wfs;     /**<Use GPU for WFS*/
    int gpu_evl;     /**<Use GPU for evaluation*/
    int gpu_tomo;    /**<Use GPU for tomography*/
    int gpu_fit;     /**<Use GPU for DM fitting*/
    int dxonedge;    /**<have points on the edge of subapertures.*/
}DBG_CFG_T;
/**
   contains input parameters for each MOAO type.
*/
typedef struct MOAO_CFG_T{
    int order;       /**<Order of this MOAO*/
    int cubic;       /**<Whether use cubic influence function*/
    double iac;      /**<Inter-actuator-coupling for cubic influence function*/
    double stroke;   /**<Stroke of the MOAO DM*/
    int actslave;    /**<Do we do actuator slaving*/
    int lrt_ptt;     /**<Piston/tip/tilt constraint*/
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
    int mvst;        /**<load MVST mvst_U and mvst_FU. see recon.c*/
    int GS0;         /**<if 1, load GS0 from powfs%d_GS0.bin.gz*/
    int tomo;        /**<if 1, load tomo matrix*/
    int fit;         /**<if 1, load fit matrix*/
    int W;           /**<if 1, load W0, W1*/
    int i0;          /**<if 1, load i0 for powfs*/
}LOAD_CFG_T;
/**
   contains input parameters for saving variables.
*/
typedef struct SAVE_CFG_T{
    int all;         /**<save absolutely everything. mainly for debugging*/
    int setup;       /**<save preparation matrices*/
    int recon;       /**<save reconstructor information. large*/
    int mvst;        /**<MVST computation intermediate matrices*/

    //run time special ones that need extra computation
    int dmpttr;      /**<save p/t/t removed dm act cmds at each time step*/

    //run time
    int atm;         /**<save atmosphere*/
    int run;         /**<save run time informaton for each time step*/
    int opdr;        /**<save reconstructed OPD on XLOC for each time step*/
    int opdx;        /**<save ATM propagated to XLOC for each time step*/
    int dm;          /**<save computed DM actuator commands for each time step*/
    int evlopd;      /**<save science OPD for each time step*/

    /*for WFS. Converted from scalar or vector input.
      Scalar input: 1: both, 2: high order only, 3: lo only
      Vector input: Equal to the number of WFS.
    */
    int *wfsopd;      /**<save WFS OPD:*/
    int *ints;        /**<save WFS subaperture image*/
    int *grad;        /**<save WFS gradients*/
    int *gradgeom;    /**<save WFS geometric gradient during physical optics simu*/
    //The following are derived from above.
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
    int *gcov;       /**<size of 2*ngcov, specifying wfs for each pair*/
}SAVE_CFG_T;
/**
   is a wrapper of all _CFG_T data types.
*/
typedef struct PARMS_T{
    ATM_CFG_T    atm;   /**<atmospheric parameters*/
    ATMR_CFG_T   atmr;  /**<information about reconstructed atm*/
    APER_CFG_T   aper;  /**<aperture parameters*/
    TOMO_CFG_T   tomo;  /**<tomography parameters*/
    FIT_CFG_T    fit;   /**<DM fit parameters*/
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
    LOAD_CFG_T   load;  /**<Specify what matrices to load for debugging*/
    SAVE_CFG_T   save;  /**<Specify what to save to file for debugging*/
    int npowfs;      /**<Number of wfs types*/
    int nwfs;        /**<Number of wfs*/
    int nwfsr;       /**<Number of wfs used in reconstruction. =npowfs in glao, =nwfs otherwise*/
    int ndm;         /**<Number of DMs*/
    int nmoao;       /**<Number of different MOAO type*/
    char **surf;     /**<OPD surfaces*/
    int nsurf;       /**<Number of OPD surfaces*/
    char **tsurf;    /**<Tilted surfaces, surface, not OPD*/
    int ntsurf;      /**<Number of tilted surfaces*/
    int *fdlock;     /**<Records the fd of the seed lock file. if -1 will skip the seed*/
    int pause;       /**<Pause at the end of every time step*/
    int force;       /**<For start, bypassing scheduler*/
}PARMS_T;
/**
   ARG_T is used for command line parsing.
*/
typedef struct ARG_T{
    int detach;      /**<Detach from the command line and run in background*/
    int force;       /**<For start, bypassing scheduler*/
    int nthread;     /**<Number of threads*/
    int pause;       /**<pause at the end of every time step*/
    int gpu;         /**<Index of GPU to use. -1 to disable*/
    char *dirout;    /**<Result output directory*/
    char *conf;      /**<master .conf file. nfiraos.conf by default. -c to change*/
    char *confcmd;   /**<Additional configuration options supplied in command line.*/
}ARG_T;
PARMS_T* setup_parms(ARG_T *arg);
void free_parms(PARMS_T *parms);
//The following are here so that we don't have to include type.h or utils.h
//convenient constants. used in utils.c
typedef enum T_TYPE{
    T_PLOC=0,
    T_ALOC,
    T_XLOC,
    T_ATM,
}T_TYPE;
void create_metapupil(const PARMS_T *parms, double ht, double dx,
		      double offset,long* nx, long* ny, double *ox, double *oy, 
		      double **map,double guard, long nin, int pad,int square);
map_t *create_metapupil_wrap(const PARMS_T *parms, double ht,double dx,
			     double offset,double guard,long nin, 
			     int pad,int square);
void plotdir(char *fig, const PARMS_T *parms, double totfov, char *format,...);
#endif

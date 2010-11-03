/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef __AOS_TYPES_H__
#define __AOS_TYPES_H__
#if USE_PTHREAD > 0
#include <pthread.h>
#endif
#include <fftw3.h>
/**
   \file maos/types.h
   A few run time structs
*/
/**
   contains the data associated with the aperture and performance
   evaluation.
 */
typedef struct APER_T{
    loc_t *locs;         /**<PLOCS in laos. the fine sampled grid on aperture
		         for peformance evaluation.*/
    locstat_t *locs_stat;/**<record the starting location of each column in
			    locs, used to accelerate accphi*/
    double *amp;         /**<amplitude map defined on locs, if exists. sum to 1. for
			    performance evaluation*/
    double *amp1;        /**<amplitude map defined o locs, maximum is 1.*/
    map_t *ampground;    /**<The input amplitude map on ground level read from file.*/
    dmat *mod;           /**<modal columne vectors if parms->evl.nmax>1*/
    dmat *mcc;           /**<piston/tip/tilt mode cross-coupling.*/
    dmat *imcc;          /**<inverse of mode cross coupling for evaluations.*/
    double ipcc;         /**<piston term in imcc.*/
    double sumamp2;      /**<sum of amplitude squared*/
    int **embed;         /**<Embedding index for PSF computing, one per wvl*/
    int *nembed;         /**<dimension of embed.*/
    double fcp;          /**<piston correction in focus term.*/
}APER_T;
/**
   contains the data associated with a detector transfer function for a
   subaperture. The PSF is computed as
   \f$\textrm{PSF}=\frac{1}{N^2\sum(\textrm{amp}^2)}|\textrm{fftshift}\mathcal{F}[A
   \exp(-\frac{2\pi}{\lambda}\textrm{opd})]|^2\f$.  The subaperture image is
   computed as
   \f$I=\textrm{si}*\mathcal{F}^{-1}[\mathcal{F}[\textrm{PSF}\times\textrm{nominal}]]\f$
*/
typedef struct DTF_T{
    ccell *nominal;      /**<The FFT of the pixel functions*/
    spcell *si;          /**<The pixel selection*/
    cmat *Ux;            /**<Special frequency vector along x*/
    cmat *Uy;            /**<Special frequency vector along y*/
    int fused;           /**<Whether the DTF has been fused to ETF*/
}DTF_T;
/**
   contains the data associated with an elongation transfer function for a
subaperture.  */
typedef struct ETF_T{
    ccell *p1;          /**<Store the ETF along radial direction when radrot==1*/
    ccell *p2;          /**<Store the 2D ETF when radrot==0*/
}ETF_T;
/**
   contains the data associated with a LLT uplink path.
 */
typedef struct LOTF_T{
    loc_t *loc;          /**<The grid that defines the LLT pupil*/
    pts_t *pts;          /**<The LLT lower left grid point location*/
    double *amp;         /**<The amplitude defined on loc*/
    dcell *mcc;          /**<modal cross coupling matrix*/
    dcell *imcc;         /**<inverse of imcc*/
}LOTF_T;
/**
   contains the intensity statistics assiciated with a certain powfs for
physical optics wfs. */
typedef struct INTSTAT_T{
    ccell *lotf;        /**<llt otf*/
    ccell **otf;        /**<short exposure OTF. time consuming to calculate. */
    dcell **sepsf;      /**<short expsoure PSF.*/
    dcell *i0;          /**<short exposure image. nsa x nllt*/
    dcell *gx;          /**<gradient of i0 along x*/
    dcell *gy;          /**<gradient of i0 along y*/
    dmat *i0sum;        /**<sum of i0*/
    dcell *mtche;       /**<mtched filter opeator too apply to gradients*/
    dcell *saneaxy;     /**<computed sanea along xy. (rad^2)*/
    dcell *saneaixy;    /**<computed sanea inverse along xy (rad^-2)*/
    dcell *saneara;     /**<computed sanea along ra (rad^2) for radian coord ccd*/
    dcell *sanea;       /**<SANEA of matched filter output. ra or xy (rad^2).*/
    int notf;           /**<number of otf; 1 unless there is ncpa.*/
    int nsepsf;         /**<number of sepsf; usually 1.*/
    int nmtche;         /**<number of matched filters. 1 or nwfs of this powfs.*/
}INTSTAT_T;
/**
   contains the data associated with a certain type of WFS. not
necessarily physical optics WFS.x */
typedef struct POWFS_T{
    loc_t *saloc;       /**<lower left corner of the subaperture*/
    pts_t *pts;         /**<records lower left-most point of each sa in a regular grid*/
    loc_t *loc;         /**<explicit pts in a regular grid.*/
    loc_t **locm;       /**<mis-registered loc, if any.*/
    dcell *ncpa;        /**<NCPA OPDs to add to WFS OPD during simulation.*/
    dcell *ncpa_grad;   /**<NCPA grads to subtract from measurement. */
    DTF_T *dtf;         /**<array of dtf for each wvl*/
    dmat *sodium;       /**<Loaded and downsampled sodium profile from the file.*/
    ETF_T *etfprep;     /**<ETF for computing short exposure matched filter.*/
    ETF_T *etfsim;      /**<ETF for simulation.*/
    dmat *focus;        /**<additional focus error. (llt->fnrange)*/
    LOTF_T *lotf;       /**<uplink aperture*/
    dcell *mcc;         /**<P/T/T cross coupling matrix with amplitude weighting.*/
    dcell *imcc;        /**<inverse of mcc.*/
    dcell *srot;        /**<subaperture rotation wrt LLT*/
    dcell *srsa;        /**<subaperture distance wrt LLT*/
    dmat *srsamax;      /**<max of srsa for each llt.*/
    dcell *sprint;      /**<which subapertures to print sanea*/
    INTSTAT_T *intstat; /**<matched filter i0 and its derivative.*/
    dsp *GS0;           /**<gtilt (average gradient)*/
    dsp *ZS0;           /**<ztilt (zernike best fit tip/tilt)*/
    double *amp;        /**<amplitude map, max at 1.*/
    double **ampm;      /**<amplitude map on misregistered grid, locm. used for simulation*/
    double areascale;   /**<1./max(area noramlized by dsa*dsa)*/
    dmat *dtheta;       /**<sampling of the imaging fft grid. wvl/(embfac*dxsa);*/
    dcell *bkgrnd;      /**<wfs background image. from parms->powfs[ipowfs].bkgrndfn.*/
    int nwfs;           /**<number of wfs belonging to this wfs*/
    int npts;           /**<number of total opd points in all subaps*/
    int namp;           /**<number of amplitude maps*/
    int pixpsax;        /**<number of detector pixels along x*/
    int pixpsay;        /**<number of detector pixels along y*/
    int ncompx;         /**<Dimension of FFT for subaperture imaging along x*/
    int ncompy;         /**<Dimension of FFT for subaperture imaging along y*/
    int nlocm;          /**<number of misregistered loc. 1 or nwfs of this powfs.*/
}POWFS_T;

/**
   NGS mode and reconstructors in ad hoc split tomography. Only useful with 1 or 2 DMs.
 */
typedef struct NGSMOD_T{
    double hs;      /**<height of LGS*/
    double ht;      /**<height of upper DM.*/
    double scale;   /**<(1-ht/hs)^-2*/
    double aper_fcp;/**<piston term in focus in plocs.*/
    dmat *MCC;      /**<cross coupling of the NGS modes. 2x2 for 1 dm. 5x5 for 2 dms*/
    dmat *MCC_OA;   /**<cross coupling of the modes for on axis direction only.*/
    dmat *IMCC_TT;  /**<inv of cross coupling of tip/tilt modes only.*/
    dmat *IMCC;     /**<inv of MCC.*/
    dcell *GM;      /**<ngsmod vector to gradient operator*/
    dcell *Rngs;    /**<NGS reconstructor from NGS grad to NGS mod vec. pinv of GM*/
    dcell *Pngs;    /**<Project DM command to NGS modes */
    dcell *Mdm;     /**<DM vector for the modes*/
    dcell *Ptt;     /**<Invidual DM tip/tilt removal.*/
    spcell *Wa;     /**<Aperture weighting. Ha'*W*Ha. It has zeros in diagonal. Add tikholnov*/
    int nmod;       /**<nmod: 5 for 2 dm, 2 for 1 dm.*/
}NGSMOD_T;
/**
   contains data for Fourier Domain Preconditioner.
 */
typedef struct FDPCG_T{
    csp *Minv;     /**<inverse of fourier domain preconditioner matrix M.*/
    ccell *Mbinv;  /**<block version of Minv. (in permuted order)*/
    long *perm;    /**<Permutation vector to get block diagonal matrix*/
    long nxtot;    /**<Total number of reconstructed phase points*/
    cmat *xhat;    /**<Intermediate matrices*/
    cmat *xhat2;   /**<Intermediate matrices*/
    ccell *xhati;  /**<Intermediate matrices*/
    ccell *xhat2i; /**<Intermediate matrices*/
    double *scale; /**<Scaling factor for each layer*/
}FDPCG_T;

/**
   contains MOAO related data
*/
typedef struct MOAO_T{
    int used;         /**<Whether this MOAO is used or not*/
    loc_t *aloc;      /**<Actuator grid*/
    spcell *HA;       /**<Propagator from this aloc to PLOC*/
    dcell *NW;        /**<null modes and constraints*/
    dmat *W1;         /**<Weighting matrix on PLOC. same as recon->W1*/
    dsp *W0;          /**<Weighting matrix on PLOC. same as recon->W0*/
    spcell *actslave; /**<Slaving operator for actuators not illuminated*/
}MOAO_T;
/**
   A convenient wrap of the data to embed into MUV_T
*/
typedef struct INVPSD_T{
    dcell *invpsd;    /**<points to recon->invpsd*/
    ccell *fftxopd;   /**<points to recon->fftxopd.*/
    long **xembed;    /**<points to recon->xloc_embed*/
}INVPSD_T;
typedef struct CN2EST_T CN2EST_T;
/**
   contains data related to wavefront reconstruction and DM fitting.  */
typedef struct RECON_T{
    double r0;         /**<r0 used in reconstruction. may get updated in cn2 estimation*/
    double l0;         /**<l0 used in reconstruction. may get updated in cn2 estimation*/
    dmat *ht;          /**<height of the layers to do tomography.*/
    dmat *wt;          /**<weight of the layers to to tomography. may get updated in cn2 estimation*/
    dmat *os;          /**<over sampling of the layers.*/
    dmat *dx;          /**<sampling in meter of the layers*/

    loc_t *ploc;       /**<reconstructed pupil*/
    long ploc_nx;      /**<Size of pmap used to build ploc*/
    long ploc_ny;      /**<size of pmap used to build ploc*/

    loc_t **xloc;      /**<reconstructed atmosphere grid.*/
    long *xloc_nx;     /**<size of xmap used to build xloc*/
    long *xloc_ny;     /**<size of xmap used to build xloc*/
    long **xloc_embed; /**<index array to embed opd on xloc to square array of xloc_nx*xloc_ny*/
    dcell *xmcc;       /**<used for tip/tilt removal from tomographic screens.*/

    loc_t **aloc;      /**<actuator grid*/
    long *aloc_nx;     /**<size of amap used to build aloc*/
    long *aloc_ny;     /**<size of amap used to build aloc*/

    dcell *aimcc;      /**<used for tip/tilt removal from DM commands.*/
    dsp  *W0;          /**<ploc weighting for circle of diam aper.d*/
    dmat *W1;          /**<ploc weighting for circle of diam aper.d*/
    dmat *fitwt;       /**<fit weighting in each direction.*/
    spcell *L2;        /**<Laplacian square regularization.*/
    spcell *L2save;    /**<save old L2 to update the tomography matrix.*/
    dcell *invpsd;     /**<inverse of atm PSD, to apply Cxx^-1 in Fourier space.*/
    ccell *fftxopd;    /**<temporary storage array to apply invpsd*/

    FDPCG_T *fdpcg;    /**<fdpcg preconditioner data.*/
    spcell *GP;        /**<Gradient operator from HXW. GX=GP*H.*/
    spcell *HXW;       /**<Ray tracing operator from xloc to ploc for all WFS.*/
    spcell *HXWtomo;   /**<Like GXtomo*/
    spcell *GX;        /**<Gradient operator for all WFS from each layer of xloc*/
    spcell *GXtomo;    /**<GX for tomography. excluding NGS in split tomography*/
    spcell *GXhi;      /**<GX for high order WFS.*/
    spcell *GXlo;      /**<GX for low order WFs*/
    spcell *GXfocus;   /**<GX used for focus tracking.*/
    dcell *GXL;        /**<dense GX for low order WFS in MV split tomography.*/
    dcell *MVRngs;     /**<NGS recon for MV split tomography*/
    dcell *MVModes;    /**<MVST Modes (svd'ed)*/
    spcell *GA;        /**<actuator to wfs grad.*/
    spcell *GAlo;      /**<GA of low order WFS.*/
    spcell *HXF;       /**<ray tracing propagator from xloc to ploc for fitting directions.*/
    spcell *HA;        /**<ray tracing from aloc to ploc for fitting directions.*/
    dcell *TT;         /**<TT modes for LGS WFS*/
    dcell *PTT;        /**<pinv of TT for tt removal from LGS gradients*/
    dcell *DF;         /**<Differential focus modes for LGS wfs*/
    dcell *TTF;        /**<Concatenation of TT and DF*/
    dcell *PTTF;       /**<pinv of TTF*/
    spcell *ZZT;       /**<single point piston constraint in tomography.*/
    dcell *NW;         /**<null modes for DM fit.*/
    spcell *actslave;  /**<force slave actuators to have similar value to active neighbor ones.*/
    double fitscl;     /**<strength of fitting FLM low rank terms (vectors)*/
    spcell *saneai;    /**<inverse of sanea^2 in radian^-2*/
    dmat *neam;        /**<subaperture averaged nea for each wfs*/
    double neamhi;     /**<average of neam for high order wfs.*/
    MUV_T RR;          /**<tomography right hand side matrix, solve RL*x=RR*y*/
    MUV_T RL;          /**<tomography left hand side matrix*/
    MUV_T FR;          /**<DM fit right hand size matrix, solve FL*x=FR*y*/
    MUV_T FL;          /**<DM fit left hand size matrix*/
    MOAO_T *moao;      /**<for MOAO DM fitting*/

    //For focus tracking.
    dcell *RFlgs;      /**<focus reconstruction from each LGS grad*/
    dcell *RFngs;      /**<focus reconstruction from NGS grad.*/
    dcell *RFtomo;     /**<focus recon from reconstructed X.*/
    int *skipwfs;      /**<skip the wfs in split tomo. (NGS)*/
    NGSMOD_T *ngsmod;  /**<ngs mod in ad hoc split tomography.*/
    CN2EST_T *cn2est;  /**<For Cn2 Estimation*/
    int *xdim;         /**<length of each xloc layer.*/
    int *gdim;         /**<length of gradient vector for each wfs.*/
    int lowfs_gtilt;   /**<=1 if any low order wfs use gtilt in recon/simu*/
    int npsr;          /**<number of reconstructor phase screens.*/
    int ndm;           /**<number of DMs;*/
    int has_ttr;       /**<whether there is any tip/tilt removed WFS*/
    int has_dfr;       /**<whether there is any differential focus removed WFS*/
    int nthread;       /**<number of threads in reconstruction.*/
    int calc_invpsd;   /**<Whether we need to compute inverve of psd*/
}RECON_T;

/**
   contains internal data for a type II integrator
*/
typedef struct TYPEII_T{
    dcell *lead;       /**<Output of lead filter*/
    dcell *errlast;    /**<Record error in last time*/
    dcell *firstint;   /**<First integrator output*/
}TYPEII_T;
typedef struct SIM_SAVE_T{
   /*cellarr to save telemetry data.*/
    //Deformable mirror.
    cellarr *dmerr_hi;
    cellarr *dmint_hi;
    cellarr *dmfit_hi;
    cellarr *dmint;
    cellarr *dmpttr;
    cellarr *dmreal;
    //Low order modes
    cellarr *Merr_lo;
    cellarr *Mint_lo;
    cellarr *opdr;
    cellarr *opdx;
    //science
    cellarr **evlopdcl;
    cellarr **evlopdol;

    cellarr **wfsopd;
    cellarr **wfslltopd;
    //gradients
    cellarr **gradcl;
    cellarr **gradgeom;
    cellarr **gradnf;
    cellarr **gradpsol;
    cellarr **intsny;
    cellarr **intsnf;
    cellarr **moao_evl;
    cellarr **moao_wfs;
}SIM_SAVE_T;
/**
   contains all the run time data struct.
*/
typedef struct SIM_T{
    map_t **atm;       /**<fine sampled simulation turbulence screens*/
    map_t **cachedm;   /**<grid cache dm actuator to a finer sampled screen. for
			  fast ray tracing to WFS and aper*/
    int (*pcachedm)[2];/**<information about cachedm struct.*/
    PROPDATA_T *cachedm_propdata; /**<wrapped data for ray tracing from aloc to cachedm*/
    thread_t *cachedm_prop; /**<wrapped cachedm_propdata for threading*/
    dmat *winddir;     /**<input wind direction*/
    dcell *surfwfs;    /**<additional OPD surface for each WFS.*/
    dcell *surfevl;    /**<additional OPD surcace for each evaluation field*/
    rectmap_t **tsurf; /**<input tilted M3 surface*/
    map_t **surf;      /**<input surface: M1, M2 or else. common to all wfs and science field.*/

    /*The following has a cell for each. wfs*/
    dcell *gradcl;     /**<cl grad output.*/
    dcell *gradpsol;   /**<WFS pseudo open loop grad output*/
    dcell *gradacc;    /**<accumulate gradident for dtrat>1*/

    dcell *opdr;       /**<reconstructed OPD defined on xloc in tomography output.*/
    dcell *opdrmvst;   /**<average of opdr for MVST*/
    dcell *opdevl;     /**<evaluation opd for selected directions.*/
    dmat  *opdevlground;/**<evaluation opd for ground layer turbulence to save ray tracing.*/
    
    ccell *opdrhat;    /**<For wind estimation (testing)*/
    ccell *opdrhatlast;/**<for wind estimation.(testing)*/
    dmat *windest;     /**<estimated wind velocity.(testing)*/
    spcell *windshift; /**<operator to do wind shift on opdr*/

    dcell **dmpsol;    /**<time averaged dm command (dtrat>1) for psol grad*/
    dcell **ints;      /**<WFS subaperture images.*/
    ccell **wfspsfout; /**<output WFS PSF history.*/
    cellarr **wfspsfoutcellarr; /**<special file to save wfs psf history*/
    cellarr **ztiltoutcellarr;  /**<special file to save zernike wfs tilt history*/

    dcell **pistatout; /**<WFS time averaged tip/tilt removed PSF*/
    dcell *sanea_sim;  /**<accumulate effective sanea during simulation.*/

    /*save results to file using mmap*/
    dcell *clep;       /**<CL error per direction.*/
    dcell *clmp;       /**<CL mode coefficient per direction.*/
    dcell *cleptomo;   /**<CL error for tomography compensation perdirection.*/
    dcell *clmptomo;   /**<CL mode coefficient for tomography compensation.*/
    dcell *olep;       /**<OL error per direction.*/
    dcell *olmp;       /**<OL mode coefficient per direction.*/
    dmat *ole;         /**<field averaged OL error*/
    dmat *cle;         /**<field averaged CL error*/
    dmat *cletomo;     /**<field averaged tomography error*/
    
    /*Only on split tomography*/
    dmat *clem;        /**<lgs/ngs mod error in split tomography*/
    dcell *clemp;      /**<lgs/ngs mod error per direction. only on-axis is computed.*/
    dmat *corrNGSm;    /**<Correction of NGS mode. (integrator output)*/
    dmat *cleNGSm;     /**<Close loop ngs mods in split tomogrpahy. */
    dcell *cleNGSmp;   /**<(M'*w*phi);*/
    dcell *res;        /**<warping of ole,cletomo,cle,clem for easy saving.*/
    dmat *gtypeII_lo;  /**<gain for type II array.*/

    /*DM commands.*/
    dcell *dmreal;     /**<real dm command*/
    dcell *dmhist;     /**<histogram of dm commands. if dbg.dmhist is 1.*/
    dcell **dmint;     /**<dm integrator. (used of fuseint==1)*/

    /*High order*/
    dcell *dmfit_hi;   /**<direct high order fit output*/
    dcell *dmerr_hi;   /**<high order dm error signal.*/
    dcell **dmint_hi;  /**< integrator for high order (only if fuseint==0)*/

    /*Low order*/
    dcell *Merr_lo;    /**<split tomography NGS mode error signal.*/
    dcell **Mint_lo;   /**<NGS integrator (only used if fuseint==0)*/
    TYPEII_T *MtypeII_lo;  /**<intermediate results for type II/lead filter*/  
    
    /*llt pointing loop*/
    dcell *upterr;     /**<uplink error*/
    dcell *upterrlast; /**<uplink error from last step*/
    dcell *uptreal;    /**<uplink real*/
    dcell **uptint;    /**<uplink integrator output.*/
    dcell *upterrs;    /**<mmaped file to store upterr history*/
    dcell *uptcmds;    /**<mmaped file to store uptcmd history*/

    dcell *focuslpf;   /**<focus tracking low pass filter*/
    dcell *focusint;   /**<focus tracking integrator*/
    
    dcell *evlpsfmean; /**<science field psf time average*/
    dcell *evlpsfolmean;/**<science field OL PSF time averging*/
    dcell *evlpsftomomean;/**<science field psf time average with direct correction from tomography*/
    cellarr **evlpsfhist;   /**<to save time history of science field psf*/
    cellarr**evlpsftomohist;/**<to save time history of science field psf with direct correction from tomography*/

    /* We maintain separate random streams for each purpose, derived from the
       same seed, so that multi-threading will produce same result */
    rand_t *wfs_rand;  /**<random stream for each wfs.*/
    rand_t *atm_rand;  /**<random stream for atmosphere turbulence generation*/
    rand_t *atmwd_rand;/**<random stream for wind direction*/
    rand_t *init;      /**<random stream to initialize other streams*/
    STATUS_T *status;  /**<status report to scheduler.*/
    
    dcell *moao_wfs;   /**<moao DM command for wfs*/
    dcell *moao_evl;   /**<moao DM command for science field*/

    dcell *gcov;       /**<covariance of psuedo open loop gradients.*/

    PROPDATA_T *wfs_propdata_dm;/**<wrap of data for ray tracing from DM in wfsgrad.c*/
    thread_t* wfs_prop_dm;  /**<wrap of wfs_propdata_dm for threaded ray tracing*/
    locstat_t *ploc_stat; /**<statistics of columns in ploc*/

    SIM_SAVE_T *save;
    
#if USE_PTHREAD > 0
    pthread_mutex_t mutex_plot;    /**<mutex for plot*/
    pthread_mutex_t mutex_wfsgrad; /**<mutex for wfsgrad*/
    pthread_mutex_t mutex_perfevl; /**<mutex for perforance evaluation*/
    pthread_mutex_t mutex_cachedm; /**<mutex for cachedm*/
#endif
    int nthread;       /**<number of threads*/
    int wfsgrad_iwfs;  /**<counter for threaded calling*/
    int perfevl_ievl;  /**<counter for threaded calling*/
    int perfevl_iground;/**<index of the layer at ground*/
    int cachedm_i;     /**<counter for threaded calling*/
    int cachedm_n;     /**<length of pcachedm array*/
    int dtrat_hi;      /**<ratio of sampling period over clock of high order wfs*/
    int dtrat_lo;      /**<dtrat of the lower order loop.*/
    int seed;          /**<current running seed.*/
    int isim;          /**<record current simulations step.*/
    //maintain pointer of other structs for use in thread functions.
    const PARMS_T *parms; /**<pointer to parms*/
    const APER_T *aper;/**<pointer to aper*/
    RECON_T *recon;    /**<pointer to recon*/
    POWFS_T *powfs;    /**<pointer to powfs*/
    double dt;         /**<System baseline clock period. 1/800 s*/
    double dtlo;       /**<low order wfs sampling period*/
    double dthi;       /**<high order wfs sampling period*/
    int has_upt;       /**<whether we have uplink pointer loop.*/
    int last_report_time;/**<The time we lasted reported status to the scheduler.*/
}SIM_T;
//convenient constants. used in utils.c
typedef enum T_TYPE{
    T_PLOC=0,
    T_ALOC,
    T_XLOC,
    T_ATM,
}T_TYPE;
#endif

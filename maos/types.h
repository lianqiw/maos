/*
  Copyright 2009,2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
    LOC_T *locs;         /**<PLOCS in laos. the fine sampled grid on aperture
		         for peformance evaluation.*/
    LOCSTAT_T *locs_stat;/**<record the starting location of each column in
			    locs, used to accelerate accphi*/
    double *amp;         /**<amplitude map defined on locs, if exists. sum to 1. for
			    performance evaluation*/
    double *amp1;        /**<amplitude map defined o locs, maximum is 1.*/
    MAP_T *ampground;    /**<The input amplitude map on ground level read from file.*/
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
    LOC_T *loc;          /**<The grid that defines the LLT pupil*/
    PTS_T *pts;          /**<The LLT lower left grid point location*/
    double *amp;         /**<The amplitude defined on loc*/
    dcell *mcc;          /**<modal cross coupling matrix*/
    dcell *imcc;         /**<inverse of imcc*/
}LOTF_T;
/**
   contains the intensity statistics assiciated with a certain powfs
*/
typedef struct INTSTAT_T{
    ccell *lotf;        /**<llt otf*/
    ccell **otf;        /**<short exposure OTF. time consuming to calculate. */
    dcell **sepsf;
    dcell *i0;/*short exposure image. nsa x nllt*/
    //dcell *i0x1, *i0x2, *i0y1, *i0y2;
    dcell *gx;
    dcell *gy;
    dmat *i0sum;
    dcell *mtche;/*mtched filter*/
    dcell *saneaxy;/*computed sanea along xy. (rad^2)*/
    dcell *saneaixy;/*computed sanea inverse along xy (rad^-2)*/
    dcell *saneara;/*computed sanea along ra (rad^2)*/
    dcell *sanea;//SANEA of mtch output. ra or xy (rad^2).
    int notf;//number of otf; 1 unless there is ncpa.
    int nsepsf;
    int nmtche;//number of matched filters. 1 or nwfs of this powfs.
}INTSTAT_T;
typedef struct CN2PAIR_T CN2PAIR_T;
typedef struct CN2EST_T CN2EST_T;
/**
   contains the data associated with a certain type of WFS. not
necessarily physical optics WFS.x */
typedef struct POWFS_T{
    LOC_T *saloc;/*lower left corner of the subaperture*/
    PTS_T *pts;/*records lower left-most point of each sa in a regular grid*/
    LOC_T *loc;/*explicit pts in a regular grid.*/
    LOC_T **locm;/*mis-registered loc*/
    dcell *ncpa;//NCPA OPDs to add to WFS OPD during simulation.
    dcell *ncpa_grad;//NCPA grads to subtract from measurement. 
    DTF_T *dtf;/*array of dtf for each wvl*/
    dmat *sodium;//Loaded and downsampled sodium profile from the file.
    ETF_T *etfprep;//ETF for computing short exposure matched filter.
    ETF_T *etfsim;//ETF for simulation.
    dmat *focus;//additional focus error. (llt->fnrange)
    LOTF_T *lotf;//uplink aperture
    dcell *mcc;//P/T/T cross coupling matrix with amplitude weighting.
    dcell *imcc;//inverse of mcc.
    dcell *srot;//subaperture rotation wrt LLT
    dcell *srsa;//subaperture distance wrt LLT
    dmat *srsamax;//max of srsa for each llt.
    dcell *sprint;//which subapertures to print sanea
    INTSTAT_T *intstat;//matched filter i0 and its derivative.
    dsp *GS0;//gtilt (average gradient)
    dsp *ZS0;//ztilt (zernike best fit tip/tilt)
    double *amp;//amplitude map, max at 1.
    double **ampm;
    double areascale;// 1./max(area noramlized by dsa*dsa)
    dmat *dtheta;//sampling of the imaging fft grid. wvl/(embfac*dxsa);
    dcell *bkgrnd;//wfs background image. from parms->powfs[ipowfs].bkgrndfn.
    int nwfs;/*number of wfs belonging to this wfs*/
    int npts;/*number of total opd points in all subaps*/
    int namp;/*number of amplitude maps*/
    int pixpsax, pixpsay;/*number of detector pixels*/
    int ncompx,ncompy;/*Dimension in dtf->nominal*/
    int nlocm;//number of misregistered loc. 1 or nwfs of this powfs.
}POWFS_T;

/**
   NGS mode in ad hoc split tomography
 */
typedef struct NGSMOD_T{
    double hs;//height of LGS
    double ht;//height of upper DM.
    double scale;//(1-ht/hs)^-2
    double aper_fcp;//piston term in focus in plocs.
    dmat *MCC;//cross coupling of the modes.
    dmat *MCC_OA;//cross coupling of the modes for on axis direction only.
    dmat *IMCC_TT;//inv of cross coupling of tip/tilt modes only.
    dmat *IMCC;//inv of MCC.
    dcell *GM;//ngsmod vector to grad
    dcell *Rngs;//NGS reconstructor from NGS grad to NGS mod vec.pinv of GM
    dcell *Pngs;//Project DM command to NGS modes 
    dcell *Mdm;//DM vector for the modes
    int nmod;//nmod: 5 for 2 dm, 2 for 1 dm.
    dcell *Ptt;//Invidual DM tip/tilt removal.
    spcell *Wa;//Aperture weighting. Ha'*W*Ha. It has zeros in diagonal. Add tikholnov
}NGSMOD_T;
/**
   contains data for Fourier Domain Preconditioner.
 */
typedef struct FDPCG_T{
    csp *Minv;//inverse of fourier domain M.
    ccell *Mbinv;//block version of Minv. (in permuted order)
    long *perm;
    long nxtot;
    cmat *xhat,*xhat2;
    ccell *xhati, *xhat2i;//points to invidual segments of xhat, xhat2
    double *scale;
}FDPCG_T;

/**
   contains MOAO related data
*/
typedef struct MOAO_T{
    int used;//used or not
    //For reconstruction
    LOC_T *aloc;
    spcell *HA;//Propagate from this DM to PLOC
    dcell *NW;//null modes and constraints
    dmat *W1; dsp *W0;//reference from recon.
    spcell *actslave;//Slaving operator
    //For WFS or Performance evaluation
    spcell *HAS;//Propagate to PLOCS
    spcell *HAW;//Propagate to POWFS->PTS
    //MUV_T FR;
    //MUV_T FL;
}MOAO_T;
typedef struct INVPSD_T{
    dcell *invpsd;//points to recon->invpsd
    ccell *fftxopd;//points to recon->fftxopd.
    long **xembed;//points to recon->xloc_embed
}INVPSD_T;
/**
   contains data related to wavefront reconstruction and DM fitting.  */
typedef struct RECON_T{
    double r0;/**<r0 used in reconstruction. may get updated in cn2 estimation*/
    double l0;/**<l0 used in reconstruction. may get updated in cn2 estimation*/
    dmat *ht; /**<height of the layers to do tomography.*/
    dmat *wt; /**<weight of the layers to to tomography. may get updated in cn2 estimation*/
    dmat *os; /**<over sampling of the layers.*/
    dmat *dx; /**<sampling in meter of the layers*/
    LOC_T *ploc;//reconstructed pupil
    long ploc_nx,ploc_ny;//size of pmap used to build ploc;
    LOC_T **xloc;//reconstructed atmosphere grid.
    long *xloc_nx,*xloc_ny;//size of xmap used to build xloc;
    long **xloc_embed;//index array to embed opd on xloc to square array of xloc_nx*xloc_ny
    dcell *xmcc;//used for tip/tilt removal from tomographic screens.
    LOC_T **aloc;//actuator loc
    long *aloc_nx, *aloc_ny;//size of amap used to build aloc;
    dcell *aimcc;//used for tip/tilt removal from DM commands.
    dsp *W0; dmat *W1;//ploc weighting for circle of diam aper.d
    dmat *fitwt;//fit weighting in each direction.
    spcell *L2;//Laplacian square regularization.
    spcell *L2save;/**<save old L2 to update the tomography matrix.*/
    dcell *invpsd;//inverse of atm PSD, to apply Cxx^-1 in Fourier space.
    ccell *fftxopd;//temporary array to apply invpsd
    //csp *fdpcg;//fdpcg preconditioner data.
    FDPCG_T *fdpcg;//fdpcg preconditioner data.
    spcell *GG;/**<Gradient operator from H0. G0=GG*H.*/
    spcell *H0;/**<Ray tracing operator for all WFS from each layer. alternative for G0*/
    spcell *H0tomo;/**<Like G0tomo*/
    spcell *G0;//Gradient operator for all WFS from each layer
    spcell *G0tomo;//G0 for tomography. excluding NGS in split tomography.
    spcell *G0hi;//G0 for high order WFS.
    spcell *G0lo;//G0 for low order WFs;
    spcell *G0focus;//G0 used for focus tracking.
    dcell *G0L;//dense G0 for low order WFS in MV split tomography.
    dcell *MVRngs;//NGS recon for MV split tomography
    dcell *MVModes;//MVST Modes (svd'ed)
    spcell *GA;//actuator to wfs grad.
    spcell *GAlo;//GA of low order WFS.
    spcell *HX;//propagator from xloc to ploc
    spcell *HA;//from aloc to ploc
    dcell *TT, *PTT;//TT modes for LGS wfs and its pseudo inverse.
    dcell *DF;//Diff focus modes for LGS wfs
    dcell *TTF, *PTTF;//combination of TT and Diff focus modes and its pseudo inverse
    spcell *ZZT;//single point piston constraint in tomography.
    dcell *NW;//null modes for DM fit.
    spcell *actslave;//force slave actuators to have similar value to active ones.
    double fitscl;//strength of fitting FLM low rank terms (vectors)
    spcell *saneai;//inverse of sanea^2
    MUV_T RR;//recon, solve RL*x=RR*y
    MUV_T RL;
    MUV_T FR;//fit, solve FL*x=FR*y
    MUV_T FL;
    //For MOAO.
    MOAO_T *moao;
    //For focus tracking.
    dcell *RFlgs;//focus reconstruction for each LGS
    dcell *RFngs;//focus recon from NGS grad.
    dcell *RFtomo;//focus recon from reconstructed X.
    int *skipwfs;//skip the wfs in split tomo. (NGS)
    NGSMOD_T *ngsmod;//ngs mod in ad hoc split tomography.
    CN2EST_T *cn2est;/**<For Cn2 Estimation*/

    //dmat *Putt;//Project t/t out from upper DM.
    //dmat *Puttwt;//Weighting used for Putt.
    int *xdim;//length of each xloc layer.
    int *gdim;//length of gradient vector for each wfs.
    int lowfs_gtilt;//=1 if any low order wfs use gtilt in recon/simu
    int npsr;//number of reconstructor phase screens.
    int ndm;//number of DMs;
    int has_ttr;//whether there is any tip/tilt removed WFS
    int has_dfr;//whether there is any differential focus removed WFS.
    int nthread;//number of threads in reconstruction.
    int calc_invpsd;
}RECON_T;

/**
   contains internal data for a type II integrator
*/
typedef struct TYPEII_T{

    dcell *lead;
    dcell *errlast;
    dcell *firstint;
}TYPEII_T;
/**
   contains all the run time data struct.
*/
typedef struct SIM_T{
    MAP_T **cachedm;/*cache dm actuator to a finer sampled
			screen. for fast ray tracing to WFS and aper*/
    thread_t *cachedm_prop;
    PROPDATA_T *cachedm_propdata;
    int (*pcachedm)[2];//information about cachedm struct.
    MAP_T **cacheatm;//cache atm for LGS. not used.
    spcell *HACT;//From aloc to cached plane (reversed)
    
    MAP_T **atm;//fine sampled simulation turbulence screens
    dmat *winddir;
    dcell *surfwfs;//additional OPD surface for WFS.
    dcell *surfevl;//additional OPD surcace for EVL.
    RECTMAP_T **tsurf;
    MAP_T **surf;
    dcell *gradcl;//cl grad output.
    dcell *gradpsol;//psol grad output
    dcell *gradacc;//accumulate gradident for dtrat>1
    dcell *opdr;//reconstructed OPD on xloc in tomography output.
    dcell *opdrmvst;//average of opdr for MVST
    dcell *opdevl;//evaluation opd for selected directions.
    dmat *opdevlground;//evaluation opd for ground layer turbulence.
    
    ccell *opdrhat, *opdrhatlast;//for wind estimation.(testing)
    dmat *windest;//estimated wind velocity.
    spcell *windshift;//operator to do wind shift on opdr
    dcell **dmpsol;//dtrat averaged dm command for psol grad
    dcell **ints;//WFS subaperture images.
    ccell **wfspsfout;//output WFS PSF history.
    cellarr **wfspsfoutcellarr;//special file to save psfout
    cellarr **ztiltoutcellarr;//special file to save ztilt;

    dcell **pistatout;//WFS time averaged tip/tilt removed PSF..
    dcell *sanea_sim;//accumulate effective sanea during simulation.

    //results
    dcell *clep;//CL error per direction.
    dcell *clmp;//CL mode coefficient per direction.
    dcell *cleptomo;//CL error for tomography compensation perdirection.
    dcell *clmptomo;//CL mode coefficient for tomography compensation.
    dcell *olep;//OL error per direction.
    dcell *olmp;//OL mode coefficient per direction.
    dmat *ole;//field averaged OL error
    dmat *cle;//field averaged CL error
    dmat *cletomo;//field averaged tomography error
    
    //Only on split tomography
    dmat *clem;// lgs/ngs mod error in split tomography
    dcell *clemp;//lgs/ngs mod error per direction. only on-axis is computed.
    dmat *corrNGSm;//Correction of NGS mode. (integrator output)
    dmat *cleNGSm;/*Close loop ngs mods in split tomogrpahy. */
    dcell *cleNGSmp;//(M'*w*phi);
    dcell *res;//warping of ole,cletomo,cle,clem for easy saving.
    dmat *gtypeII_lo;//gain for type II array.

    //DM commands.
    dcell *dmreal;//real dm command
    dcell *dmhist;//histogram of dm commands. if dbg.dmhist is 1.
    dcell **dmint;//dm integrator. (used of fuseint==1)

    //High order
    dcell *dmfit_hi; //direct high order fit output
    dcell *dmerr_hi; //high order dm error signal.
    dcell **dmint_hi; // integrator for high order (only if fuseint==0)

    //Low order
    dcell *Merr_lo; //split tomography NGS mode error signal.
    dcell **Mint_lo;//NGS integrator (only used if fuseint==0)
    //The following several are only used in type II/lead filter
    TYPEII_T *MtypeII_lo;    

    //llt pointing loop
    dcell *upterr;//uplink error
    dcell *upterrlast;//uplink error from last step
    dcell *uptreal;//uplink real
    dcell **uptint;//uplink integrator output.

    //output files
    dcell *upterrs;
    dcell *uptcmds;

    dcell *focuslpf;//focus tracking low pass filter
    dcell *focusint;//focus tracking integrator
    
    dcell *evlpsfmean, *evlpsftomomean, *evlpsfolmean;//science PSF time averging. CL and OL.
    cellarr **evlpsfhist, **evlpsftomohist;//science PSF time history
    /* We maintain separate random streams for each purpose so
       that multi-threading will produce same result */
    struct_rand *wfs_rand;//random stream for each wfs.
    struct_rand atm_rand, atmwd_rand;//random sream for atmosphere, wind direction
    struct_rand init;
    STATUS_T *status;//status report to scheduler.
    
    //Only on MOAO mode:
    dcell *moao_wfs;
    dcell *moao_evl;

    dcell *gcov;//covariance of psol gradients.

    thread_t* wfs_prop_dm;/**<wrap of data for ray tracing from DM in wfsgrad.c*/
    PROPDATA_T *wfs_propdata_dm;
    LOCSTAT_T *ploc_stat;
    int nthread;
#if USE_PTHREAD > 0
    pthread_mutex_t mutex_plot;
    pthread_mutex_t mutex_wfsgrad;
    pthread_mutex_t mutex_perfevl;
    pthread_mutex_t mutex_cachedm;
#endif
    //counters for threaded calling.
    int wfsgrad_iwfs;
    int perfevl_ievl;
    int perfevl_iground;
    int cachedm_i;
    int cachedm_n;//length of pcachedm.
    int dtrat_hi, dtrat_lo;//dtrat of the lower order loop. 
    int seed;//current seed.
    int isim;//current simulations step.
    //maintain pointer of other structs for use in thread functions.
    const PARMS_T *parms;
    RECON_T *recon;
    POWFS_T *powfs;
    const APER_T *aper;
    double dt,dtlo,dthi;
    int has_upt;//whether this is uplink pointer loop.
    int last_report_time;
}SIM_T;
//convenient constants. used in utils.c
typedef enum T_TYPE{
    T_PLOC=0,
    T_ALOC,
    T_XLOC,
    T_ATM,
}T_TYPE;
#endif

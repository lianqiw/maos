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
#ifndef __AOS_TYPES_H__
#define __AOS_TYPES_H__
#include <pthread.h>
#include "../lib/aos.h"
/**
   \file maos/types.h
   A few run time structs
*/
/**
   contains the data associated with the aperture and performance
   evaluation.
 */
typedef struct aper_t{
    loc_t *locs;         /**<Fine sampled grid on aperture for peformance evaluation.*/
    loccell *locs_dm;    /**<Distorted locs when mapped onto DM*/
    dmat *amp;           /**<amplitude map defined on locs, if exists. sum to 1. for performance evaluation*/
    dmat *amp1;          /**<amplitude map defined on locs, maximum is 1. use for plotting.*/
    dmat *mod;           /**<modal column vectors if parms->evl.nmax>1*/
    dmat *mcc;           /*piston/tip/tilt mode cross-coupling for evaluations.*/
    dmat *imcc;          /**<inverse of piston/tip/tilt mode cross-coupling for evaluations.*/
    real ipcc;           /**<piston term in imcc.*/
    real sumamp2;        /**<sum of amplitude squared*/
    locfft_t *locfft;     /**<For computing PSF with FFT*/
    real fcp;            /**<piston correction in focus term.*/
    dcell *opdadd;       /**<All OPD surface for each evaluation direction, used for ray tracing*/
    dcell *opdbias;      /**<NCPA OPD surface for each NCPA direction, used for dm_ncpa calibration.*/
    dsp *hpetal;         /**<Mapping from petal p/t/t to phase at aperture.*/
}aper_t;

/**
   contains the data associated with a LLT uplink path.
 */
typedef struct llt_t{
    loc_t *loc;          /**<The grid that defines the LLT pupil*/
    pts_t *pts;          /**<The LLT lower left grid point location*/
    dmat *amp;          /**<The amplitude defined on loc*/
    dcell *mcc;          /**<modal cross coupling matrix*/
    dcell *imcc;         /**<inverse of imcc*/
    dcell *ncpa;         /**<The LLT surface error*/
}llt_t;
/**
   contains the intensity statistics assiciated with a certain powfs for
physical optics wfs. */
typedef struct intstat_t{
    cccell *fotf;       /**<The final optf before fft and multiply with si to get i0. Used for MAP tracking and wfslinearity.*/
    dccell *sepsf;      /**<short expsoure PSF.*/
    dcell *i0;          /**<short exposure image. nsa x nllt*/
    dcell *gx;          /**<gradient of i0 along x*/
    dcell *gy;          /**<gradient of i0 along y*/
    dmat *i0sum;        /**<sum of i0 for each subaperture*/
    dmat *i0sumsum;     /**<sum of i0sum for all subapertures*/
    dcell *mtche;       /**<mtched filter operator along x/y, even if radpix=1*/
	dcell *cogmask;     /**<mask for cog calculation*/
}intstat_t;

/**
   contains the data associated with a certain type of WFS. not
necessarily physical optics WFS.x */
typedef struct powfs_t{
    /*Parameters about subapertures. */
    loc_t *saloc;       /**<lower left corner of the subaperture*/
    pts_t *pts;         /**<records lower left-most point of each sa in a regular grid.*/
	loc_t *loc;         /**<concatenated points for all subapertures.*/
	loccell *loc_dm;    /**<distorted loc mapped onto DM. size: (nwfs, ndm)*/
	loccell *loc_tel;   /**<distorted loc mapped onto pupil. size: (nwfs, 1). used for mkamp and ray trace. do not use in mkg*/

	dcell *amp;     /**<The real (after misregisteration/distortion) amplitude map defined on loc. nx is powfs.nwfs or 1.*/
	dcell *saa;     /**<The real (after misregisteration/distortion) subaperture area. same dimension as amp*/
	dmat *sumamp;       /**<sum of amp*/
	dmat *sumamp2;      /**<sum of amp.^2*/
    dmat* saasum;       /**<sum of saa*/
	dmat* saamax;		/**<Max of subaperture area across all wfs in this powfs. */
	dmat *saamin;       /**<Min of subaperture area across all wfs in this powfs. */
    real areascale;   /**<1./max(area noramlized by dsa*dsa)*/
    /*NCPA */
    dcell *opdadd;      /**<opdadd includes both common and NCPA OPDS. It is used for ray tracing*/
    dcell *opdbias;     /**<opdbias is used to compute gradient offset. It includes contributions only from NCPA OPD.*/
    dcell *gradoff;    /**<Offset to grads due to ncpa, sodium profile, etc. Copied to simu->gradoff for runtime update*/
    /*Physical optics */
    dtf_t *dtf;         /**<array of dtf for each wvl*/
    /*LGS Physical Optics */
    dcell *sodium;      /**<Loaded and downsampled sodium profile from the file for simulation.*/
    dcell* sodiumprep;  /**<Loaded and downsampled sodium profile from the file for i0 if not NULL.*/
    etf_t *etfprep;     /**<ETF for computing short exposure matched filter.*/
    etf_t *etfsim;      /**<ETF for simulation.*/
    etf_t *etfsim2;     /**<Second ETF for interpolation during simulation.*/
    dmat *focus;        /**<additional focus error. (llt->fnrange)*/
    llt_t *llt;         /**<uplink aperture parameters*/
    dccell *saimcc;     /**<inverse of p/t/t model cross coupling matrix for each subaps to compute ztilt.*/
    dcell *srot;        /**<subaperture rotation wrt LLT*/
    dcell *srsa;        /**<subaperture distance wrt LLT*/
    dmat *srsamax;      /**<max of srsa for each llt.*/
    dmat *pixoffx;      /**<Actual pixel offset wrt lenslet. Along r/a if radpix is true.*/
    dmat *pixoffy;      /**<Actual pixel offset wrt lenslet*/
    /*Geometric gradient */
    dspcell *GS0;        /**<gtilt (average gradient) on ampm*/
    dcell *neasim;      /**<LL' decomposition of nea covariance matrix, with
			   unit in radian, at dtrat, to be used in simulation
			   for geometric wfs model.*/
    /*Matched filter */
    lcell *sprint;      /**<indices of subapertures to print sanea*/
    intstat_t *intstat; /**<matched filter i0 and its derivative.*/
    dmat *dtheta;       /**<sampling of the imaging fft grid. wvl/(embfac*dxsa);*/
    dcell *bkgrnd;      /**<wfs background image. from parms->powfs[ipowfs].bkgrndfn.*/
    dcell *bkgrndc;     /**<wfs background image calibration. from parms->powfs[ipowfs].bkgrndfnc.*/
    //subaperture noise equivalent angles.
    dcell *sanea;       /**<computed sanea (nsa*3) for each WFS. (rad^2)*/
    //int namp;         /**<number of amplitude maps*/
    int pixpsax;        /**<number of detector pixels along x*/
    int pixpsay;        /**<number of detector pixels along y*/
    int notfx;          /**<PSF is extended to this side before computing OTF*/
    int notfy;          /**<PSF is extended to this side before computing OTF*/

    locfft_t **fieldstop;/**<For computing field stop (aka focal plane mask, spatial filter)*/
    struct pywfs_t *pywfs;/**<For pyramid WFS*/
    struct fit_t *fit;  /**<Fit turbulence to lenslet grid. For aliasing computation.*/
}powfs_t;

/**
   NGS mode and reconstructors in ad hoc split tomography. Only useful with 1 or 2 DMs.
 */
typedef struct ngsmod_t{
    real hs;      /**<height of LGS*/
    real hdm;      /**<height of upper DM.*/
    real scale;   /**<(1-ht/hs)^-2*/
    real aper_fcp;/**<piston term in focus in plocs.*/
    real lp2;     /**<LPF coefficient for P(Rngs,1)*/
    dcell *MCCP;    /**<cross coupling of the NGS modes for each direction. Hm'*W*Hm*/
    dmat *MCC;      /**<cross coupling of the NGS modes. 2x2 for 1 dm. 5x5 for 2 dms*/
    dmat *IMCC_TT;  /**<inv of cross coupling of tip/tilt modes only.*/
    dmat *IMCC;     /**<inv of MCC.*/
    dmat *IMCC_F;   /**<inv of MCC only including focus mode*/
    dmat *MCCu;     /**<MCCu'*MCCu=MCC*/
    dcell *GM;      /**<ngsmod vector to gradient operator*/
    dccell *Rngs;    /**<NGS reconstructor from NGS grad to NGS mod vec. pinv of GM*/
    dcell *Pngs;    /**<Reconstruct DM command to NGS modes with weighting specified by tomo.ahst_wt*/
    dcell *Modes;   /**<DM vector for the modes*/
	dcell *Pngs2;   /**<Reconstruct DM command to Modes2*/
	dcell *Modes2;   /**<DM vector for the modes (only for Pngs)*/
    dspcell *Wa;    /**<Aperture weighting. Ha'*W*Ha. It has zeros in diagonal. Add tikholnov*/
    lmat *modvalid; /**<Flag of valid modes that has multi-rate control*/
    int nmod;       /**<nmod: 5 for 2 dm, 2 for 1 dm.*/
    int ahstfocus;  /**<records parms->tomo.ahst_focus*/
    int indfocus;  /**<Include focus in NGS controlled modes. Records the index*/
    int indps;     /**<Include plate scale in NGS controlled modes. Records the index*/
    int indastig;  /**<Include astigmatism mode in NGS controlled modes. Records the index*/
}ngsmod_t;
/**
   contains data for Fourier Domain Preconditioner.
 */
typedef struct fdpcg_t{
    loccell *xloc;  /**<record recon->xloc*/
    csp *Minv;     /**<inverse of fourier domain preconditioner matrix M.*/
    ccell *Mbinv;  /**<block version of Minv. (in permuted order)*/
    lmat *perm;    /**<Permutation vector to get block diagonal matrix*/
    lmat *permhf;  /**<Permutation vector to be used when complex2real fft is used. Size is about half of perm.*/
    long nxtot;    /**<Total number of reconstructed phase points*/
    int square;    /**<Whether xloc is square*/
    int scale;     /**<Do we need to scale after fft.*/
    long nbx;      /**<Basic frequency range in x*/
    long nby;      /**<Basic frequency range in y nb=nbx*nby.*/
    long bs;       /**<Size of each diagonal block.*/
    /*Temperary Data used during execution*/
    /*Cannot save xhat here because we may be calling fdpcg from multiple
     * threads simultanously during assembly of control matrix*/
    //ccell *xhat, *xhat2, /**<FFT of each xout and multiplied. */
}fdpcg_t;

/**
   contains MOAO related data
*/
typedef struct moao_t{
    int used;         /**<Whether this MOAO is used or not*/
    loccell *aloc;      /**<Actuator grid*/
    mapcell *amap;      /**<Points to aloc->map*/
    dspcell *HA;       /**<Propagator from this aloc to PLOC*/
    dcell *NW;        /**<null modes and constraints*/
    dmat *W1;         /**<Weighting matrix on PLOC. same as recon->W1*/
    dsp *W0;          /**<Weighting matrix on PLOC. same as recon->W0*/
    dcell *actcpl;    /**<actuator coupling factor. 0 means actuator is outside
			 of FoV and need to be slaved.*/
    dspcell *actslave; /**<Slaving operator for actuators not illuminated*/
    dmat *aimcc;      /**<used for tip/tilt removal from DM commands.*/
    lcell *actstuck;  /**<stuck actuators*/
    lcell *actfloat;  /**<floating actuators*/
}moao_t;
/**
   A convenient wrap of the data to embed into muv_t for applying invpsd to opds defined on xloc.
*/
typedef struct invpsd_t{
    dcell *invpsd;    /**<inverse of the turbulence PSF*/
    ccell *fftxopd;   /**<temporary array to apply inverse PSD in Fourier domain.*/
    loccell *xloc;     /**<points to recon->xloc*/
    int     square;    /**<whether opd is on square xloc or not.*/
}invpsd_t;
/**
   A convenient wrap of the data to embed into muv_t for applying fractal
   regularization to opds defined on xloc.
 */
typedef struct fractal_t{
    dcell *xopd;      /**<A square array to embed x on xloc into for fractal */
    loccell *xloc;      /**<points to recon->xloc*/
    const real *wt;  /**<weight of each layer*/
    real  r0;        /**<The Fried parameter*/
    real  L0;        /**<The outer scale*/
    real  scale;     /**<An additional scaling factor*/
    long   ninit;      /**<The initial size to do with covariance matrix. 2 is minimum*/
}fractal_t;
/**
   Holding parameters for DM fitting.
 */
typedef struct fit_t{
    //Input data from recon and parms. Do not free.
    loccell *xloc;     /**<Input grid for DM fitting*/
    loc_t   *floc;     /**<intermediate pupil plane grid. */
    loccell *aloc;     /**<Destination grid (DM actuator)*/
	dcell *amod;	   /**<recon->amod.*/
    dsp  *W0;          /**<floc weighting for circle of diam aper.d*/
    dmat *W1;          /**<floc weighting for circle of diam aper.d*/
    dmat *thetax;      /**<DM fitting directions*/
    dmat *thetay;      /**<DM fitting directions*/
    dmat *wt;          /**<DM fitting weights*/
    dmat *hs;          /**<DM fitting GS height*/
    lcell *actstuck;
    lcell *actfloat;
    //Flags
    char **misreg;
    fit_cfg_t flag;
    int notrecon;      /**<Not for reconstruction*/
    int isref;         /**<Do not free generated data if isref is set*/
	int modal;		   /**<Modal reconstruction*/
    //Generated data
    dspcell *HXF;      /**<ray tracing propagator from xloc to floc for fitting directions.*/
    cell *HA;       /**<ray tracing from aloc or amod to floc for fitting directions.*/
    dcell *actcpl;
    dspcell* actextrap; /**<actuator interpolation*/
    dspcell *actslave;  /**<force slave actuators to have similar value to active neighbor ones.*/
    dcell *NW;         /**<null modes for DM fit.*/

    muv_t FR;          /**<DM fit right hand size matrix, solve FL*x=FR*y*/
    muv_t FL;          /**<DM fit left hand size matrix*/

}fit_t;

/**
   contains data related to wavefront reconstruction and DM fitting. */
typedef struct recon_t{
    real r0;         /**<r0 used in reconstruction. may get updated in cn2 estimation*/
    real L0;         /**<L0 used in reconstruction. may get updated in cn2 estimation*/
    dmat *ht;          /**<height of the layers to do tomography.*/
    dmat *wt;          /**<weight of the layers to to tomography. may get updated in cn2 estimation*/
    dmat *os;          /**<over sampling of the layers.*/
    dmat *dx;          /**<sampling in meter of the layers*/
	wfsr_cfg_t *wfsr;  /**<references parms->wfsr */
    loccell *saloc;    /**<referenced from powfs.saloc*/
    loc_t *ploc;       /**<Grid on pupil for tomography*/
    map_t *pmap;       /**<square grid of ploc.*/

    loccell *xloc;      /**<reconstructed atmosphere grid.*/
    mapcell *xmap;      /**<The map of xloc (only if tomo.square is true)*/
    mapcell *xcmap;     /**<The map of xloc on non-cone coordinate, with floc sampling.*/
    dcell *xmcc;       /**<used for tip/tilt removal from tomographic screens.*/
    lmat *xnx;         /**<A convenient reference of the size of each xloc*/
    lmat *xny;         /**<A convenient reference of the size of each xloc*/
    lmat *xnloc;       /**<A convenient reference of the size of each xloc*/
    map_t *fmap;       /**<Grid on pupil for DM fitting*/
    loc_t *floc;       /**<Grid on pupil for DM fitting. */

    loccell *aloc;      /**<actuator grid*/
    mapcell *amap;      /**<square grid of actuators*/
    mapcell *acmap;     /**For caching DM to intermediate plane during fitting (GPU)*/
    lmat *anx;        /**<Size of each amap*/
    lmat *any;        /**<Size of each amap*/
    lmat *anloc;      /**<Size of each aloc*/
    lmat *ngrad;      /**<Size of each grad for each wfs*/
    lcell *actfloat;   /**<floating actuators*/
    lcell *actstuck;   /**<stuck actuators*/
    dcell *amod;       /**<ndmx1. Zernike/KL modes defined on aloc for modal control*/
	dcell *amodpinv;   /**<Pinv of amod */
    lmat *anmod;       /**<Sizeof of amod*/

    fit_t *fit;        /**<Holding data and parameters for DM fitting.*/

    dcell *aimcc;      /**<used for tip/tilt removal from DM commands.*/
    dsp *W0;          /**<floc weighting for circle of diam aper.d*/
    dmat *W1;          /**<floc weighting for circle of diam aper.d*/
    //dmat *fitwt;       /**<fit weighting in each direction.*/
    dspcell *L2;        /**<Laplacian square regularization.*/
    dspcell *L2save;    /**<save old L2 to update the tomography matrix.*/
    invpsd_t *invpsd;  /**<data to apply inverse of psf to opd on xloc*/
    fractal_t *fractal;/**<data to apply fractal regularization on opd on xloc*/
    fdpcg_t *fdpcg;    /**<fdpcg preconditioner data.*/
    dspcell *GP;        /**<Gradient operator from HXW. GX=GP*H for each wfs..*/
    dspcell *HXW;       /**<Ray tracing operator from xloc to ploc for all WFS.*/
    dspcell *HXWtomo;   /**<Like GXtomo*/
    dspcell *GX;        /**<Gradient operator for all WFS from each layer of xloc*/
    dspcell *GXtomo;    /**<GX for tomography right hand side. excluding NGS in split tomography*/
	dspcell *GXhi;      /**<GX for high order WFS.*/
    dspcell *GXlo;      /**<GX for low order WFs*/
    dcell *GXL;        /**<dense GX for low order WFS in MV split tomography.*/
    dcell *MVRngs;     /**<NGS recon for MV split tomography*/
    dcell *MVModes;    /**<MVST Modes (svd'ed)*/
    dcell *MVGM;       /**<NGS WFS gradient operator from MVST Modes.*/
    dcell *MVFM;       /**<NGS Focus reconstructed from MVST Modes.*/
    cell *GA;        	/**<actuator or modes to wfs grad.*/
    cell *GAlo;         /**<GA or GM of low order WFS.*/
    cell *GAhi;      	/**<GA of high order WFS.*/
    //dcell *GM;          /**<GM for all WFS.*/
    //dcell *GMhi;        /**<GM for high order WFS.*/
    dspcell *HA_ncpa;   /**<ray tracing from aloc to floc for NCPA directions*/
    dcell *TT;         /**<TT modes for LGS WFS*/
    dcell *PTT;        /**<pinv of TT for tt removal from LGS gradients*/
    dcell *FF;         /**<Global focus or Differential focus modes for LGS wfs low rank*/
    dcell *PFF;        /**<pinv of FF. to use in RTC.*/
    dcell *TTF;        /**<Concatenation of TT and DF*/
    dcell *PTTF;       /**<pinv of TTF*/
	dcell *TTF_LHS;    /**<Concatenation of TT and DF for tomography left hand side*/
	dcell *PTTF_LHS;   /**<pinv of TTF for tomography left hand side*/
    dspcell *ZZT;       /**<single point piston constraint in tomography.*/
    dcell *DMTT;       /**<DM tip/tilt mode.*/
    dcell *DMPTT;      /**<DM tip/tilt reconstructor.*/
    dcell *actcpl;     /**<actuator coupling factor. 0 means actuator is outside of FoV and need to be slaved.*/
    dspcell *actextrap; /**<Interpolation operator for floating actuators and edge actuators. Slaving does not work well in CG. */
    dspcell *sanea;     /**<Measurement noise covairance, sanea^2 for each wfs in radian^2*/
    dspcell *saneal;    /**<cholesky decomposition L of sanea^2 for each wfs to compute noise propagation*/
    dspcell *saneai;    /**<inverse of sanea^2 in radian^-2 for each wfs*/
    dcell *ecnn;       /**<covairance of Hx*(E*Cnn*E^t)*Hx^t: noise propagation to science.*/
    dmat *neam;        /**<subaperture averaged nea for each wfs*/
    real neamhi;     /**<average of neam for high order wfs.*/
    real sigmanlo;   /**<Wavefront error (m^2) due to noise for lo order.*/
    real sigmanhi;   /**<Wavefront error (m^2) due to noise for high order. Needs to multiply with simu->gradscale^2.*/

    muv_t RR;          /**<tomography right hand side matrix, solve RL*x=RR*y*/
    muv_t RL;          /**<tomography left hand side matrix*/
    muv_t LR;          /**<least square reconstructor rhs*/
    muv_t LL;          /**<least square reconstructor lhs. solve LL*x=LR*y*/
    dmat *MVM;        /**<Matrix vector multiply*/
    dcell *MVA;        /**<Correction to MVM*g by (MVA-I)*a for PSOL.*/
    moao_t *moao;      /**<for MOAO DM fitting*/
    /*For focus tracking. */
    dcell *GFlgs;      /**<Focus to LGS gradients*/
    dcell *GFngs;      /**<Focus to NGS gradients*/
    dcell *GFall;      /**<Focus to WFS Gradients.*/
    dcell *RFlgsg;     /**<focus reconstruction for each LGS from grad*/
    dcell *RFngsg;     /**<focus reconstruction for all TTF NGS from grad.*/
    /*For focus offloading*/
    dcell *RFdm;       /**<Focus from DM commands. For telescope offloading*/

    dcell *GRall;      /**<Truth zernike modes to gradient. Was only radial modes but now all modes (depends on recon.twfs_radonly*/
    dcell* GRtwfs;     /**<Truth zernike modes to gradient for twfs*/
    dcell* RRtwfs;     /**<Truth zernike modes reconstruction from twfs grads*/
    dcell* GRlgs;      /**<Truth zernike modes to gradient for LGS (sodium fit)*/
    dcell* RRlgs;      /**<Truth zernike modes reconstruction from LGS grad adjustments (sodium fit)*/

	//for petal mode mitigation using phase retrieval
	petal_t **petal;   /**<Petaling mode reconstruction tools.*/
	dsp *apetal;       /**<Petal mode defined at the ground DM*/

    //For common path dithering
    dcell *dither_m;   /**<The dither mode added to DM command (ndm*1)*/
    int dither_npoint; /**<The dither period*/
    int dither_dtrat;  /**<The dtrat of the powers that requests dithering*/
    int dither_md;     /**<multi-mode dithering bin size in amod*/
    dcell *dither_rg;  /**<The dither mode recon from grads (nwfs*nwfs)*/
    dcell *dither_ra;  /**<The dither mode recon from dm commands (ndm*ndm)*/
    ngsmod_t *ngsmod;  /**<ngs mod in ad hoc split tomography.*/
    cn2est_t *cn2est;  /**<For Cn2 Estimation*/
    dcell *dm_ncpa;    /**<NCPA calibration for DM. add to dmreal.*/

    int lowfs_gtilt;   /**<=1 if any low order wfs use gtilt in recon/simu*/
    int npsr;          /**<number of reconstructor phase screens.*/
    int nthread;       /**<number of threads in reconstruction.*/
    int cxxalg;        /**<records parms->tomo.cxxalg*/
	int nwfsr;		   /**<references parms->nwfsr */
    //For Error PSD computation
    cell *Herr;      /**<Ray tracing from DM along science directions for a few points. dcell for modal control. sparse for zonal control*/
}recon_t;

typedef struct sim_save_t{
    /*zfarrs to save telemetry data.*/
    zfarr **wfspsfout; /**<special file to save wfs psf history*/
    zfarr **ztiltout;  /**<special file to save zernike wfs tilt history*/
    /*Evaluation directions PSF. */
    zfarr * evlpsfolmean;  /**<science field psf OL time average*/
    zfarr **evlpsfmean;    /**<science field psf CL time average*/
    zfarr **evlpsfhist;    /**<to save time history of science field psf*/
    zfarr **evlopdcov;     /**<science field OPD covariance*/
    zfarr **evlopdmean;    /**<science field OPD mean*/
    zfarr *evlopdcovol;    /**<science field OPD covariance (open loop)*/
    zfarr *evlopdmeanol;   /**<science field OPD mean (open loop)*/
    zfarr **evlpsfmean_ngsr;    /**<science field psf CL time average with NGS mode removed*/
    zfarr **evlpsfhist_ngsr;    /**<to save time history of science field psf with NGS mode removed*/
    zfarr **evlopdcov_ngsr;     /**<science field OPD covariance with NGS mode removed*/
    zfarr **evlopdmean_ngsr;    /**<science field OPD mean with NGS mode removed.*/
    zfarr **ecovxx;     	/**<the time history of xx used to calculate ecov.*/
    /*Deformable mirror. */
    zfarr *dmerr;
    zfarr *dmint;
    zfarr *dmrecon;
    zfarr *dmreal;
    dmat *ttmreal;
    zfarr *dmcmd;
    zfarr *dmproj;
    /*Low order modes */
    zfarr *Merr_lo;
    zfarr *Mint_lo;
    zfarr *opdr;
    zfarr *opdx;
    /*science */
    zfarr **evlopdcl;
    zfarr **evlopdol;
    zfarr **wfsopd;
    zfarr **wfsopdol;
    zfarr **wfslltopd;
    /*gradients */
    zfarr **gradcl;
    zfarr **gradnf;
    zfarr **gradgeom;
    zfarr **gradol;
    zfarr **intsny;
    zfarr **intsnf;
    zfarr **dm_evl;
    zfarr** dm_wfs;
    /*psds*/
    zfarr* psdcl;
    zfarr* psdol;
    zfarr* psdcl_lo;
    zfarr* psdol_lo;
    //
    zfarr* restwfs;    /**<Truth wfs output*/
    dcell* fsmerrs;    /**< file to store fsmerr history*/
    dcell* fsmcmds;    /**< file to store fsmcmd history*/
    dcell* ltpm_real; /**<file to store ltpm_real history*/
	dcell* gain; /**<gain update time history. */
}sim_save_t;
/*
  data wrap for wfsints.
*/
typedef struct wfsints_t{
    dcell *ints;
    ccell *psfout;
    dcell *pistatout;
    const dmat *gradref;
    const dmat *opd;
    const dmat *lltopd;
    int iwfs;
    int isim;
}wfsints_t;
/**
  data for dithering statistics collection
*/
typedef struct dither_t{
    //For PLL
    real delta; /**<PLL estimation of servo lag (only) at every time step*/
    real deltam;/**<Output of PLL*/
    real deltao;/**<Offset of delta from outer loop*/
    real delay; /**<Diference of delay from 2 frame due to beam propagation*/
    real a2m;    /**<input dither amplitude*/
    real a2me;   /**<measured dither amplitude*/
    dmat *a2mv;  /**<input dither amplitude per mode*/
    dmat *a2mev;  /**<measured dither amplitude per mode*/
    //For matched filter
    dcell *imb;   /**<accumulated im*/
    dcell *imx;   /**<accumulated cos()*im */
    dcell *imy;   /**<accumulated sin()*im */
    dcell *i0;    /**<accumulated imb for matched filter*/
    dcell *gx;    /**<accumulated imx for matched filter*/
    dcell *gy;    /**<accumulated imy for matched filter*/
    //For CoG
    dmat *ggm;    /**<Accumulated cos()*tt_x and sin()*tt_y*/
    dmat *gg0;    /**<Averaged ggm*/
    //For DM dithering
    dcell *mr;    /**<Mode reconstructed from actuator and gradients*/
}dither_t;
/**
 * A collection of flags during simulation for each powfs.
 * */
typedef struct {
    int do_phy;   /**<Do physical optics*/
    int do_pistat;/**<Collect pixel intensity statistics*/
    int gradout;  /**<Gradient count or 0 if no output*/
    int pllout;   /**<PLL output count or 0 if no update*/
    //int ogcount;   /**<Number of optical gain accumulations*/
    int ogacc;    /**<OG accumulation count or 0 if no action*/
    int ogout;    /**<OG output count or 0 if no update*/
	//int zoomout;  /**<Trombone zoom output*/
}wfsflags_t;
/**
   contains all the run time data struct.
*/
typedef struct sim_t{
    /*A few data structs generated in beginning of simulation*/
    /*random stream. We maintain separate random streams for each purpose,
       derived from the same seed, so that multi-threading will produce same
       result */
    rand_t *wfs_rand;  /**<random stream for each wfs.*/
    rand_t *atm_rand;  /**<random stream for atmosphere turbulence generation*/
    rand_t *atmwd_rand;/**<random stream for wind direction*/
    rand_t *telws_rand;/**<random stream for wind shake*/
    rand_t *init_rand; /**<random stream to initialize other streams*/
    rand_t *misc_rand; /**<For misc purposes*/
    /*Atmosphere*/
    genatm_t *atmcfg;
    mapcell *atm;       /**<fine sampled simulation turbulence screens*/
    dmat *atmscale;     /**<Scale atmosphere during ray tracing to simulate r0 variation.*/
    mapcell *cachedm;   /**<grid cache dm actuator to a finer sampled screen. for
			   fast ray tracing to WFS and aper*/
    dmat *winddir;     /**<input wind direction*/

    /*Optional surface errors in M1, M2, or M3*/
    rmapcell *tsurf; /**<input tilted M3 surface read from parms->tsurf*/
    mapcell *surf;      /**<input surface: M1, M2 or else. common to all wfs and science field.*/

    /*Telescope windshake time series and direction.*/
    dmat *telws;      /**<Telescope wind shake time series, along ground layer wind direction*/

    /*WFS data for each time step. Each has a cell for each wfs*/
    dcell *wfsopd;     /**<WFS Ray tracing result*/
    dccell *ints;      /**<WFS subaperture images.*/
    dccell *intsout;   /**<Time averaged subaperture image*/
    cccell *wfspsfout; /**<output WFS PSF history.*/
    dccell *pistatout; /**<WFS time averaged tip/tilt removed PSF*/
    dcell *gradcl;     /**<cl grad output at step isim.*/
    dcell *gradacc;    /**<accumulate gradident for dtrat>1*/
    dcell *gradlastcl; /**<cl grad from last time step, for reconstructor*/
    dcell *gradlastol; /**<psol grad from last time step, for reconstructor*/
    dcell *cn2res;     /**<Cn2 Estimation Result*/
    dcell *gradoff;    /**<Offset to grads to subtract from measurement. */
    dcell *gradoffacc; /**<Accumulates gradoff to determine its average*/
    dcell *gradoffdrift;/**<for drift control*/
    int gradoffnacc; /**<gradoffacc counter*/
    int gradoffisim; /**<last isim the new gradoff is activated*/
    int gradoffisim0; /**<last isim the gradoffacc is reset*/
    real eptwfs;     /**<Twfs reference vector servo gain.*/
    /*CoG gain adjustment*/
    dcell *gradscale;  /**<Gain adjustment for cog and pywfs.*/
    dcell *gradscale2; /**<Gain scaling for other dithering modes.*/
    /*LLT*/
    dcell *llt_ws;     /**<LLT uplink jitter (per LLT)*/
    dcell *ltpm_lpf;  /**<For center launch only: LLT common path pointing mirror LPF*/
    dcell *ltpm_cmd;  /**<For center launch only: LLT common path pointing mirror command*/
    dcell *ltpm_real; /**<For center launch only: LLT common path pointing mirror state after sho filtering*/
    sho_t** ltpm_sho; /**<For center launch only: LLT common path pointing mirror state*/
	/*petaling mode*/
	dccell *petal_i0;   /**<WFS accumulated i0.*/
	dcell *petal_m;   /**<Petal estimator output.*/
    /*Tomography*/
    dcell *opdr;       /**<reconstructed OPD defined on xloc in tomography output.*/
    dcell *gngsmvst;   /**<opdr to NGS gradient.*/
    dcell *cgres;      /**<CG residuals for tomography and fit*/
    /*Only on split tomography*/
    dmat *clem;        /**<lgs/ngs mod error in split tomography*/
    dcell *clemp;      /**<lgs/ngs mod error per direction. only on-axis is computed.*/
    dmat *corrNGSm;    /**<Correction of NGS mode. (integrator output)*/
    dmat *cleNGSm;     /**<Close loop ngs mods in split tomogrpahy. */
    dmat *oleNGSm;     /**<Open loop ngs mods in split tomogrpahy. */
    dcell *cleNGSmp;   /**<(M'*w*phi);*/
    dcell *oleNGSmp;   /**<(M'*w*phi); for OL*/
    dcell *res;        /**<warping of ole,cletomo,cle,clem for easy saving.*/
    //dmat *timing;      /**<Timing and memory using for each step*/

    dcell *resdither;   /**<Phase and amplitude estimation of dithering*/

    /*DM commands.*/
    dcell *dmpsol;     /**<DM command for PSOL feedback*/
    dcell *dmtmp;      /**<Holds a temporary dm vector. Maybe in zonal or modal space.*/
    dcell *dmtmp2;     /**<Holds a temporary dm vector. Always in zonal space.*/
    dcell *dmcmd;      /**<This is the final DM command send to DME.*/
    dcell *dmreal;     /**<This is the actual position of DM actuators after
			  receiving command dmcmd. Should only be used in
			  system, not in reconstruction since it is unknown.*/
    dmat *ttmreal;     /**<TT mirror command*/
    mapcell *dmrealsq;  /**<dmreal embeded into an square map, zero padded.*/
    dcell *dmproj;     /**<The projection of atm onto DM space directly.*/
    mapcell *dmprojsq;  /**<dmproj embeded into square map, zero padded.*/
    dccell *wfspsol;    /**<time averaged dm command (dtrat>1) for psol grad*/
    dcell *dmhist;     /**<histogram of dm commands. if dbg.dmhist is 1.*/
    hyst_t**hyst;      /**<Hysterisis computation stat*/
    dcell *dmadd;      /**<cell of dm vector to simulate turbulence (added to
			  integrator output).  A cell for each DM. In each cell,
			  one colume is for each time step. Wraps over in the
			  end*/
    /*High order*/
    servo_t *dmint;    /**<dm integrator. (used of fuseint==1)*/
    dcell *dmrecon;      /**<direct high order fit output*/
    dcell *dmerr;      /**<equals to dmerr_store when there is output.*/
    dcell *dmerr_store; /**<high order dm error signal.*/
    /*Low order*/
    dcell *Merr_lo,*Merr_lo_store;    /**<split tomography NGS mode error signal.*/
    dcell *Merr_lo2;   /**<Saves LPF of Merr_lo result*/
    servo_t *Mint_lo;  /**<intermediate results for type II/lead filter*/
    dcell *Mngs;       /**<Temporary: NGS mode in DM commands*/
    /*llt pointing loop*/
    dcell *fsmerr,*fsmerr_store;     /**<uplink error*/
    dcell *fsmerr_drift;/**<Drift control of uplink*/

    servo_t **fsmint;    /**<uplink integrator output.*/
    sho_t **fsmsho;      /**<FSM sho response*/
    dcell *fsmcmd;      /**<FSM command*/
    dcell *fsmreal;     /**<FSM real position*/

    /*focus tracking loop*/
    dcell *LGSfocus;  /**<LGS focus error*/
    dcell *LGSfocus_drift;  /**<LGS focus drift error*/
    dcell *LGSfocusts; /**<Time history of focus error*/
    dmat *lgsfocuslpf;/**<low pass filtered individual LGS focus*/
    real ngsfocuslpf;/**<low pass filtered NGS focus*/
    //dmat *zoomerr;    /**<Trombone error signal from zoomavg*/
    dmat *zoomdrift; /**<Trombone error signal from i0/ib drift control*/
    lmat *zoomdrift_count;
    dmat *zoomint;    /**<Trombone integrator*/
    dmat *zoomavg;    /**<Trombone averager from gradients*/
    lmat *zoomavg_count;/**<Count of zoomavg accumulation*/
    dcell *zoompos;    /**<Trombone position history. for saving*/
    long zoompos_icol;  /**<Current column*/
    /*focus offloading*/
    dcell *telfocusint; /**<Telescope focus integrated*/
    dcell *telfocusreal;/**<Telescope focus integrated*/
    /*PSD computation*/
    dcell *dmerrts;    /**<Herr*dmerr time history*/
    dmat *Merrts;      /**<Time history of low order mode.*/
    /*science evaluation*/
    dcell *evlopd;     /**<Save science ifeld opd for use in perfevl_mean().*/
    dcell *evlopdmean;    /**<science field opd mean*/
    dmat  *evlopdground;  /**<evaluation opd for ground layer turbulence to save ray tracing.*/
    dcell *evlpsfmean;    /**<science field psf time average*/
    dcell *evlpsfolmean;  /**<science field OL PSF time averging*/
    dcell *evlopdcov;     /**<science field opd covariance*/

    dmat *evlopdcovol;   /**<science field opd covariance (open loop)*/
    dmat *evlopdmeanol;  /**<science field opd mean (open loop)*/
    dcell *evlpsfmean_ngsr;    /**<science field psf time average with NGS mode removed.*/
    dcell *evlopdcov_ngsr;     /**<science field opd covariance with NGS mode removed.*/
    dcell *evlopdmean_ngsr;    /**<science field opd mean with NGS mode removed.*/

    /*Optinal telemetry saving for PSF reconstruction.*/
    dcell *ecov;       /**<covariance of Hx*x-Ha*a for science directions.*/
    dcell *gcov;       /**<covariance of psuedo open loop gradients.*/

    /*save performance results to file using mmap*/
    dcell* resp;       /**<Stores olmp, clmp, olep, clep*/
    dcell* olmp;       /**<OL mode coefficient per direction.*/
    dcell* clmp;       /**<CL mode coefficient per direction.*/
    dcell* olep;       /**<OL error per direction.*/
    dcell* clep;       /**<CL error per direction.*/

    dmat* ole;         /**<field averaged OL error*/
    dmat *cle;         /**<field averaged CL error*/

    //Temporary
    dmat *ngsmodlpf;  /**<For removal low frequency component of ngsmod from LGS recon*/

    /*MOAO*/
    dcell *dm_wfs;   /**<moao DM command computed for wfs*/
    dcell *dm_evl;   /**<moao DM command computed for science field*/
    //Timing
	real tk_s0;		/**<first seed simulation start time.*/
	real tk_si;		/**<current seed simulation start time.*/
    real tk_istart;  /**<Start time of each isim*/
	real tk_iend;    /**<End time of each isim*/
	real tk_i1;		/**<2nd step start time*/
    real tk_eval;    /**<time spent in perfevl in this step*/
    real tk_recon;   /**<time spent in reconstruct in this step*/
    real tk_cache;   /**<time spent in cachedm in this step*/
    real tk_wfs;     /**<time spent in wfsgrad in this step*/

    /*A few data wraps for multi-threading*/
    propdata_t *cachedm_propdata; /**<wrapped data for ray tracing from aloc to cachedm*/
    propdata_t *wfs_propdata_dm;  /**<wrap of data for ray tracing from DM in wfsgrad.c*/
    propdata_t *wfs_propdata_atm; /**<wrap of data for ray tracing from ATM in wfsgrad.c*/
    propdata_t *evl_propdata_atm;
    propdata_t *evl_propdata_dm;
    thread_t **cachedm_prop;   /**<wrapped cachedm_propdata for threading*/
    thread_t **wfs_prop_dm;   /**<wrap of wfs_propdata_dm for threaded ray tracing*/
    thread_t **wfs_prop_atm;  /**<wrap of wfs_propdata_atm for threaded ray tracing*/
    thread_t **evl_prop_atm;
    thread_t **evl_prop_dm;

    wfsints_t *wfs_intsdata;  /**<wrap of data for wfsints.c*/
    thread_t **wfs_ints;     /**<wrap of wfs_intsdata for threaded processing*/

    thread_t *wfsgrad_pre;  /**to call wfsgrad_iwfs or gpu_wfsgrad_queue in threads.*/
    thread_t *wfsgrad_post; /**to call wfsgrad_post in threads.*/
    thread_t *perfevl_pre;  /**to call perfevl_ievl or gpu_perfevl_queue in threads.*/
    thread_t *perfevl_post; /**to call gpu_perfevl_sync in threads.*/

    sim_save_t *save;  /**<Telemetry output*/
    status_t *status;  /**<status report to scheduler.*/
    dither_t **dither;
    /**For testing*/
    ccell *opdrhat;    /**<For wind estimation (testing)*/
    ccell *opdrhatlast;/**<for wind estimation.(testing)*/
    /**Saved for plotting*/
    const char **plot_legs;///legend
    dccell *plot_res;    ///results array
    int plot_isim;      ///previous plotted isum;
    /*A few indicators*/
    int wfsints_isa;   /**<sa counter for wfsints*/
    int perfevl_iground;/**<index of the layer at ground*/
    int seed;          /**<current running seed.*/
    int iseed;         /**<index of current running seed.*/
    int wfsisim;       /**<record current simulations step for wfs.*/
    int perfisim;      /**<record current simulations step for pefevl.*/
    int reconisim;     /**<The time step for the gradlast data struct. =isim for OL, =isim-1 for CL*/
    wfsflags_t *wfsflags;/**<Runtime flags for each wfs*/
    /*maintain pointer of other structs for use in thread functions.*/
    const parms_t *parms; /**<pointer to parms*/
    const aper_t *aper;/**<pointer to aper*/
    recon_t *recon;    /**<pointer to recon*/
    powfs_t *powfs;    /**<pointer to powfs*/
    real last_report_time;/**<The time we lasted reported status to the scheduler.*/
    int tomo_update;   /**<Triggering setup_recon_tomo_upate*/
    int pause;         /**<pause simulation every this many steps. Copies from sim.pause*/
    //For synchronization. perfevl and wfsgrad waist for dmreal to be updated.
    //reconstruct waist for gradients to be available.
    int dmreal_isim;//which isim this dmreal is valid for
    unsigned int dmreal_count;//how many time this dmreal has been consumed (parms->evl.nevl + parms->nwfs)
    pthread_cond_t dmreal_condr;
    pthread_cond_t dmreal_condw;
    pthread_mutex_t dmreal_mutex;
    int wfsgrad_isim;//which isim this wfsgrad is from
    unsigned int wfsgrad_count;//how many times this wfsgrad has been consumed.
    pthread_cond_t wfsgrad_condr;
    pthread_cond_t wfsgrad_condw;
    pthread_mutex_t wfsgrad_mutex;

}sim_t;
#define CHECK_SAVE(start,end,now,every) ((now)>=(start) && (((every)>1 && ((now)+1-(start))%(every)==0) || (now)+1==(end)))

typedef struct global_t{
    const parms_t *parms;
    powfs_t *powfs;
    aper_t *aper;
    recon_t *recon;
    sim_t *simu;
    int iseed;
    int setupdone;
}global_t;

#endif

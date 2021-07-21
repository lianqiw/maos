/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#ifndef AOS_CUDA_CUWFS_H
#define AOS_CUDA_CUWFS_H


/**
   This file is work in progress and is not used.
 */


namespace cuda_wfs{
/*Basic geometric information for WFS.*/
struct wfscfg_t;
class cuwfs_info{
public:
    culoc_t *loc;   /**<location of OPD points*/
    curmat *amp;    /**<Amplitude map defined on loc*/  
    Array<int> *embed;
    int embednx, embedny;
    int iwfs;
    int igpu;
    int nsa;
    cuwfs_dbg(const parms_t *parms, const powfs_t *powfs, int _iwfs, int _igpu);
};
/*For field stop implementation*/
class cufieldstop_t{
    curmat *fieldmask;      /**<mask for field stop computation*/
    cucmat *wvf;
    cufftHandle plan;  /**<FFT plan*/
    Real *wvl;
    int nwvl;
    cudaStream_t stream;
public:
    cufieldstop_t(dmat *_mask, double *_wvl, int _nwvl, cudaStream_t _stream)
	:fieldmask(0), nwvl(_nwvl),stream(_stream){
	cp2gpu(&fieldmask, _mask);
	DO(cufftPlan2d(&plan, fieldmask->nx, fieldmask->ny, FFT_T_C2C));
	cufftSetStream(plan, stream);
	wvf=new cucmat(fieldmask->nx, fieldmask->ny);
	wvl=new Real[nwvl];
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    wvl[iwvl]=_wvl[iwvl];
	}
    }
    void apply(cuwfs_info *wfsinfo, curmat &opd, cudaStream_t stream);
    ~cufieldstop_t(){
	delete fieldmask;
	delete wvf;
	cufftDestroy(plan);
    }
};
class curand_t{
public:
    int nb;
    int nt;
    curandState *rstat;
    curand_t(int seed, int n);
};
/*Abstract class for WFS*/
class cuwfs_base{
protected:
    curand_t *rand;
    stream_t& stream;
    const int dtrat;
public:
    cuwfs_base(wfscfg_t *wfscfg);
    virtual void seeding(int seed);
    virtual void initsim();
    virtual void acc(curmat *opd);
    virtual void addnoise();
    virtual void output();
};
class cushwfs_t:public cuwfs_base{
protected:
    const int nsa;
    const Real dxsa;
    const int nxsa;
public:
    cushwfs_t(wfscfg_t *wfscfg);
    virtual void seeding(int seed){
	delete rand;
	rand=new curand_t(seed, nsa);
    }
};
class cushgeom_t:public cushwfs_t{
protected: 
    curmat *gradacc;
    curmat *gradcalc;
    curmat *nea;
public:
    cushgeom_t(wfscfg_t*);
    virtual void initsim(){
	gradacc->zero();
    }
    virtual void addnoise();
    virtual void output();
    virtual void calcg(curmat &opd, Real ratio);
    virtual void acc(curmat *opd);
};

//LLT configuration. Only useful in SHWFS Phy mode.
class cullt_t{
    cupts_t *pts;/**<LLT pupil sampling parameters*/
    culoc_t *loc;/**<location of OPD points for laser launch telescope*/
    curmat *amp;
    curmat *ncpa;/**<Aberration*/
    Real (**imcc)[3];/**<Ztilt*/
    cufftHandle plan_wvf, plan_otf;
public:
    cullt_t(wfscfg_t*);
};
/*SHWFS G tilt*/
class cushg_t:public cushgeom_t{
    cusp *GS0;
public:
    cushg_t(wfscfg_t* wfscfg);
    virtual void calcg(curmat* opd, Real ratio);
};
/*SHWFS Z tilt*/
class cushz_t:public cushgeom_t{
    Real (**imcc)[3];/**<For ztilt.*/
public:
    cushz_t(wfscfg_t* wfscfg);
    virtual void calcg(curmat* opd, Real ratio);
};

/*SHWFS Phy */
class cushphy_t:public cushwfs_t{
    cufftHandle plan1, plan2, plan3;
    int nsa;//total subapertures
    int msa;//each time handle this subapertures
    cuccell *dtf;
    cuccell *etf;
    int etfis1d;
    curmat *srot;//subaperture gradient vector rotation
    //noise
    Real bkgrnd, bkgrndc;
    curmat *bkgrnd2; curmat *bkgrnd2c;
    //gradient operator
    curmat *mtche;
    curmat *i0sum;
    //runtime data
    curcell *ints;
    curcell *pistatout;
public:
    cushphy_t(wfscfg_t*);
    ~cushphy_t(){
	delete dtf;
	delete etf;
	delete srot;
	delete bkgrnd2;
	delete bkgrnd2c;
	delete mtche;
	delete i0sum;
	delete ints;
	delete pistatout;
    }
    virtual void initsim(){
	ints->zero();
    }
    virtual void acc(curmat*);
    virtual void addnoise();
    virtual void output();
};
/*Pyramid WFS */
class cupywfs_t:public cuwfs_base{
    
};

/*SHWFS. If some components share the same value as previous WFS in the same GPU, share the data.*/
class cuwfs_t{
    cuwfs_info *wfsinfo;
    //pupil plane OPD grid
    cufieldstop_t *fieldstop; /**<Implementation of field stop*/
    //Geometric wfs
    cuwfs_base *geom;
    //Physical Optics WFS
    cuwfs_base *phy;
    //laser launch telescope
    cullt_t *llt;
    
    //Runtime data
    curmat *opd;
    curmat *opdadd;
    curmat *gradoff;
    stream_t stream;
public:
    cuwfs_t(const parms_t *parms, const powfs_t *powfs, int iwfs, int igpu);
    void initsim();
    void seeding(int seed);
};

};//namespace
#endif

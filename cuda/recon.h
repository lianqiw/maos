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
#ifndef AOS_CUDA_RECON_H
#define AOS_CUDA_RECON_H

#include "gpu.h"
#include "solve.h"
#include "recon_geom.h"
typedef struct cumoao_t{
    curcell *fitNW;
    cusp *actslave;
    curmat *cubic_cc; 
    cugrid_t amap;/**<Info of amap*/
    cugrid_t acmap;/**<Info of acmap*/
    cugrid_t fmap;/**<Info of fmap*/
    /*aperture weighting*/
    W01_T *W01;
    curcell *rhs;/*right hand side.*/
    float *pis; /*a temporary variable*/
    curmat *xp; /*temporary array*/
    curmat *xp2;/*temporary array*/
    curmat *tmpNW;/*temporary array*/
    cumoao_t(){
	memset(this, 0, sizeof(*this));
    }
}cumoao_t;

class curecon_t{
public:
    curecon_geom *grid;
    cusolve_r *FR;
    cusolve_l *FL;
    cusolve_r *RR;
    cusolve_l *RL;
    cusolve_l *MVM;
private:
    curcell *gradin; /**< The grad to operator on*/
   
    curcell *opdr;  /**<Reconstructed atm on xloc. Don't free to have warm restart. Free with new seed*/
    curcell *opdr_vec; /**<Referencing opdr in vector form*/
    curcell *dmfit;
    curcell *dmfit_vec;
    stream_t cgstream;
    curcell *opdr_save;/*for debugging*/
    curcell *dmfit_save;/*for debugging*/
    curcell *RFdfx;
    curcell *GXL;

    curcell *gngsmvst;
    curcell *deltafocus;
    
    friend void gpu_setup_moao(const PARMS_T *parms, RECON_T *recon);
    friend void gpu_moao_filter(SIM_T *simu);
    friend void gpu_moao_recon(SIM_T *simu);
    friend void gpu_moao_FitR(curcell **xout, float beta, SIM_T *simu, cumoao_t *cumoao, float thetax, float thetay, float hs, const float alpha);
    friend void gpu_moao_FitL(curcell **xout, float beta, const void *A, const curcell *xin, float alpha, stream_t &stream);
    
public:
    cumoao_t *moao;/**<moao configurations for GPU*/
    curcell **dm_wfs;/**<moao results for wfs for warm restart*/
    curcell **dm_evl;/**<moao results for evl for warm restart*/
    stream_t *moao_stream;
    CGTMP_T cgtmp_moaowfs;
    CGTMP_T cgtmp_moaoevl;
    /*the following data reside in the gpu memory*/
    friend void gpu_update_recon(const PARMS_T *parms, RECON_T *recon);
public:
    void init(const PARMS_T *parms, POWFS_T *powfs, RECON_T *recon);
    curecon_t(const PARMS_T *parms=0, POWFS_T *powfs=0, RECON_T *recon=0){ 
	grid=0;
	FR=0; FL=0; RR=0; RL=0; MVM=0;
	gradin=0;opdr=0;opdr_vec=0;opdr_save=0;
	dmfit=0; dmfit_vec=0;dmfit_save=0;
	RFdfx=0;GXL=0;gngsmvst=0; deltafocus=0;
	moao=0;
	dm_wfs=0;dm_evl=0;moao_stream=0;
	if(!parms) return;
	init(parms, powfs, recon);
    }
    ~curecon_t(){
	delete grid;
	if(FL!=dynamic_cast<cusolve_l*>(FR)) delete FL;
	delete FR;
	if(RL!=dynamic_cast<cusolve_l*>(RR)) delete RL;
	delete RR;
	delete gradin;

	delete opdr;
	delete opdr_vec; 
	delete opdr_save;

	delete dmfit;
	delete dmfit_vec; 
	delete dmfit_save;

	delete moao;
	delete dm_wfs;
	delete dm_evl;
	delete moao_stream;
	delete RFdfx;
	delete GXL;
    }
    void reset(const PARMS_T *parms);
    float tomo(dcell **_opdr, dcell **gngsmvst, dcell **dfocus, const dcell *_gradin);
    float fit(dcell **_dmfit, dcell *_opdr);
    void mvm(dcell **_dmerr, dcell *_gradin);
    void tomo_test(SIM_T *simu);
    void fit_test(SIM_T *simu);
};

__global__ void apply_W_do(float *restrict out, const float *restrict in, const int *W0f, 
			   float alpha, int nx, int n);



void gpu_TomoR(curcell **xout, float beta, const void *A, const curcell *grad, float alpha, stream_t &stream);
void gpu_TomoRt(curcell **gout,float beta, const void *A, const curcell *xin, float alpha, stream_t &stream);
void gpu_TomoL(curcell **xout, float beta, const void *A, const curcell *xin, float alpha, stream_t &stream);
void gpu_FitR (curcell **xout, float beta, const void *A, const curcell *xin, float alpha, stream_t &stream);
void gpu_FitRt(curcell **xout, float beta, const void *A, const curcell *xin, float alpha, stream_t &stream);
void gpu_FitL (curcell **xout, float beta, const void *A, const curcell *xin, float alpha, stream_t &stream);
void gpu_Tomo_fdprecond(curcell **xout, const void *A, const curcell *xin, stream_t &stream);

void cumuv(curcell **out, float beta, cumuv_t *A, const curcell *in, float alpha, stream_t &stream);
void cumuv_trans(curcell **out, float beta, cumuv_t *A, const curcell *in, float alpha, stream_t &stream);

void cuchol_solve(float *restrict out, cusp *Cl, int *Cp, const float *restrict in, 
		  cudaStream_t stream);

void gpu_tomo_test(SIM_T *simu);
void gpu_fit_test(SIM_T *simu);

void gpu_setup_recon_mvm_trans(const PARMS_T *parms, RECON_T *recon, POWFS_T *powfs);
void gpu_setup_recon_mvm_direct(const PARMS_T *parms, RECON_T *recon, POWFS_T *powfs);
void gpu_recon_free_do();
double gpu_fit_do(const PARMS_T *parms,const RECON_T *recon, curcell *fitr, curcell *opdr, stream_t &stream);
double gpu_tomo_do(const PARMS_T *parms,const RECON_T *recon, curcell *opdr, curcell *opdx, curcell *grad, stream_t &stream);
#endif

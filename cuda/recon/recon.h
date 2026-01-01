/*
  Copyright 2009-2026 Lianqi Wang
  
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
#include "../math/cumath.h"
#include "solve.h"
#include "recon_geom.h"
#include "moao.h"
#include "gpu_recon.h"
class curecon_t:public nonCopyable{
public:
    curecon_geom *grid;
    cusolve_r *FR;
    cusolve_l *FL;
    cusolve_r *RR;
    cusolve_l *RL;
    cusolve_l *MVM;
private:
    curcell gradin; /**< The grad to operator on*/
   
    curcell opdr;  /**<Reconstructed atm on xloc. Don't free to have warm restart. Free with new seed*/
    curcell opdr_vec; /**<Referencing opdr in vector form*/
	cumapcell opdr_map;/**<Referencing opdr in cumap form */
    curcell tomo_rhs;
    curcell dmrecon;
    curcell dmrecon_vec;
	curcell dmrecon_zonal;//if dmrecon is modal, convert to zonal.
	curcell dmrecon_zonal_vec;//if dmrecon is modal, convert to zonal.
    curcell fit_rhs;
    stream_t cgstream;
    curcell opdr_save;/*for debugging*/
    curcell dmrecon_save;/*for debugging*/
    curcell GXL;

    curcell gngsmvst;
    curcell deltafocus;
    curcell amod;//if modal control, DM modes.
public:
    int nmoao;
    cumoao_t **moao;/**<moao configurations for GPU*/
    X(mat) *moao_gwfs, *moao_gevl;
    curcccell dm_moao;/**<moao output*/
    curcell dm_wfs;/**<moao results for wfs for warm restart*/
    curcell dm_evl;/**<moao results for evl for warm restart*/
    /*the following data reside in the gpu memory*/
    friend void gpu_update_recon_control(const parms_t *parms, recon_t *recon);
    friend void gpu_update_recon_cn2(const parms_t *parms, recon_t *recon);
public:
    curecon_t(const parms_t *parms=0, recon_t *recon=0);
    void reset_fit(){
        if(FL||FR){
            dbg("reset DM fitting in GPU %d.\n", current_gpu());
            if(FL!=dynamic_cast<cusolve_l *>(FR)) delete FL;
            FL=0;
            delete FR; FR=0;
        }
        dmrecon.Zero();
    }
    void reset_tomo(){
        if(RL||RR){
            dbg("reset tomography in GPU %d.\n", current_gpu());
            if(RL!=dynamic_cast<cusolve_l *>(RR)) delete RL;
            RL=0;
            delete RR; RR=0;
        }
        //opdr.Zero();//no need here. first_run is automatically set in tomo pcg.
    }
    void reset_moao(){
        if(moao){
            dbg("reset MOAO in GPU %d.\n", current_gpu());
            for(int im=0; im<nmoao; im++){
                delete moao[im];
            }
            delete[] moao;
            moao=0; nmoao=0;
        }
    }
    void reset_mvm(){
        if(MVM){
            dbg("reset MVM in GPU %d.\n", current_gpu());
            delete MVM; MVM=0;
        }
    }
    ~curecon_t(){
	    reset_fit();
        reset_tomo();
        reset_moao();
        reset_mvm();
	    delete grid;
    }
    void reset_runtime();
    void init_mvm(const parms_t* parms, recon_t* recon);
    void init_fit(const parms_t* parms, recon_t* recon);
    void init_tomo(const parms_t* parms, recon_t* recon);
    void init_moao(const parms_t* parms, recon_t* recon);
    void update_cn2(const parms_t *parms, recon_t *recon);
    Real tomo(dcell **_opdr, dcell **gngsmvst, const dcell *_gradin);
    Real fit(dcell **_dmrecon, dcell *_opdr);
    Real moao_recon(dcell *_dmrecon, dcell *_opdr);
    void moao_filter(dcell *_dm_wfs, dcell *_dm_evl);
    void mvm(dcell **_dmerr, dcell *_gradin);
    void tomo_test(sim_t *simu);
    void fit_test(sim_t *simu);
	void copy_opdr();
};

curcell new_opdr(const parms_t *parms, const recon_t *recon, Real *p=NULL);
curcell new_dmr(const parms_t *parms, const recon_t *recon, int zonal=0, Real *p=NULL);
void gpu_setup_recon_mvm_trans(const parms_t *parms, recon_t *recon);
void gpu_setup_recon_mvm_direct(const parms_t *parms, recon_t *recon);


#endif

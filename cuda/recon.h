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
#ifndef AOS_CUDA_RECON_H
#define AOS_CUDA_RECON_H

#include "solve.h"
#include "recon_base.h"
#include "moao.h"
namespace cuda_recon{
class curecon_t:public nonCopiable{
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
    curcell tomo_rhs;
    curcell dmfit;
    curcell dmfit_vec;
    curcell fit_rhs;
    stream_t cgstream;
    curcell opdr_save;/*for debugging*/
    curcell dmfit_save;/*for debugging*/

    curcell GXL;

    curcell gngsmvst;
    curcell deltafocus;
    
public:
    int nmoao;
    cumoao_t **moao;/**<moao configurations for GPU*/
    X(mat) *moao_gwfs, *moao_gevl;
    curcccell dm_moao;/**<moao output*/
    curcell dm_wfs;/**<moao results for wfs for warm restart*/
    curcell dm_evl;/**<moao results for evl for warm restart*/
    /*the following data reside in the gpu memory*/
    friend void gpu_update_recon(const PARMS_T *parms, POWFS_T*powfs, RECON_T *recon);
    friend void gpu_update_recon_cn2(const PARMS_T *parms, RECON_T *recon);
public:
    curecon_t(const PARMS_T *parms=0, POWFS_T *powfs=0, RECON_T *recon=0);
    void reset_config();
    ~curecon_t(){
	reset_config();
    }
    void reset_runtime();
    void update(const PARMS_T *parms, POWFS_T*powfs, RECON_T *recon);
    void update_cn2(const PARMS_T *parms, RECON_T *recon);
    Real tomo(dcell **_opdr, dcell **gngsmvst, const dcell *_gradin);
    Real fit(dcell **_dmfit, dcell *_opdr);
    Real moao_recon(dcell *_dmfit, dcell *_opdr);
    void moao_filter(dcell *_dm_wfs, dcell *_dm_evl);
    void mvm(dcell **_dmerr, dcell *_gradin);
    void tomo_test(SIM_T *simu);
    void fit_test(SIM_T *simu);
};

}//namespace

void gpu_setup_recon_mvm_trans(const PARMS_T *parms, RECON_T *recon, POWFS_T *powfs);
void gpu_setup_recon_mvm_direct(const PARMS_T *parms, RECON_T *recon, POWFS_T *powfs);


#endif

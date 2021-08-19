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
/**
   \file maos/utils.h
*/

#ifndef __AOS_UTILS_H
#define __AOS_UTILS_H
#include "common.h"

#if USE_STATIC
extern char _binary____config_tar_gz_start;
extern char _binary____config_tar_gz_end;
#endif

void plotloc(const char *fig, const parms_t *parms, 
	     loc_t *loc, real ht, const char *format,...) CHECK_ARG(5);
void plotdir(const char* fig, const parms_t* parms, real totfov, const char* format, ...) CHECK_ARG(4);
void rename_file(int sig);
int maos_signal_handler(int sig);
arg_t* parse_args(int argc, const char *argv[]);
char *evl_header(const parms_t *parms, const aper_t *aper, int ievl, int iwvl, int isim);
void apply_fieldstop(dmat *opd, const dmat *amp, const lmat *embed, long nembed, const dmat* fieldstop, real wvl);
void display_server(int sock);
void plot_setup(const parms_t *parms, const powfs_t *powfs, const aper_t *aper, const recon_t *recon);
dmat *mkamp(loc_t *loc, map_t *ampground, real misregx, real misregy, real D, real Din);
void maxapriori(real *g, const dmat *ints, const parms_t *parms, 
		const powfs_t *powfs, int iwfs, int isa, int noisy,
		real bkgrnd, real rne);
void wfslinearity(const parms_t *parms, powfs_t *powfs, const int iwfs);
void lgs_wfs_sph_psd(const parms_t *parms, powfs_t *powfs, recon_t *recon, const int iwfs);
real wfsfocusadj(sim_t *simu, int iwfs);
void dither_position(real *cs, real *ss, int alfsm, int dtrat, int npoint, int isim, real deltam);   
void shwfs_grad(dmat **pgrad, dmat *ints[], const parms_t *parms, const powfs_t *powfs, int iwfs, int phytype);
dcell *dcellread_prefix(const char *file, const parms_t *parms, int ipowfs);

/**
   Create first order low pass filter coeffcient from cross over frequency and sampling rate.
*/
static inline real fc2lp(real fc, real dt){
    if(fc*dt>=1){
	return 1;
    }else{
	return 1-exp(-2*M_PI*fc*dt);
    }
}
real average_powfs(dmat* A, lmat* wfsindex, int replace);
#endif

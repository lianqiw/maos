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
#ifndef AOS_CUDA_TOMO_H
#define AOS_CUDA_TOMO_H
#include "solve.h"
#include "recon_geom.h"
#include "fdpcg.h"

typedef struct gpu_gp_t{
	int ipowfs;
	int nwfs; //number of wfs in this group
	int jwfs; //wfs index in this group.
	int(*saptr)[2];//index in ploc for lower left corner of subaperture
	Real* PTTF;
	Real* PFF;
	Real dsa;
	int nsa;
	short2* GPp;
	Real GPscale;
	int pos;
	int nxp;
	Real dxp, dyp;/*pmap dx*/
	Real oxp, oyp;/*pmap origin*/
	const Real(*neai)[3];
	//Must have empty constructor since the array may be in GPU memory.
	gpu_gp_t(){};//:ipowfs(0),nwfs(0),jwfs(0),saptr(0),PTT(0),PDF(0),PDFTT(0),dsa(0),nsa(0),GPp(0),GPscale(0),pos(0),nxp(0),dxp(0),dyp(0),oxp(0),oyp(0){}
}gpu_gp_t;
struct lap_t{
	int nxps, nyps;
	Real l2c;
	int zzi;
	Real zzv;
};
class cutomo_grid:public cusolve_r, public cusolve_cg{
	const curecon_geom* grid;
	curmat wfsrot;
	/*Temporary data*/
	curcell opdwfs;
	curcell opdwfs2;
	curcell grad;
	curmat ttf;

	/*Configuration data*/
	curcell neai;
	curcell PTTF;  /**< Global tip/tilt and optionally focus*/
	curcell PFF;  /**<Focal removal*/
	NumCell<int, Gpu> saptr;
	NumCell<short2, Gpu> GPp;
	Array<Real> GPscale;
	cuspcell GP;
	int rhs_nttf=0;		/**<Number of modes in TTF projection (2 for TT only and 3 with focus)*/
	int lhs_nttf=0;   /**<LHS flag for L(). =0 or nttf*/
	int lhs_skip0=0;	/**<LHS skip powfs[0].wfs[0] to avoid rank deficiency */
	int nwfs=0;
	map2map hx;
	Array<gpu_gp_t, Gpu>gpdata;
	Array<lap_t, Gpu> lap;

	void do_gp(curcell& grad, const curcell& opdwfs, int ptt, int skip0, stream_t& stream);
	void do_gpt(curcell& opdwfs, curcell& grad, int ptt, stream_t& stream);
public:
	cutomo_grid(const parms_t* parms, const recon_t* recon, const curecon_geom* _grid);
	void init_hx(const parms_t* parms, const recon_t* recon);
	void update_fdpcg(fdpcg_t* fdpcg){
		dynamic_cast<cufdpcg_t*>(precond)->update(fdpcg);
	}
	void HX(const curcell& xin, Real alpha, stream_t& stream);
	void HXT(curcell& xout, Real alpha, stream_t& stream);
	virtual void R(curcell& out, Real beta,
		curcell& xin, Real alpha, stream_t& stream);
	virtual void L(curcell& out, Real beta,
		const curcell& xin, Real alpha, stream_t& stream);
	virtual void Rt(curcell& out, Real beta,
		const curcell& xin, Real alpha, stream_t& stream);
	~cutomo_grid(){}
};


#endif

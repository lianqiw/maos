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
/*
  Tests of atmospheric turbulence screen generation (genscreen, genatm_simple,
  genscreen_str). Ported from test_atm.c, which generated a large screen and
  wrote it to file.
*/
#include "unittest.h"

#if CPU_SINGLE
#define TOL 1e-5
#else
#define TOL 1e-10
#endif

/*Configuration of the old test, with a smaller screen.*/
typedef struct{
	rand_t rstat;
	real wt[16];
	real L0[16];
	genatm_t data;
}atmcfg_t;
static void cfg_init(atmcfg_t *cfg, int seed, long nlayer, long nx){
	memset(cfg, 0, sizeof(*cfg));
	rand_t bstat;
	seed_rand(&bstat, seed);
	seed_rand(&cfg->rstat, lrand(&bstat));
	const real wt[7]={0.2887, 0.17795, 0.06602, 0.07833, 0.1405, 0.1216, 0.1269};
	for(int i=0; i<16; i++){
		cfg->wt[i]=wt[i%7];
		cfg->L0[i]=30;
	}
	cfg->data.rstat=&cfg->rstat;
	cfg->data.nx=nx;
	cfg->data.ny=nx;
	cfg->data.dx=1./64.;
	cfg->data.r0=0.1987;
	cfg->data.L0=cfg->L0;
	cfg->data.wt=cfg->wt;
	cfg->data.nlayer=nlayer;
	cfg->data.slope=-11./3.;
}
static real map_mean(const map_t *map){
	return dsum(map->dmat)/PN(map);
}
static real map_var(const map_t *map){
	real mean=map_mean(map);
	return dsumsq(map->dmat)/PN(map)-mean*mean;
}
/*Maximum absolute difference between scaled map A and map B.*/
static real map_diff(const map_t *A, real scale, const map_t *B){
	real maxd=0;
	for(long i=0; i<PN(A); i++){
		real d=fabs(P(A, i)*scale-P(B, i));
		if(d>maxd) maxd=d;
	}
	return maxd;
}

/** genscreen returns the requested number of layers with the right size, sampling and zero mean. */
static void test_layers(void **state){
	(void)state;
	for(int ifractal=0; ifractal<2; ifractal++){
		atmcfg_t cfg;
		const long nx=ifractal?129:128;
		cfg_init(&cfg, 1, 2, nx);
		if(ifractal){
			cfg.data.slope=-1;
			cfg.data.ninit=2;
		}
		mapcell *map=genscreen(&cfg.data);
		assert_non_null(map);
		assert_int_equal(NX(map), 2);
		assert_int_equal(NY(map), 1);
		for(int il=0; il<2; il++){
			const map_t *m=P(map, il);
			assert_non_null(m);
			assert_int_equal(NX(m), nx);
			assert_int_equal(NY(m), nx);
			assert_float_equal(m->dx, 1./64., 0);
			assert_float_equal(m->dy, 1./64., 0);
			const real rms=sqrt(map_var(m));
			assert_true(rms>1e-9 && rms<1e-5);//order of magnitude of a turbulence screen in meter
			assert_float_equal(map_mean(m), 0, 1e-10*rms);//piston is removed
		}
		//layers are different
		assert_true(map_diff(P(map, 0), 1, P(map, 1))>0);
		cellfree(map);
	}
}
/** The screens are reproducible for the same seed and differ for different seeds. */
static void test_seed(void **state){
	(void)state;
	for(int ifractal=0; ifractal<2; ifractal++){
		mapcell *map[3];
		for(int i=0; i<3; i++){
			atmcfg_t cfg;
			cfg_init(&cfg, i<2?1:2, 2, ifractal?65:64);
			if(ifractal){
				cfg.data.slope=-1;
				cfg.data.ninit=2;
			}
			map[i]=genscreen(&cfg.data);
		}
		for(int il=0; il<2; il++){
			assert_float_equal(map_diff(P(map[0], il), 1, P(map[1], il)), 0, 0);
			real rms=sqrt(map_var(P(map[0], il)));
			assert_true(map_diff(P(map[0], il), 1, P(map[2], il))>0.1*rms);
		}
		for(int i=0; i<3; i++){
			cellfree(map[i]);
		}
	}
}
/** FFT screen amplitude scales with sqrt(wt) and r0^(-5/6) exactly for the same random numbers. */
static void test_scaling(void **state){
	(void)state;
	const long nx=64;
	atmcfg_t cfg1, cfg2, cfg3;
	cfg_init(&cfg1, 1, 3, nx);
	cfg_init(&cfg2, 1, 3, nx);
	cfg_init(&cfg3, 1, 3, nx);
	cfg1.wt[0]=1;   cfg1.wt[1]=1;    cfg1.wt[2]=1;
	cfg2.wt[0]=4;   cfg2.wt[1]=0.25; cfg2.wt[2]=0.5;
	cfg3.wt[0]=1;   cfg3.wt[1]=1;    cfg3.wt[2]=1;
	cfg3.data.r0=cfg1.data.r0*2;
	mapcell *map1=genscreen(&cfg1.data);
	mapcell *map2=genscreen(&cfg2.data);
	mapcell *map3=genscreen(&cfg3.data);
	for(int il=0; il<3; il++){
		const real rms=sqrt(map_var(P(map1, il)));
		assert_float_equal(map_diff(P(map1, il), sqrt(cfg2.wt[il]), P(map2, il)), 0, TOL*rms);
		assert_float_equal(map_diff(P(map1, il), pow(2, -5./6.), P(map3, il)), 0, TOL*rms);
	}
	cellfree(map1);
	cellfree(map2);
	cellfree(map3);
}
/** Fractal screen variance scales with wt (statistically, as the initial coarse grid factorization depends on r0). */
static void test_scaling_fractal(void **state){
	(void)state;
	const long nx=65;
	const long nlayer=16;
	atmcfg_t cfg1, cfg2;
	cfg_init(&cfg1, 1, nlayer, nx);
	cfg_init(&cfg2, 1, nlayer, nx);
	for(int il=0; il<nlayer; il++){
		cfg1.wt[il]=1;
		cfg2.wt[il]=4;
	}
	cfg1.data.slope=cfg2.data.slope=-1;
	cfg1.data.ninit=cfg2.data.ninit=2;
	mapcell *map1=genscreen(&cfg1.data);
	mapcell *map2=genscreen(&cfg2.data);
	real var1=0, var2=0;
	for(int il=0; il<nlayer; il++){
		var1+=map_var(P(map1, il));
		var2+=map_var(P(map2, il));
	}
	assert_float_equal(var2/var1, 4, 0.4);
	cellfree(map1);
	cellfree(map2);
}
/** The per-pixel variance of FFT screens agrees with the sum of the PSD and the layers are uncorrelated. */
static void test_variance(void **state){
	(void)state;
	const long nlayer=16;
	const long nx=128;
	atmcfg_t cfg;
	cfg_init(&cfg, 1, nlayer, nx);
	//small outer scale so that many independent modes contribute.
	for(int il=0; il<nlayer; il++){
		cfg.L0[il]=il<8?0.25:0.5;//two different L0 exercises the unpaired layer generation
		cfg.wt[il]=1;
	}
	mapcell *map=genscreen(&cfg.data);
	real neff[2];//effective number of independent modes
	for(int iL0=0; iL0<2; iL0++){
		dmat *spect=turbpsd(nx, nx, cfg.data.dx, cfg.data.r0, iL0?0.5:0.25, cfg.data.slope, 0.5);
		const real var_theory=dsumsq(spect);
		dcwpow(spect, 2);
		neff[iL0]=var_theory*var_theory/dsumsq(spect);
		dfree(spect);
		real var=0;
		for(int il=iL0*8; il<iL0*8+8; il++){
			var+=map_var(P(map, il));
		}
		var/=8;
		assert_float_equal(var/var_theory, 1, 0.1);
	}
	//correlation between layers is consistent with zero (standard deviation is about 1/sqrt(neff))
	for(int il=0; il+1<nlayer; il++){
		real cc=ddot(P(map, il)->dmat, P(map, il+1)->dmat)/sqrt(dsumsq(P(map, il)->dmat)*dsumsq(P(map, il+1)->dmat));
		assert_float_equal(cc, 0, 4/sqrt(neff[il<8?0:1]));
	}
	cellfree(map);
}
/** genatm_simple and genscreen_str generate the same screen as genscreen with the same seed. */
static void test_simple(void **state){
	(void)state;
	const real r0=0.2, L0=30, dx=1./16.;
	const long nx=64;
	const int seed=5;
	map_t *atm=genatm_simple(r0, L0, -11./3., dx, nx, seed);
	rand_t rstat;
	seed_rand(&rstat, seed);
	real wt=1;
	real L0v=L0;
	genatm_t data;
	memset(&data, 0, sizeof(data));
	data.rstat=&rstat; data.wt=&wt; data.r0=r0; data.L0=&L0v; data.dx=dx;
	data.slope=-11./3.; data.nx=nx; data.ny=nx; data.nlayer=1;
	mapcell *map=genscreen(&data);
	assert_float_equal(map_diff(atm, 1, P(map, 0)), 0, 0);
	mapcell *map2=genscreen_str("r0=0.2;L0=30;dx=0.0625;nx=64;seed=5;ht=1000;vx=10;vy=-5");
	assert_non_null(map2);
	assert_int_equal(PN(map2), 1);
	assert_int_equal(NX(P(map2, 0)), nx);
	assert_float_equal(P(map2, 0)->ht, 1000, 0);
	assert_float_equal(P(map2, 0)->vx, 10, 0);
	assert_float_equal(P(map2, 0)->vy, -5, 0);
	assert_float_equal(map_diff(atm, 1, P(map2, 0)), 0, 0);
	cellfree(map2);
	cellfree(map);
	mapfree(atm);
}
/** Spatially varying r0 (r0logpsds) multiplies the screen by a positive smooth field. */
static void test_r0log(void **state){
	(void)state;
	atmcfg_t cfg1, cfg2;
	cfg_init(&cfg1, 3, 2, 64);
	cfg_init(&cfg2, 3, 2, 64);
	dmat *r0logpsds=dnew(2, 1);
	P(r0logpsds, 0)=-2;
	P(r0logpsds, 1)=1e-3;
	cfg2.data.r0logpsds=r0logpsds;
	mapcell *map1=genscreen(&cfg1.data);
	mapcell *map2=genscreen(&cfg2.data);
	for(int il=0; il<2; il++){
		real rmin=INFINITY, rmax=-INFINITY;
		const map_t *m1=P(map1, il), *m2=P(map2, il);
		for(long i=0; i<PN(m1); i++){
			if(fabs(P(m1, i))>1e-3*sqrt(map_var(m1))){
				real ratio=P(m2, i)/P(m1, i);
				assert_true(ratio>0);
				if(ratio<rmin) rmin=ratio;
				if(ratio>rmax) rmax=ratio;
			}
		}
		assert_true(rmax>rmin*1.01);//the scaling varies across the screen
	}
	cellfree(map1);
	cellfree(map2);
	dfree(r0logpsds);
}

UNITTEST_MAIN(
	cmocka_unit_test(test_layers),
	cmocka_unit_test(test_seed),
	cmocka_unit_test(test_scaling),
	cmocka_unit_test(test_scaling_fractal),
	cmocka_unit_test(test_variance),
	cmocka_unit_test(test_simple),
	cmocka_unit_test(test_r0log)
)

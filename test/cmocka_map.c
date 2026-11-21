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
  Unit tests for map_t and rmap_t (math/map.c) and the conversion between maps
  and loc_t in math/loc.c (loc_create_map, loc_embed, loc_extract,
  loc_from_map, ...).
*/
#include <unistd.h>
#include "unittest.h"

/** mapnew, mapnew2, mapref, mapdup and map_convert. */
static void map_create(void **state){
	(void)state;
	map_t *A=mapnew(10, 7, 0.5, 0.25);
	assert_non_null(A);
	assert_int_equal(A->id, M_MAP);
	assert_int_equal(NX(A), 10);
	assert_int_equal(NY(A), 7);
	assert_float_equal(A->ox, -2.5, 0);
	assert_float_equal(A->oy, -0.75, 0);
	assert_float_equal(A->dx, 0.5, 0);
	assert_float_equal(A->dy, 0.25, 0);
	assert_float_equal(dsum(A->dmat), 0, 0);
	A->ht=1000; A->vx=3; A->vy=-4; A->iac=0.3; A->dratio=1.5; A->aoi=0.1;
	A->ox=-1.25;
	P(A, 3, 2)=7;
	//mapnew2 copies the metadata only
	map_t *B=mapnew2(A);
	assert_int_equal(NX(B), 10);
	assert_float_equal(B->ox, -1.25, 0);
	assert_float_equal(B->ht, 1000, 0);
	assert_float_equal(B->vy, -4, 0);
	assert_float_equal(B->iac, 0.3, 0);
	assert_float_equal(B->dratio, 1.5, 0);
	assert_float_equal(B->aoi, 0.1, 0);
	assert_float_equal(dsum(B->dmat), 0, 0);
	//mapref shares the data
	map_t *C=mapref(A);
	assert_ptr_equal(P(C), P(A));
	assert_float_equal(C->vx, 3, 0);
	assert_float_equal(C->dy, 0.25, 0);
	assert_int_equal(mem_nref(A->mem), 2);
	//mapdup copies the data
	map_t *D=mapdup(A);
	assert_ptr_not_equal(P(D), P(A));
	assert_float_equal(P(D, 3, 2), 7, 0);
	assert_float_equal(D->ht, 1000, 0);
	P(D, 3, 2)=1;
	assert_float_equal(P(A, 3, 2), 7, 0);
	assert_null(mapnew2(NULL));
	assert_null(mapref(NULL));
	assert_null(mapdup(NULL));
	map_t *E=mapnew(0, 0, 1, 1);
	assert_null(mapnew2(E));
	mapfree(E);
	//map_convert only accepts real matrices
	assert_null(map_convert(NULL));
	assert_null(rmap_convert(NULL));
	cmat *cm=cnew(2, 2);
	assert_null(map_convert((dmat *)cm));
	assert_null(rmap_convert((dmat *)cm));
	cfree(cm);
	mapfree(A); mapfree(B); mapfree(C); mapfree(D);
	assert_null(A);
	//mapcell
	mapcell *mc=mapcellnew(2, 1);
	P(mc, 0)=mapnew(4, 4, 1, 1);
	P(mc, 1)=mapnew(2, 2, 1, 1);
	cellfree(mc);
}

/** mapcircle area and map_d_din. */
static void map_circle(void **state){
	(void)state;
	const real dx=1./16;
	map_t *A=mapnew(64, 64, dx, dx);
	real R=1.5;
	mapcircle(A, R, 2);
	assert_float_equal(dsum(A->dmat)*dx*dx, 2*M_PI*R*R, 0.01*M_PI*R*R);
	//annulus
	mapcircle(A, 0.5, -2);
	real d=0, din=0;
	map_d_din(A, &d, &din);
	assert_float_equal(d, 2*R, 4*dx);
	assert_float_equal(din, 1, 4*dx);
	//invalid input is ignored
	d=din=-1;
	map_d_din(NULL, &d, &din);
	map_d_din(A, NULL, &din);
	assert_float_equal(d, -1, 0);
	assert_float_equal(din, -1, 0);
	mapcircle(NULL, 1, 1);
	mapfree(A);
}

/** map_make_keywords and map_parse_keywords round trip and defaults. */
static void map_keywords(void **state){
	(void)state;
	map_t *A=mapnew(8, 6, 0.2, 0.3);
	A->ox=-0.7; A->oy=0.4; A->ht=500; A->vx=10; A->vy=-5;
	assert_non_null(A->make_keywords);
	A->make_keywords((cell *)A);
	assert_non_null(A->keywords);
	assert_non_null(strstr(A->keywords, "ox=-0.7\n"));
	assert_non_null(strstr(A->keywords, "h=500\n"));
	//regenerate replaces the old string
	A->vx=11;
	map_make_keywords((cell *)A);
	assert_non_null(strstr(A->keywords, "vx=11\n"));
	//convert a copy with keywords back
	dmat *B=dnew(8, 6);
	B->keywords=strdup(A->keywords);
	map_t *Bm=map_convert(B);
	assert_float_equal(Bm->ox, -0.7, 1e-15);
	assert_float_equal(Bm->oy, 0.4, 1e-15);
	assert_float_equal(Bm->dx, 0.2, 1e-15);
	assert_float_equal(Bm->dy, 0.3, 1e-15);
	assert_float_equal(Bm->ht, 500, 0);
	assert_float_equal(Bm->vx, 11, 0);
	assert_float_equal(Bm->vy, -5, 0);
	assert_float_equal(Bm->iac, 0, 0);
	mapfree(Bm);
	//D is used for dx when dx is missing; ox, oy are centered with offset.
	dmat *C=dnew(10, 4);
	C->keywords=strdup("D=2\nh=100\niac=0.3\ndratio=2\naoi=0.5\noffset=0.5\n");
	map_t *Cm=map_convert(C);
	assert_float_equal(Cm->dx, 0.2, 1e-15);
	assert_float_equal(Cm->dy, 0.2, 1e-15);
	assert_float_equal(Cm->ox, (-5+0.5)*0.2, 1e-15);
	assert_float_equal(Cm->oy, (-2+0.5)*0.2, 1e-15);
	assert_float_equal(Cm->ht, 100, 0);
	assert_float_equal(Cm->iac, 0.3, 0);
	assert_float_equal(Cm->dratio, 2, 0);
	assert_float_equal(Cm->aoi, 0.5, 0);
	//neither dx nor D: default to 1/64
	free(Cm->keywords);
	Cm->keywords=strdup("dy=0.1\n");
	map_parse_keywords(Cm);
	assert_float_equal(Cm->dx, 1./64, 1e-15);
	assert_float_equal(Cm->dy, 0.1, 1e-15);
	assert_float_equal(Cm->ox, -5./64, 1e-15);
	assert_float_equal(Cm->oy, -0.2, 1e-15);
	//Nothing to do without keywords or for other types
	free(Cm->keywords); Cm->keywords=NULL;
	Cm->dx=3;
	map_parse_keywords(Cm);
	assert_float_equal(Cm->dx, 3, 0);
	map_parse_keywords(NULL);
	map_make_keywords(NULL);
	rmap_make_keywords((cell *)Cm);//wrong type
	assert_null(Cm->keywords);
	mapfree(Cm);
	mapfree(A);
}

/** DM grid parameters iac, dratio, and aoi survive writing and reading a map (bin and fits). */
static void map_keywords_file(void **state){
	(void)state;
	const char *suffixes[]={".bin", ".fits"};
	for(int is=0; is<2; is++){
		for(int it=0; it<2; it++){
			map_t *A=mapnew(4, 3, 0.5, 0.5);
			A->ht=1000; A->vx=2;
			if(it){
				A->iac=0.3; A->dratio=1.5; A->aoi=0.25;
			}
			int fds[2];
			assert_int_equal(pipe(fds), 0);
			file_t *fp=zfdopen_as(fds[1], "w", suffixes[is]);
			writedata(fp, (cell *)A, 0);//data is small enough to fit in the pipe buffer
			zfclose(fp);
			if(!it){//zero values are not written
				assert_null(strstr(A->keywords, "iac"));
				assert_null(strstr(A->keywords, "dratio"));
				assert_null(strstr(A->keywords, "aoi"));
			}
			fp=zfdopen_as(fds[0], "r", suffixes[is]);
			map_t *B=map_convert(dmat_cast(readbin_fp(fp, M_REAL, -1)));
			zfclose(fp);
			assert_non_null(B);
			assert_float_equal(B->dx, A->dx, 1e-15);
			assert_float_equal(B->ox, A->ox, 1e-15);
			assert_float_equal(B->ht, A->ht, 1e-15);
			assert_float_equal(B->vx, A->vx, 1e-15);
			assert_float_equal(B->iac, A->iac, 1e-15);
			assert_float_equal(B->dratio, A->dratio, 1e-15);
			assert_float_equal(B->aoi, A->aoi, 1e-15);
			mapfree(B);
			mapfree(A);
		}
	}
}

/** rmap keywords round trip and the required fields. */
static void map_rmap(void **state){
	(void)state;
	dmat *A0=dnew(5, 3);
	A0->keywords=strdup("ox=-1\noy=-2\ndx=0.5\ndy=0.25\ntxdeg=30\ntydeg=90\nftel=100\nfexit=90\nfsurf=10\n");
	rmap_t *A=rmap_convert(A0);
	assert_int_equal(A->id, M_RMAP);
	assert_float_equal(A->ox, -1, 0);
	assert_float_equal(A->dy, 0.25, 0);
	assert_float_equal(A->txdeg, 30, 0);
	assert_float_equal(A->tydeg, 90, 0);
	assert_float_equal(A->ftel, 100, 0);
	assert_float_equal(A->fexit, 90, 0);
	assert_float_equal(A->fsurf, 10, 0);
	A->fsurf=12;
	A->make_keywords((cell *)A);
	assert_non_null(strstr(A->keywords, "fsurf=12\n"));
	rmap_parse_keywords(A);
	assert_float_equal(A->fsurf, 12, 0);
	assert_float_equal(A->fexit, 90, 0);
	map_make_keywords((cell *)A);//wrong type: no change
	assert_non_null(strstr(A->keywords, "fsurf=12\n"));
	//missing a required field is reported by error()
	free(A->keywords);
	A->keywords=strdup("ox=1\noy=1\ndx=1\ndy=1\n");
	rmap_parse_keywords(A);
	assert_true(isnan(A->fsurf));
	rmap_parse_keywords(NULL);
	rmapfree(A);
	rmapcell *rc=rmapcellnew(1, 1);
	cellfree(rc);
}

/*Check that the map contains (value 1) all points within radius r-margin of (cx, cy).*/
static void assert_covers(const map_t *map, real cx, real cy, real r){
	long nin=0;
	for(long iy=0; iy<NY(map); iy++){
		real y=map->oy+iy*map->dy;
		for(long ix=0; ix<NX(map); ix++){
			real x=map->ox+ix*map->dx;
			if((x-cx)*(x-cx)+(y-cy)*(y-cy)<r*r){
				assert_float_equal(P(map, ix, iy), 1, 0);
				nin++;
			}
		}
	}
	assert_true(nin>0);
}
/** create_metapupil for the different grid types and options. */
static void map_metapupil(void **state){
	(void)state;
	const real D=2, dx=0.25;
	map_t *map=NULL;
	long nx=0, ny=0;
	//ground layer, minimal grid
	create_metapupil(&map, &nx, &ny, NULL, D, 0, dx, dx, 0, 0, 0, 0, 0, 0);
	assert_int_equal(nx, 9);
	assert_int_equal(ny, 9);
	assert_int_equal(NX(map), 9);
	assert_float_equal(map->ox, -1, 0);
	assert_float_equal(map->oy, -1, 0);
	assert_float_equal(map->ht, 0, 0);
	assert_covers(map, 0, 0, D/2-dx);
	for(long i=0; i<PN(map); i++){
		assert_true(P(map, i)==0||P(map, i)==1);
	}
	assert_float_equal(P(map, 0, 0), 0, 0);//corner is outside
	mapfree(map);
	//full square grid
	create_metapupil(&map, &nx, &ny, NULL, D, 0, dx, dx, 0, 0, 0, 0, 0, 1);
	assert_float_equal(dsum(map->dmat), 81, 0);
	mapfree(map);
	//circular grid
	create_metapupil(&map, &nx, &ny, NULL, D, 0, dx, dx, 0, 0, 0, 0, 0, 2);
	assert_float_equal(dsum(map->dmat)*dx*dx, M_PI, 0.1*M_PI);
	mapfree(map);
	//half pixel offset adds one point
	create_metapupil(&map, &nx, &ny, NULL, D, 0, dx, dx, 0.5, 0, 0, 0, 0, 0);
	assert_int_equal(nx, 10);
	assert_float_equal(map->ox, -1.125, 1e-15);
	mapfree(map);
	//offset is used modulo 1
	create_metapupil(&map, &nx, &ny, NULL, D, 0, dx, dx, 2.5, 0, 0, 0, 0, 0);
	assert_int_equal(nx, 10);
	mapfree(map);
	//guard and power of 2 padding
	create_metapupil(&map, &nx, &ny, NULL, D, 0, dx, dx, 0, 0.5, 0, 0, 1, 0);
	assert_int_equal(nx, 16);
	assert_int_equal(ny, 16);
	assert_float_equal(map->ox, -1.5-(16-13)/2*dx, 1e-15);
	assert_covers(map, 0, 0, D/2-dx);//guard is not used without directions
	mapfree(map);
	//requested size larger (and smaller with a warning)
	create_metapupil(NULL, &nx, &ny, NULL, D, 0, dx, dx, 0, 0, 20, 12, 0, 0);
	assert_int_equal(nx, 20);
	assert_int_equal(ny, 12);
	create_metapupil(&map, NULL, NULL, NULL, D, 0, dx, dx, 0, 0, 4, 4, 0, 1);
	assert_int_equal(NX(map), 4);
	assert_int_equal(NY(map), 4);
	mapfree(map);

	//two directions at altitude: the map covers both footprints
	dmat *dirs=dnew(4, 2);
	const real ht=10000;
	P(dirs, 0, 0)=20/206265.; P(dirs, 1, 0)=0; P(dirs, 2, 0)=INFINITY; P(dirs, 3, 0)=0;
	P(dirs, 0, 1)=-10/206265.; P(dirs, 1, 1)=15/206265.; P(dirs, 2, 1)=90000; P(dirs, 3, 1)=0;
	create_metapupil(&map, &nx, &ny, dirs, D, ht, dx, dx, 0, 0, 0, 0, 0, 0);
	assert_float_equal(map->ht, ht, 0);
	for(int idir=0; idir<2; idir++){
		real xc=P(dirs, 0, idir)*ht;
		real yc=P(dirs, 1, idir)*ht;
		real RR=(1-ht/P(dirs, 2, idir))*D/2;
		assert_covers(map, xc, yc, RR-dx);
		assert_true(map->ox<=xc-RR+1e-12);
		assert_true(map->ox+(nx-1)*dx>=xc+RR-1e-12);
		assert_true(map->oy<=yc-RR+1e-12);
		assert_true(map->oy+(ny-1)*dx>=yc+RR-1e-12);
	}
	mapfree(map);
	//square=1 and square=2 with directions
	create_metapupil(&map, &nx, &ny, dirs, D, ht, dx, dx, 0, 0, 0, 0, 0, 1);
	assert_int_equal(NX(map), NY(map));
	assert_float_equal(dsum(map->dmat), PN(map), 0);
	mapfree(map);
	create_metapupil(&map, &nx, &ny, dirs, D, ht, dx, dx, 0, 0, 0, 0, 0, 2);
	assert_true(dsum(map->dmat)>0);
	assert_true(dsum(map->dmat)<PN(map));
	mapfree(map);
	//invalid square: error and empty map
	create_metapupil(&map, &nx, &ny, NULL, D, 0, dx, dx, 0, 0, 0, 0, 0, 3);
	assert_float_equal(dsum(map->dmat), 0, 0);
	mapfree(map);
	//altitude without directions: error, falls back to ground layer
	create_metapupil(NULL, &nx, &ny, NULL, D, ht, dx, dx, 0, 0, 0, 0, 0, 0);
	assert_int_equal(nx, 9);
	dfree(dirs);
}

/*Fill a map with value v.*/
static map_t *fillmap(long nx, long ny, real v, real vx, real vy){
	map_t *A=mapnew(nx, ny, 0.5, 0.5);
	dset(A->dmat, v);
	A->vx=vx;
	A->vy=vy;
	return A;
}
/** map_blend blends the overlapping region according to the wind direction. */
static void map_blend_test(void **state){
	(void)state;
	const long nx=8, ny=6, over=3;
	//wind along +x: atm2 is to the left of atm1
	for(int dir=0; dir<5; dir++){
		real vx=0, vy=0;
		switch(dir){
		case 0: vx=1; break;
		case 1: vx=-1; break;
		case 2: vy=1; break;
		case 3: vy=-1; break;
		case 4: break;
		}
		map_t *a1=fillmap(nx, ny, 1, vx, vy);
		map_t *a2=fillmap(nx, ny, 3, vx, vy);
		map_blend(a1, a2, over, over);
		if(dir<2){
			int ca=vx>0?-1:1;
			assert_float_equal(a2->ox, a1->ox+(nx-over)*ca*a1->dx, 1e-15);
			assert_float_equal(a2->oy, a1->oy, 0);
			long off1=ca<0?0:nx-over;//columns of a1 in the overlap
			long off2=ca<0?nx-over:0;
			for(long iy=0; iy<ny; iy++){
				for(long ix=0; ix<over; ix++){
					real w=fabs((ca<0?1:0)-(real)(ix+1)/over);
					real expect=(1-w)*1+w*3;
					assert_float_equal(P(a1, ix+off1, iy), expect, 1e-12);
					assert_float_equal(P(a2, ix+off2, iy), expect, 1e-12);
				}
			}
			//outside of the overlap region unchanged
			assert_float_equal(P(a1, ca<0?nx-1:0, 0), 1, 0);
			assert_float_equal(P(a2, ca<0?0:nx-1, 0), 3, 0);
		}else{
			int sa=vy>0?-1:(vy<0?1:0);
			assert_float_equal(a2->oy, a1->oy+(ny-over)*sa*a1->dy, 1e-15);
			assert_float_equal(a2->ox, a1->ox, 0);
			long off1=sa<0?0:ny-over;
			long off2=sa<0?ny-over:0;
			for(long iy=0; iy<over; iy++){
				real w=fabs((sa<0?1:0)-(real)(iy+1)/over);
				real expect=(1-w)*1+w*3;
				for(long ix=0; ix<nx; ix++){
					assert_float_equal(P(a1, ix, iy+off1), expect, 1e-12);
					assert_float_equal(P(a2, ix, iy+off2), expect, 1e-12);
				}
			}
			assert_float_equal(P(a1, 0, sa<0?ny-1:0), 1, 0);
		}
		mapfree(a1);
		mapfree(a2);
	}
	//diagonal wind is not supported: no change
	map_t *a1=fillmap(nx, ny, 1, 1, 1);
	map_t *a2=fillmap(nx, ny, 3, 1, 1);
	map_blend(a1, a2, over, over);
	assert_float_equal(dsum(a1->dmat), nx*ny, 0);
	assert_float_equal(dsum(a2->dmat), 3*nx*ny, 0);
	mapfree(a1);
	mapfree(a2);
}

/** mksqloc_map, loc_from_map, loc_create_map and loc_map_get. */
static void map_loc(void **state){
	(void)state;
	map_t *A=mapnew(12, 10, 0.5, 0.25);
	A->ox=-2; A->oy=1;
	loc_t *sq=mksqloc_map(A);
	assert_int_equal(sq->nloc, 120);
	assert_float_equal(sq->locx[0], -2, 0);
	assert_float_equal(sq->locy[0], 1, 0);
	assert_float_equal(sq->locx[13], -1.5, 1e-15);
	assert_float_equal(sq->locy[13], 1.25, 1e-15);
	assert_float_equal(sq->dx, 0.5, 0);
	assert_float_equal(sq->dy, 0.25, 0);
	assert_null(mksqloc_map(NULL));
	locfree(sq);
	//loc_from_map collects each row from the first to the last valid point
	A->ht=300; A->iac=0.2;
	P(A, 2, 1)=1; P(A, 7, 1)=0.5;//row 1: 2..7
	P(A, 5, 4)=0.2;//below threshold
	P(A, 0, 6)=2; P(A, 11, 6)=1;//row 6: all
	loc_t *loc=loc_from_map(A, 0.2);
	assert_int_equal(loc->nloc, 6+12);
	assert_float_equal(loc->locx[0], -2+2*0.5, 1e-15);
	assert_float_equal(loc->locy[0], 1+0.25, 1e-15);
	assert_float_equal(loc->locx[5], -2+7*0.5, 1e-15);
	assert_float_equal(loc->locy[6], 1+6*0.25, 1e-15);
	assert_float_equal(loc->ht, 300, 0);
	assert_float_equal(loc->iac, 0.2, 0);
	locfree(loc);
	loc=loc_from_map(A, 0);//all positive points
	assert_int_equal(loc->nloc, 6+1+12);
	locfree(loc);
	assert_null(loc_from_map(NULL, 0));
	map_t *Z=mapnew(3, 3, 1, 1);
	loc=loc_from_map(Z, 0);//no valid point
	assert_int_equal(loc->nloc, 0);
	locfree(loc);
	mapfree(Z);
	mapfree(A);

	//loc_create_map stores index+1 at loc points
	loc_t *ann=mkannloc(2, 0.6, 0.25, 0.5);
	loc_create_map(ann);
	assert_non_null(ann->map);
	map_t *m=ann->map;
	for(long iloc=0; iloc<ann->nloc; iloc++){
		long ix=(long)round((ann->locx[iloc]-m->ox)/m->dx);
		long iy=(long)round((ann->locy[iloc]-m->oy)/m->dy);
		assert_int_equal(loc_map_get(m, ix, iy), iloc+1);
	}
	assert_int_equal(loc_map_get(m, -1, 0), 0);
	assert_int_equal(loc_map_get(m, 0, NY(m)), 0);
	//vacant points next to valid points are extended with negative index of a neighbor
	long nneg=0;
	for(long i=0; i<PN(m); i++){
		if(P(m, i)<0){
			nneg++;
			long iphi=(long)(-P(m, i))-1;
			long ix=i%NX(m), iy=i/NX(m);
			real x=m->ox+ix*m->dx, y=m->oy+iy*m->dy;
			assert_true(fabs(ann->locx[iphi]-x)<=m->dx*1.01);
			assert_true(fabs(ann->locy[iphi]-y)<=m->dy*1.01);
		}
	}
	assert_true(nneg>0);
	const long nx0=NX(m), ny0=NY(m);
	const real ox0=m->ox, oy0=m->oy;
	loc_create_map(ann);//existing map is kept
	assert_ptr_equal(ann->map, m);
	loc_create_map_npad(ann, 2, 0, 0);//cannot redo with larger padding
	assert_ptr_equal(ann->map, m);
	loc_free_map(ann);
	assert_null(ann->map);
	//padding and requested size
	loc_create_map_npad(ann, 1, 0, 0);
	assert_int_equal(ann->npad, 1);
	assert_int_equal(NX(ann->map), nx0+2);
	assert_int_equal(NY(ann->map), ny0+2);
	assert_float_equal(ann->map->ox, ox0-0.25, 1e-12);
	loc_free_map(ann);
	loc_create_map_npad(ann, 0, 15, 13);
	assert_int_equal(NX(ann->map), 15);
	assert_int_equal(NY(ann->map), 13);
	assert_int_equal(ann->npad, MIN((15-nx0)/2, (13-ny0)/2));
	assert_float_equal(ann->map->ox, ox0-(15-nx0)/2*0.25, 1e-12);
	assert_float_equal(ann->map->oy, oy0-(13-ny0)/2*0.25, 1e-12);
	loc_free_map(ann);
	locfree(ann);
	loc_create_map(NULL);
}

/** loc_embed, loc_extract, loc_embed2 and loc_embed_cell. */
static void map_embed(void **state){
	(void)state;
	loc_t *loc=mkannloc(2, 0.5, 0.25, 0.5);
	dmat *opd=dnew(loc->nloc, 1);
	for(long i=0; i<loc->nloc; i++){
		P(opd, i)=loc->locx[i]*2+loc->locy[i]*loc->locy[i]+1;
	}
	//embed into an empty array: resized to the map
	dmat *dest=dnew(0, 0);
	loc_embed(dest, loc, opd, 0, 0);
	assert_non_null(loc->map);
	map_t *m=loc->map;
	assert_int_equal(NX(dest), NX(m));
	assert_int_equal(NY(dest), NY(m));
	long nvalid=0;
	for(long i=0; i<PN(m); i++){
		if(P(m, i)>0){
			assert_float_equal(P(dest, i), P(opd, (long)P(m, i)-1), 0);
			nvalid++;
		}else{
			assert_true(isnan(P(dest, i)));
		}
	}
	assert_int_equal(nvalid, loc->nloc);
	//add mode doubles the valid points
	loc_embed(dest, loc, opd, 0, 1);
	for(long i=0; i<PN(m); i++){
		if(P(m, i)>0) assert_float_equal(P(dest, i), 2*P(opd, (long)P(m, i)-1), 1e-15);
	}
	//expanded fills ghost points with the neighbor's value
	loc_embed(dest, loc, opd, 1, 0);
	for(long i=0; i<PN(m); i++){
		if(P(m, i)!=0){
			assert_float_equal(P(dest, i), P(opd, (long)fabs(P(m, i))-1), 0);
		}else{
			assert_true(isnan(P(dest, i)));
		}
	}
	//extract is the inverse
	map_t *in=mapnew2(m);
	loc_embed(in->dmat, loc, opd, 0, 1);
	dmat *opd2=dnew(loc->nloc, 1);
	loc_extract(opd2, loc, in);
	assert_dmat_near(opd2, opd, 0);
	loc_extract(NULL, loc, in);//no-op
	//larger and smaller destination are concentric
	dmat *big=dnew(NX(m)+4, NY(m)+2);
	loc_embed(big, loc, opd, 0, 0);
	loc_embed(dest, loc, opd, 0, 0);
	assert_false(isnan(P(dest, 1, 4)));
	assert_float_equal(P(big, 2+1, 1+4), P(dest, 1, 4), 0);
	dmat *small=dnew(NX(m)-2, NY(m)-2);
	loc_embed(small, loc, opd, 0, 0);
	assert_float_equal(P(small, 0, 3), P(dest, 1, 4), 0);
	loc_embed(small, NULL, opd, 0, 0);//no-op
	//loc_embed2 embeds each column
	dmat *opd3=dnew(loc->nloc, 3);
	for(long ic=0; ic<3; ic++){
		for(long i=0; i<loc->nloc; i++) P(opd3, i, ic)=P(opd, i)*(ic+1);
	}
	dcell *e2=loc_embed2(loc, opd3);
	assert_int_equal(NX(e2), 1);
	assert_int_equal(NY(e2), 3);
	for(long ic=0; ic<3; ic++){
		for(long i=0; i<PN(m); i++){
			if(P(m, i)>0) assert_float_equal(P(P(e2, ic), i), P(opd, (long)P(m, i)-1)*(ic+1), 1e-12);
		}
	}
	dcellfree(e2);
	//row vector input
	dmat *opdt=dtrans(opd);
	e2=loc_embed2(loc, opdt);
	assert_int_equal(NX(e2), 1);
	for(long i=0; i<PN(m); i++){
		if(P(m, i)>0) assert_float_equal(P(P(e2, 0), i), P(opd, (long)P(m, i)-1), 0);
	}
	dcellfree(e2);
	assert_null(loc_embed2(NULL, opd));
	//loc_embed_cell places references in a map shaped cell
	dcell *inc=dcellnew(loc->nloc, 1);
	for(long i=0; i<loc->nloc; i++){
		P(inc, i)=dnew(1, 1);
		P(P(inc, i), 0)=(real)i;
	}
	dcell *outc=NULL;
	loc_embed_cell(&outc, loc, inc);
	assert_int_equal(NX(outc), NX(m));
	for(long i=0; i<PN(m); i++){
		if(P(m, i)>0){
			assert_ptr_equal(P(P(outc, i)), P(P(inc, (long)P(m, i)-1)));
		}else{
			assert_null(P(outc, i));
		}
	}
	dcellfree(outc);
	dcellfree(inc);
	dfree(opdt); dfree(opd3); dfree(small); dfree(big); dfree(opd2);
	mapfree(in);
	dfree(dest); dfree(opd);
	locfree(loc);
}

/** loc_coord2ind maps coordinates to loc index. */
static void map_coord2ind(void **state){
	(void)state;
	loc_t *loc=mksqloc(5, 4, 0.5, 0.5, -1, -1);
	dmat *dead=dnew(3, 2);
	P(dead, 0, 0)=-1; P(dead, 0, 1)=-1;//index 0
	P(dead, 1, 0)=0.5; P(dead, 1, 1)=0.02;//index 3+2*5=13
	P(dead, 2, 0)=10; P(dead, 2, 1)=10;//outside
	lmat *ind=loc_coord2ind(loc, dead);
	assert_int_equal(PN(ind), 20);
	assert_int_equal(lsum(ind), 2);
	assert_int_equal(P(ind, 0), 1);
	assert_int_equal(P(ind, 13), 1);
	lfree(ind);
	//with values (in m) converted to nm
	dmat *dead3=dnew(1, 3);
	P(dead3, 0, 0)=0; P(dead3, 0, 1)=0.5; P(dead3, 0, 2)=2e-9;
	ind=loc_coord2ind(loc, dead3);
	assert_int_equal(P(ind, 2+3*5), 2);
	assert_int_equal(lsum(ind), 2);
	lfree(ind);
	//a single row vector is a single point
	dmat *v=dnew(2, 1);
	P(v, 0)=1; P(v, 1)=0.5;
	ind=loc_coord2ind(loc, v);
	assert_int_equal(P(ind, 4+3*5), 1);
	assert_int_equal(lsum(ind), 1);
	lfree(ind);
	dfree(v); dfree(dead3); dfree(dead);
	locfree(loc);
}

UNITTEST_MAIN(
	cmocka_unit_test(map_create),
	cmocka_unit_test(map_circle),
	cmocka_unit_test(map_keywords),
	cmocka_unit_test(map_keywords_file),
	cmocka_unit_test(map_rmap),
	cmocka_unit_test(map_metapupil),
	cmocka_unit_test(map_blend_test),
	cmocka_unit_test(map_loc),
	cmocka_unit_test(map_embed),
	cmocka_unit_test(map_coord2ind),
)

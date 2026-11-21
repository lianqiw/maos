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
  Test the bin and fits file formats, plain and compressed. Data is written by
  a thread into a pipe and read back from the other end, so that no temporary
  file is used.
*/
#include <pthread.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include "unittest.h"

static const char *suffixes[]={".bin", ".bin.gz", ".fits", ".fits.gz",
#if HAVE_LIBZSTD
	".bin.zst", ".fits.zst",
#endif
};
#define NSUFFIX (int)(sizeof(suffixes)/sizeof(suffixes[0]))
static const char *keywords="dx=1/64\ndy=1/64\nLong=kaflasdjfl ksflksjlfjasdflkj asldfj als fals fl asfjalsdf jalsf djasldf jal fsalfkjasdflkj asldfkj asldf \nadf\nwvf=1.e-4";

typedef struct{
	int fd;
	const char *suffix;
	const cell *A;
}writer_t;
/**
   Thread that writes A to the write end of the pipe and closes it.
*/
static void *writer(void *data){
	writer_t *w=(writer_t *)data;
	file_t *fp=zfdopen_as(w->fd, "w", w->suffix);
	if(fp){
		writedata(fp, w->A, 0);
		zfclose(fp);
	}else{
		close(w->fd);
	}
	return NULL;
}
/**
   Write A in the format given by suffix and read it back as type id (0 to keep the type).
*/
static cell *roundtrip(const_anyarray A, const char *suffix, M_ID id){
	int fds[2];
	if(pipe(fds)){
		fail_msg("pipe failed: %s", strerror(errno));
	}
	writer_t w={fds[1], suffix, A.c};
	pthread_t thread;
	pthread_create(&thread, NULL, writer, &w);
	file_t *fp=zfdopen_as(fds[0], "r", suffix);
	cell *B=readbin_fp(fp, id, -1);
	//drain any remaining data so that the writer is not blocked.
	char buf[1024];
	if(fp) while(zfread(buf, 1, sizeof(buf), fp)==0);
	pthread_join(thread, NULL);
	zfclose(fp);
	return B;
}
/**
   Check that the keywords are preserved. The bin format keeps them exactly. The
   fits format stores them as upper case cards (KEY = VALUE), so the values are
   compared ignoring case.
*/
static void assert_keywords(const char *kb, const char *suffix){
	assert_non_null(kb);
	if(!strstr(suffix, ".fits")){
		assert_string_equal(kb, keywords);
	}else{
		char *lower=strdup(kb);
		for(char *p=lower; *p; p++) *p=(char)tolower((unsigned char)*p);
		assert_non_null(strstr(lower, "dx      = 1/64"));
		assert_non_null(strstr(lower, "dy      = 1/64"));
		assert_non_null(strstr(lower, "wvf     = 1.e-4"));
		//long values are continued in CONTINUE cards without losing characters
		assert_non_null(strstr(lower, "long    = kaflasdjfl ksflksjlfjasdflkj asldfj als fals fl asfjalsdf jalsf djasldf jal fsalfkjasdflkj asldfkj asldf\n"));
		//text without = is kept as a separate line
		assert_non_null(strstr(kb, "\nadf\n"));
		free(lower);
	}
}
/** zfdopen_as rejects invalid arguments. */
static void fits_zfdopen(void **state){
	(void)state;
	assert_null(zfdopen_as(-1, "r", ".bin"));
	int fds[2];
	assert_int_equal(pipe(fds), 0);
	assert_null(zfdopen_as(fds[0], "x", ".bin"));
	assert_null(zfdopen_as(fds[0], NULL, ".bin"));
	file_t *fp=zfdopen_as(fds[0], "r", ".fits.gz");
	assert_non_null(fp);
	assert_int_equal(zfisfits(fp), 1);
	zfclose(fp);
	fp=zfdopen_as(fds[1], "w", NULL);
	assert_non_null(fp);
	assert_int_equal(zfisfits(fp), 0);
	zfclose(fp);
}
/** dmat with keywords is unchanged after writing and reading in all formats. */
static void fits_dmat(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	dmat *A=dnew(4, 5);
	drandn(A, 1, &rstat);
	A->keywords=strdup(keywords);
	for(int is=0; is<NSUFFIX; is++){
		dmat *B=dmat_cast(roundtrip(A, suffixes[is], M_REAL));
		assert_non_null(B);
		assert_dmat_near(A, B, 0);
		assert_keywords(B->keywords, suffixes[is]);
		dfree(B);
	}
	dfree(A);
}
/** Keywords longer than a fits card are continued without loss, and separate comment lines are kept. */
static void fits_longkey(void **state){
	(void)state;
	char value[200]={0};
	for(int i=0; i<18; i++) strcat(value, "abcdefghi ");//a space falls at the end of each card
	value[strlen(value)-1]='\0';
	char *kw=stradd("first comment\nk1=", value, "\nsecond comment;k2=3", NULL);
	dmat *A=dnew(2, 2);
	A->keywords=kw;
	char *expect=stradd("k1      = ", value, NULL);
	for(char *p=expect; *p; p++) *p=(char)toupper((unsigned char)*p);
	for(int is=0; is<NSUFFIX; is++){
		dmat *B=dmat_cast(roundtrip(A, suffixes[is], M_REAL));
		assert_non_null(B);
		if(strstr(suffixes[is], ".fits")){
			assert_non_null(strstr(B->keywords, expect));
			assert_non_null(strstr(B->keywords, "first comment\nK1"));
			assert_non_null(strstr(B->keywords, "\nsecond comment\nK2      = 3"));
			assert_float_equal(search_keyword_num(B->keywords, "K2"), 3, 0);
		}else{
			assert_string_equal(B->keywords, A->keywords);
		}
		dfree(B);
	}
	free(expect);
	dfree(A);
}
/** smat, cmat, and lmat are unchanged after writing and reading in all formats. */
static void fits_types(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 2);
	smat *As=snew(3, 7);
	srandn(As, 1, &rstat);
	cmat *Ac=cnew(5, 3);
	crandn(Ac, 1, &rstat);
	lmat *Al=lnew(6, 2);
	for(long i=0; i<PN(Al); i++){
		P(Al, i)=(i-5)*1000000007L;
	}
	for(int is=0; is<NSUFFIX; is++){
		smat *Bs=smat_cast(roundtrip(As, suffixes[is], M_FLT));
		assert_non_null(Bs);
		assert_int_equal(NX(Bs), NX(As));
		assert_int_equal(NY(Bs), NY(As));
		for(long i=0; i<PN(As); i++){
			assert_true(P(Bs, i)==P(As, i));
		}
		sfree(Bs);
		cmat *Bc=cmat_cast(roundtrip(Ac, suffixes[is], M_COMP));
		assert_non_null(Bc);
		assert_int_equal(NX(Bc), NX(Ac));
		assert_int_equal(NY(Bc), NY(Ac));
		assert_float_equal(creal(cdiff(Ac, Bc)), 0, 0);
		cfree(Bc);
		lmat *Bl=lmat_cast(roundtrip(Al, suffixes[is], M_LONG));
		assert_non_null(Bl);
		assert_int_equal(PN(Bl), PN(Al));
		for(long i=0; i<PN(Al); i++){
			assert_true(P(Bl, i)==P(Al, i));
		}
		lfree(Bl);
	}
	sfree(As);
	cfree(Ac);
	lfree(Al);
}
/** Reading converts between real types: dmat read as smat and smat read as dmat. */
static void fits_convert(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 3);
	dmat *A=dnew(4, 4);
	drandn(A, 1, &rstat);
	for(int is=0; is<NSUFFIX; is++){
		smat *B=smat_cast(roundtrip(A, suffixes[is], M_FLT));
		assert_non_null(B);
		assert_int_equal(PN(B), PN(A));
		for(long i=0; i<PN(A); i++){
			assert_float_equal(P(B, i), P(A, i), 1e-6*fabs(P(A, i))+1e-30);
		}
		dmat *C=dmat_cast(roundtrip(B, suffixes[is], M_DBL));
		assert_non_null(C);
		for(long i=0; i<PN(A); i++){
			assert_true(P(C, i)==(double)P(B, i));
		}
		sfree(B);
		dfree(C);
	}
	dfree(A);
}
/** dcell with keywords and an empty element is unchanged after writing and reading in all formats. */
static void fits_cell(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	dcell *A=dcellnew(2, 2);
	P(A, 0)=dnew(4, 4);
	P(A, 2)=dnew(3, 5);
	P(A, 3)=dnew(1, 1);
	drandn(P(A, 0), 1, &rstat);
	drandn(P(A, 2), 1, &rstat);
	drandn(P(A, 3), 1, &rstat);
	A->keywords=strdup(keywords);
	for(int is=0; is<NSUFFIX; is++){
		dcell *B=dcell_cast(roundtrip(A, suffixes[is], M_REAL));
		assert_non_null(B);
		if(strstr(suffixes[is], ".fits")){//fits keeps the elements in order but not the cell shape.
			assert_int_equal(PN(B), PN(A));
		}else{
			assert_int_equal(NX(B), NX(A));
			assert_int_equal(NY(B), NY(A));
		}
		assert_true(!P(B, 1)||PN(P(B, 1))==0);
		for(int i=0; i<PN(A); i++){
			if(P(A, i)) assert_dmat_near(P(A, i), P(B, i), 0);
		}
		//fits copies the cell keywords to each element.
		assert_keywords(strstr(suffixes[is], ".fits")?P(B, 0)->keywords:B->keywords, suffixes[is]);
		dcellfree(B);
	}
	dcellfree(A);
}
/** Nested cell arrays of different types are unchanged in the bin format. */
static void fits_nested(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 4);
	cell *A=cellnew(3, 1);
	dcell *A0=dcellnew_same(2, 1, 3, 2);
	for(int i=0; i<PN(A0); i++) drandn(P(A0, i), 1, &rstat);
	ccell *A1=ccellnew_same(1, 2, 2, 2);
	for(int i=0; i<PN(A1); i++) crandn(P(A1, i), 1, &rstat);
	P(A, 0)=(cell *)A0;
	P(A, 1)=(cell *)A1;
	for(int is=0; is<NSUFFIX; is++){
		if(strstr(suffixes[is], ".fits")) continue;//fits does not support nested cells
		cell *B=roundtrip(A, suffixes[is], 0);
		assert_non_null(B);
		assert_int_equal(PN(B), 3);
		assert_null(P(B, 2));
		assert_float_equal(dcelldiff(A0, dcell_cast(P(B, 0))), 0, 0);
		ccell *B1=ccell_cast(P(B, 1));
		assert_non_null(B1);
		assert_int_equal(PN(B1), PN(A1));
		for(int i=0; i<PN(A1); i++){
			assert_float_equal(creal(cdiff(P(A1, i), P(B1, i))), 0, 0);
		}
		cellfree(B);
	}
	cellfree(A);
}
/** Sparse matrix is unchanged after writing and reading in the bin format. */
static void fits_sparse(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 5);
	dsp *A=dspnewrandu(20, 15, 1, 0.2, &rstat);
	dmat *Af=NULL;
	dspfull(&Af, A, 'n', 1);
	for(int is=0; is<NSUFFIX; is++){
		if(strstr(suffixes[is], ".fits")) continue;//fits does not support sparse
		dsp *B=dsp_cast(roundtrip(A, suffixes[is], 0));
		assert_non_null(B);
		assert_int_equal(NX(B), NX(A));
		assert_int_equal(NY(B), NY(A));
		dmat *Bf=NULL;
		dspfull(&Bf, B, 'n', 1);
		assert_dmat_near(Af, Bf, 0);
		dfree(Bf);
		dspfree(B);
	}
	dfree(Af);
	dspfree(A);
}
/** Empty arrays and cells are read back as empty. */
static void fits_empty(void **state){
	(void)state;
	for(int is=0; is<NSUFFIX; is++){
		cell *B=roundtrip((cell *)NULL, suffixes[is], 0);
		assert_null(B);
		dmat *A=dnew(0, 0);
		B=roundtrip(A, suffixes[is], M_REAL);
		assert_true(!B||PN(B)==0);
		cellfree(B);
		dfree(A);
		dcell *Ac=dcellnew(2, 1);//cell with only empty elements
		B=roundtrip(Ac, suffixes[is], M_REAL);
		if(B){
			assert_int_equal(PN(B), 2);
			assert_null(P(B, 0));
			assert_null(P(B, 1));
			cellfree(B);
		}
		dcellfree(Ac);
	}
}
/**
   Map the shared memory object fn read only and check its content against A
   using the bin format header parser.
*/
static void check_shm(const char *fn, const dcell *A, int iscellarr){
	int fd=shm_open(fn, O_RDONLY, 0600);
	assert_true(fd>=0);
	struct stat st;
	assert_int_equal(fstat(fd, &st), 0);
	void *p0=mmap(NULL, st.st_size, PROT_READ, MAP_SHARED, fd, 0);
	close(fd);
	assert_true(p0!=MAP_FAILED);
	char *p=(char *)p0;
	uint32_t magic;
	long nx, ny;
	const char *kw;
	if(iscellarr){
		mmap_read_header(&p, &magic, &nx, &ny, &kw);
		assert_true(iscell(&magic));
		assert_int_equal(nx, NX(A));
		assert_int_equal(ny, NY(A));
	}
	for(long ic=0; ic<PN(A); ic++){
		mmap_read_header(&p, &magic, &nx, &ny, &kw);
		const dmat *Ai=P(A, ic);
		assert_int_equal(nx*ny, PN(Ai));
		if(nx>0&&ny>0){
			assert_int_equal(magic, M_REAL);
			assert_int_equal(nx, NX(Ai));
			assert_int_equal(memcmp(p, P(Ai), nx*ny*sizeof(real)), 0);
			if(!iscellarr){
				assert_non_null(kw);
				assert_float_equal(search_keyword_num(kw, "dx"), 1, 0);
			}
		}
		p+=nx*ny*sizeof(real);
	}
	munmap(p0, st.st_size);
}
/**
   Arrays created with mmap in shared memory are stored in the bin format and
   can be read back with mmap or normal reading.
*/
static void fits_mmap(void **state){
	(void)state;
#if defined(__linux__)
	char fn[64];
	snprintf(fn, sizeof(fn), "/maos_cmocka_%d.bin", (int)getpid());
	rand_t rstat;
	seed_rand(&rstat, 6);
	dmat *A=dnew_mmap(20, 10, "dx=1", "%s", fn);
	assert_non_null(A);
	drandn(A, 1, &rstat);
	dcell *Ac=dcellnew(1, 1);
	P(Ac, 0)=dref(A);
	check_shm(fn, Ac, 0);
	dmat *B=dread_mmap("%s", fn);
	assert_non_null(B);
	assert_dmat_near(A, B, 0);
	assert_float_equal(search_keyword_num(B->keywords, "dx"), 1, 0);
	dfree(B);
	B=dread("%s", fn);
	assert_non_null(B);
	assert_dmat_near(A, B, 0);
	dfree(B);
	dcellfree(Ac);
	dfree(A);
	shm_unlink(fn);

	long nnx[6]={2, 3, 4, 5, 6, 0};
	long nny[6]={3, 4, 2, 5, 6, 3};
	dcell *C=dcellnew_mmap(2, 3, nnx, nny, NULL, "%s", fn);
	assert_non_null(C);
	assert_null(P(C, 5));//zero size block
	for(int ix=0; ix<PN(C); ix++){
		if(P(C, ix)){
			assert_int_equal(NX(P(C, ix)), nnx[ix]);
			assert_int_equal(NY(P(C, ix)), nny[ix]);
			assert_float_equal(dsumabs(P(C, ix)), 0, 0);//initialized to zero
			drandn(P(C, ix), 1, &rstat);
		}
	}
	check_shm(fn, C, 1);
	dcell *D=dcellread_mmap("%s", fn);
	assert_non_null(D);
	assert_int_equal(PN(D), PN(C));
	for(int ix=0; ix<PN(C); ix++){
		if(P(C, ix)) assert_dmat_near(P(C, ix), P(D, ix), 0);
	}
	dcellfree(D);
	dcellfree(C);
	shm_unlink(fn);
	//Anonymous memory when the name starts with -
	A=dnew_mmap(5, 5, NULL, "-");
	assert_non_null(A);
	assert_float_equal(dsumabs(A), 0, 0);
	dset(A, 2);
	assert_float_equal(dsum(A), 50, 0);
	dfree(A);
#else
	skip();
#endif
}

/*Raw arrays of each element type for fits_raw.*/
static const int8_t raw_i8[6]={-128, -1, 0, 1, 64, 127};
static const int16_t raw_i16[6]={-32768, -300, 0, 1, 300, 32767};
static const int32_t raw_i32[6]={INT32_MIN, -70000, 0, 1, 70000, INT32_MAX};
static const int64_t raw_i64[6]={INT64_MIN, -5000000000LL, 0, 1, 5000000000LL, INT64_MAX};
static const float raw_zmp[12]={1, -2, 3.5f, -4.5f, 5, 6, 7, 8, 9, 10, 11, 12};
static const double raw_dbl[6]={-1e300, -1.5, 0, 1e-300, 2.5, 1e300};
typedef struct{
	uint32_t magic;
	size_t size;
	const void *p;
	uint64_t nx, ny;
}raw_t;
static const raw_t raws[]={
	{M_INT16, 2, raw_i16, 3, 2},
	{M_INT32, 4, raw_i32, 2, 3},
	{M_INT8, 1, raw_i8, 2, 3},//not last: checks that the fits block padding after 1 byte data is skipped
	{M_INT64, 8, raw_i64, 6, 1},
	{M_ZMP, 8, raw_zmp, 3, 2},
	{M_DBL, 8, raw_dbl, 1, 6},
};
#define NRAW (int)(sizeof(raws)/sizeof(raws[0]))
static void *raw_writer(void *data){
	writer_t *w=(writer_t *)data;
	file_t *fp=zfdopen_as(w->fd, "w", w->suffix);
	for(int i=0; i<NRAW; i++){
		writearr(fp, 0, raws[i].size, raws[i].magic, i==1?"key=1":NULL, raws[i].p, raws[i].nx, raws[i].ny);
	}
	if(!zfisfits(fp)){
		uint32_t magic=M_SKIP;
		zfwrite(&magic, sizeof(uint32_t), 1, fp);//skipped by the reader
		magic=M_EOD;
		zfwrite(&magic, sizeof(uint32_t), 1, fp);
		magic=0x1234;//invalid
		zfwrite(&magic, sizeof(uint32_t), 1, fp);
		uint64_t dims[2]={1, 1};
		zfwrite(dims, sizeof(uint64_t), 2, fp);
	}
	zfclose(fp);
	return NULL;
}
/** writearr() and read_header()/zfread() keep integer, complex float, and double arrays in bin and fits formats. */
static void fits_raw(void **state){
	(void)state;
	assert_int_equal(writearr(NULL, 0, 8, M_DBL, NULL, raw_dbl, 6, 1), -1);
	for(int is=0; is<NSUFFIX; is++){
		int fds[2];
		assert_int_equal(pipe(fds), 0);
		writer_t w={fds[1], suffixes[is], NULL};
		pthread_t thread;
		pthread_create(&thread, NULL, raw_writer, &w);
		file_t *fp=zfdopen_as(fds[0], "r", suffixes[is]);
		assert_null(zfname(fp));
		const int isfits=zfisfits(fp);
		for(int i=0; i<NRAW; i++){
			header_t header={0, 0, 0, NULL};
			assert_int_equal(read_header(&header, fp), 0);
			assert_int_equal(header.id, raws[i].magic);
			assert_int_equal(header.nx, raws[i].nx);
			assert_int_equal(header.ny, raws[i].ny);
			if(i==1){
				assert_non_null(header.str);
				assert_non_null(strcasestr(header.str, "key"));
			}
			free(header.str);
			char buf[256]={0};
			//complex numbers are read as pairs of real numbers (see readvec)
			const int nc=magic_iscomplex(raws[i].magic)?2:1;
			assert_int_equal(zfread(buf, raws[i].size/nc, raws[i].nx*raws[i].ny*nc, fp), 0);
			assert_memory_equal(buf, raws[i].p, raws[i].size*raws[i].nx*raws[i].ny);
		}
		assert_int_equal(zferrno(fp), 0);
		if(!isfits){
			header_t header={0, 0, 0, NULL};
			assert_int_equal(read_header(&header, fp), 1);//end of data
			assert_int_not_equal(read_header(&header, fp), 0);//invalid magic
			assert_int_equal(header.id, 0);
			assert_int_equal(zferrno(fp), 2);
			assert_int_equal(read_header(&header, fp), -1);//previous error
			char c;
			assert_int_equal(zfread(&c, 1, 1, fp), -1);
		}
		//drain any remaining data so that the writer is not blocked.
		char buf[1024];
		if(fp) while(zfread(buf, 1, sizeof(buf), fp)==0);
		pthread_join(thread, NULL);
		zfclose(fp);
	}
	assert_int_equal(zferrno(NULL), -1);
	assert_null(zfname(NULL));
}
/**
   Read exactly n bytes from fd.
*/
static void read_full(int fd, unsigned char *buf, long n){
	while(n>0){
		ssize_t nr=read(fd, buf, n);
		assert_true(nr>0);
		buf+=nr; n-=nr;
	}
}
/** Convert big endian bytes to native order. */
static void from_big(void *out, const unsigned char *in, int size){
	unsigned char *po=(unsigned char *)out;
	for(int i=0; i<size; i++){
		po[i]=in[BIGENDIAN?i:(size-1-i)];
	}
}
/**
   Complex numbers are stored in fits files as big endian real followed by
   big endian imaginary part (as read by astropy/cfitsio), and round trip.
*/
static void fits_complex(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 5);
	zmat *Az=znew(3, 2);
	zrandn(Az, 1, &rstat);
	cmat *Ac=cnew(2, 3);
	crandn(Ac, 1, &rstat);
	for(int it=0; it<2; it++){
		int fds[2];
		assert_int_equal(pipe(fds), 0);
		file_t *fp=zfdopen_as(fds[1], "w", ".fits");
		writedata(fp, it==0?(cell*)Az:(cell*)Ac, 0);
		zfclose(fp);//data is small enough to fit in the pipe buffer
		unsigned char buf[2880*2];
		read_full(fds[0], buf, sizeof(buf));//one header block and one data block
		close(fds[0]);
		assert_memory_equal(buf, "SIMPLE  =", 9);
		if(it==0){
			float re, im;
			from_big(&re, buf+2880, 4);
			from_big(&im, buf+2884, 4);
			assert_true(re==crealf(P(Az, 0)));
			assert_true(im==cimagf(P(Az, 0)));
		}else{
			double re, im;
			from_big(&re, buf+2880, 8);
			from_big(&im, buf+2888, 8);
			assert_true(re==creal(P(Ac, 0)));
			assert_true(im==cimag(P(Ac, 0)));
		}
	}
	for(int is=0; is<NSUFFIX; is++){
		zmat *Bz=zmat_cast(roundtrip(Az, suffixes[is], M_ZMP));
		assert_non_null(Bz);
		assert_int_equal(PN(Bz), PN(Az));
		assert_memory_equal(P(Bz), P(Az), sizeof(fcomplex)*PN(Az));
		zfree(Bz);
		//converted on reading
		cmat *Bc=cmat_cast(roundtrip(Az, suffixes[is], M_CMP));
		assert_non_null(Bc);
		for(long i=0; i<PN(Az); i++){
			assert_true(P(Bc, i)==(dcomplex)P(Az, i));
		}
		cfree(Bc);
	}
	zfree(Az);
	cfree(Ac);
}
/**
   Asynchronous writing (as used by dnew_file and dcellnew_file) is not
   supported for fits files: the data is written in full when the array is
   freed.
*/
static void fits_async(void **state){
	(void)state;
	for(int it=0; it<2; it++){
		int fds[2];
		assert_int_equal(pipe(fds), 0);
		file_t *fp=zfdopen_as(fds[1], "w", ".fits");
		dmat *A=NULL;
		dcell *C=NULL;
		if(it==0){
			A=dnew(3, 4);
			A->fp=fp;
		}else{
			C=dcellnew_same(2, 1, 3, 4);
			C->fp=fp;
			A=P(C, 1);
		}
		writedata(fp, it==0?(cell*)A:(cell*)C, -1);//no async for fits
		assert_null(A->async);
		for(long iy=0; iy<NY(A); iy++){
			for(long ix=0; ix<NX(A); ix++){
				P(A, ix, iy)=ix+iy*10+1;
			}
			writecell_async(it==0?(cell*)A:(cell*)C, iy+1);//no-op
		}
		dmat *A0=ddup(A);
		if(it==0){
			dfree(A);//writes the data and closes fp
		}else{
			dcellfree(C);
		}
		file_t *fpr=zfdopen_as(fds[0], "r", ".fits");
		cell *B=readbin_fp(fpr, 0, -1);
		zfclose(fpr);
		assert_non_null(B);
		dmat *B1=NULL;
		if(it==0){
			B1=dmat_cast(B);
		}else{
			assert_true(iscell(B));
			assert_int_equal(PN(B), 2);
			B1=dmat_cast(P(B, 1));
			assert_int_equal(NX(P(B, 0)), 3);
		}
		assert_non_null(B1);
		assert_dmat_near(A0, B1, 0);
		cellfree(B);
		dfree(A0);
	}
}

UNITTEST_MAIN(
	cmocka_unit_test(fits_zfdopen),
	cmocka_unit_test(fits_dmat),
	cmocka_unit_test(fits_longkey),
	cmocka_unit_test(fits_types),
	cmocka_unit_test(fits_convert),
	cmocka_unit_test(fits_cell),
	cmocka_unit_test(fits_nested),
	cmocka_unit_test(fits_sparse),
	cmocka_unit_test(fits_empty),
	cmocka_unit_test(fits_mmap),
	cmocka_unit_test(fits_raw),
	cmocka_unit_test(fits_complex),
	cmocka_unit_test(fits_async),
)

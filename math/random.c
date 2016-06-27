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
#include <tgmath.h>
#include "random.h"
/**
   Routine to generate random numbers.
*/
/**
   Create a new random stream, seeded with seed.
*/
rand_t *new_rand(int seed){
    rand_t *out=(rand_t*)mycalloc(1,rand_t);
    seed_rand(out, seed);
    return out;
}
/**
   Save random stream to file.
 */
void writerand(rand_t *rstat, const char *format, ...){
    format2fn;
    FILE *fp=fopen(fn,"w");
    if(fwrite(rstat, 1, sizeof(rand_t), fp)!=1) 
	error("Error writing\n");
    fclose(fp);
}
/**
   Seed a random stream*/
void readrand(rand_t *rstat, const char *format,...){
    format2fn;
    if(rstat){
	FILE *fp=fopen(fn,"r");
	if(fread(rstat, 1, sizeof(rand_t), fp)!=1)
	    error("Error reading\n");
	fclose(fp);
    }
}


/* 
   extracted from mtwist.h/c
   The Web page on the Mersenne Twist algorithm is at:
   www.math.keio.ac.jp/~matumoto/emt.html 
   These functions were written by Geoffrey H. Kuenning, Claremont, CA.
   
   * This software is based on LGPL-ed code by Takuji Nishimura.  It has
   * also been heavily influenced by code written by Shawn Cokus, and
   * somewhat influenced by code written by Richard J. Wagner.  It is
   * therefore also distributed under the LGPL:
   *
   * This library is free software; you can redistribute it and/or
   * modify it under the terms of the GNU Library General Public License
   * as published by the Free Software Foundation; either version 2 of
   * the License, or (at your option) any later version.
   */
#define RECURRENCE_OFFSET 397
#define MATRIX_A	0x9908b0df
#define BIT_WIDTH	32		/* Work with 32-bit words */
#define UPPER_MASK	0x80000000	/* Most significant w-r bits */
#define LOWER_MASK	0x7fffffff	/* Least significant r bits */

#define COMBINE_BITS(x, y) \
			(((x) & UPPER_MASK) | ((y) & LOWER_MASK))
#define MATRIX_MULTIPLX(original, new) \
			((original) ^ ((new) >> 1) \
			  ^ matrix_decider[(new) & 0x1])
#define KNUTH_MULTIPLIER_OLD \
			69069

/*
 * Parameters of Knuth's PRNG (p. 106 of "The Art of Computer
 * Programming, Vol. 2, 3rd ed).
 */
#define KNUTH_MULTIPLIER_NEW 1812433253ul
#define KNUTH_SHIFT	30		/* Even on a 64-bit machine! */

/*
 * Default 32-bit random seed if mts_seed32 wasn't called
 */
#define DEFAULT_SEED32_OLD 4357
#define DEFAULT_SEED32_NEW 5489ul

/*
 * Where to get random numbers
 */
#define DEVRANDOM	"/dev/random"
#define DEVURANDOM	"/dev/urandom"

static mt_u32bit_t	matrix_decider[2] =
			  {0x0, MATRIX_A};
/**
   Mark state as initialized.
 */
static void mts_mark_initialized(mt_state* state /**< State vector to mark initialized */
				 ){
    state->initialized = 1;
}
/**
   Refresh the state for next set of random numbers.
 */
void mts_refresh(register mt_state* state /**< State for the PRNG */
		 ){
    register int	i;		/* Index into the state */
    register mt_u32bit_t*
			state_ptr;	/* Next place to get from state */
    register mt_u32bit_t
			value1;		/* Scratch val picked up from state */
    register mt_u32bit_t
			value2;		/* Scratch val picked up from state */

    /*
     * Start by making sure a random seed has been set.  If not, set
     * one.
     */
    if (!state->initialized)
	{
	mts_seed32(state, DEFAULT_SEED32_OLD);
	return;				/* Seed32 calls us recursively */
	}

    state_ptr = &state->statevec[MT_STATE_SIZE - 1];
    value1 = *state_ptr;
    for (i = (MT_STATE_SIZE - RECURRENCE_OFFSET) / 2;  --i >= 0;  )
	{
	state_ptr -= 2;
	value2 = state_ptr[1];
	value1 = COMBINE_BITS(value1, value2);
	state_ptr[2] =
	  MATRIX_MULTIPLX(state_ptr[-RECURRENCE_OFFSET + 2], value1);
	value1 = state_ptr[0];
	value2 = COMBINE_BITS(value2, value1);
	state_ptr[1] =
	  MATRIX_MULTIPLX(state_ptr[-RECURRENCE_OFFSET + 1], value2);
	}
    value2 = *--state_ptr;
    value1 = COMBINE_BITS(value1, value2);
    state_ptr[1] =
      MATRIX_MULTIPLX(state_ptr[-RECURRENCE_OFFSET + 1], value1);

    for (i = (RECURRENCE_OFFSET - 1) / 2;  --i >= 0;  )
	{
	state_ptr -= 2;
	value1 = state_ptr[1];
	value2 = COMBINE_BITS(value2, value1);
	state_ptr[2] =
	  MATRIX_MULTIPLX(state_ptr[MT_STATE_SIZE - RECURRENCE_OFFSET + 2],
	    value2);
	value2 = state_ptr[0];
	value1 = COMBINE_BITS(value1, value2);
	state_ptr[1] =
	  MATRIX_MULTIPLX(state_ptr[MT_STATE_SIZE - RECURRENCE_OFFSET + 1],
	    value1);
	}

    /*
     * The final entry in the table requires the "previous" value
     * to be gotten from the other end of the state vector, so it
     * must be handled specially.
     */
    value1 = COMBINE_BITS(value2, state->statevec[MT_STATE_SIZE - 1]);
    *state_ptr =
      MATRIX_MULTIPLX(state_ptr[MT_STATE_SIZE - RECURRENCE_OFFSET], value1);

    /*
     * Now that refresh is complete, reset the state pointer to allow more
     * pseudorandom values to be fetched from the state array.
     */
    state->stateptr = MT_STATE_SIZE;
    }

/**
   See the state vector.
 */
void mts_seed32(
    mt_state*		state,		/* State vector to initialize */
    unsigned long	seed)		/* 32-bit seed to start from */
    {
	seed=~seed;/*added by lianqiw. */
    int			i;		/* Loop index */

    if (seed == 0)
	seed = DEFAULT_SEED32_OLD;

    /*
     * Fill the state vector using Knuth's PRNG.  Be sure to mask down
     * to 32 bits in case we're running on a machine with 64-bit
     * longs.
     */
    state->statevec[MT_STATE_SIZE - 1] = seed & 0xffffffff;
    for (i = MT_STATE_SIZE - 2;  i >= 0;  i--)
        state->statevec[i] =
          (KNUTH_MULTIPLIER_OLD * state->statevec[i + 1]) & 0xffffffff;

    state->stateptr = MT_STATE_SIZE;
    mts_mark_initialized(state);

    /*
     * Matsumoto and Nishimura's implementation refreshes the PRNG
     * immediately after running the Knuth algorithm.  This is
     * probably a good thing, since Knuth's PRNG doesn't generate very
     * good numbers.
     */
    mts_refresh(state);
    }

/**
   Get a poisson random number.
*/
long randp(rand_t *rstat, double xm){
    double g, t;
    const long thres=200;/*use gaussian distribution instead*/
    double xu,xmu;
    long x;
    x=0;
    if(xm>thres){
	x=(long)round(xm+randn(rstat)*sqrt(xm));
    }else{
	while(xm>0){
	    xmu = xm > 12 ? 12 : xm;
	    xm-=xmu;
	    g=exp(-xmu);
	    xu=-1;
	    t=1;
	    while(t>g){
		xu++;
		t*=randu(rstat);
	    }
	    x+=xu;
	}
    }
    return x;
}

/*The following is ziggurate normal number generator */

#define LEVELS 256 /*must be 2^n*/
#define LEVELS_1 (LEVELS-1) /*LEVELS-1*/
#define RANGE 1. /*range of random numbber from random number genrator*/

static  double  xx[LEVELS];
static  double  ytab[LEVELS];
static  double  ktab[LEVELS];
static  double  wtab[LEVELS];
static  double  x1,x11;
static  double  v;

/**
   the target density */
INLINE double ff (double xr){
    return  exp(-xr*xr*0.5);
}
/**
   the inverse of the target density */
INLINE double f_inv (double y){
    return  sqrt(-2.*log(y));
}
/**
   try_r_value
*/
static double try_r_value (double r){
    int  i;
    v = r*ff(r) + exp(-0.5*r*r)/r;
    xx[LEVELS-1] = r;
    for (i=LEVELS-1; i>1; --i) {
	xx[i-1] = f_inv(v/xx[i]+ff(xx[i]));
    }
    return  xx[1]*(1-ff(xx[1]))/v;
}
/**
   Initialize the data to do norm distribution.
 */
static __attribute__((constructor)) void initialize(){
    double  a, b, aa, bb, r;
    double  q;
    double RANGE_R=1./RANGE;
    int  i;

    a=0;
    b=10;
    do {
	aa=a, bb=b;
	r=.5*(a+b);
	q=try_r_value(r);
	if (q>1) { 
	    b = r;
	} else {
	    a = r;
	}
    } while (aa<r && r<bb);
    xx[0] = 0;
    x1=r;
    x11=1./x1;
    for(i=0; i<LEVELS; i++){
	ytab[i]=ff(xx[i]);
    }
    
    for(i=0; i<LEVELS-1; i++){
	ktab[i]=(xx[i]/xx[i+1])*RANGE;
    }
    ktab[LEVELS-1]=x1*ff(x1)*RANGE/v;

    for(i=0; i<LEVELS-1; i++){
	wtab[i]=xx[i+1]*RANGE_R;
    }
    wtab[LEVELS-1]=v*RANGE_R/ff(x1);
}
/**
   Get normal distributed random number.
*/
double randn(rand_t *rstat){
    double x=0;
    double U;
    unsigned int l;
    int i;
    while(1){
	U = 2*randu(rstat)-1;
	l = lrand(rstat);
	i = l & LEVELS_1;/*8 bit for addressing*/
	x = U*wtab[i];
	if(fabs(U)<ktab[i]) break;
	if(i==LEVELS_1){/*tail.*/
	    double x0, y;
	    do{
		x0=log(randu(rstat))*x11;
		y=-log(randu(rstat));
	    }while(y+y<x0*x0);
	    x=x>0?(x1-x0):(x0-x1);
	    break;
	}else{
	    double yy0, yy1, y;
	    yy0=ytab[i];
	    yy1=ytab[i+1];
	    y = yy1+(yy0-yy1)*randu(rstat);
	    if(y<ff(x)) break;
	}
    }
    return x;
}

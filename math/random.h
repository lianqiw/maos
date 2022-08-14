/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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


#ifndef __AOS_RANDOM_H
#define __AOS_RANDOM_H
#include "../sys/sys.h"
/*
   Random number generators.
*/
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


#ifndef MT_MACHINE_BITS
#include <limits.h>
#if LONG_MAX == 2147483647 /*changed to LONG by lianqiw */
#define MT_MACHINE_BITS	32
#else /* INT_MAX */
#define MT_MACHINE_BITS	64
#endif /* INT_MAX */
#endif /* MT_MACHINE_BITS */

/*
 * Define an unsigned type that is guaranteed to be 32 bits wide.
 */
#if MT_MACHINE_BITS == 32
typedef unsigned long	mt_u32bit_t;
#else /* MT_MACHINE_BITS */
typedef unsigned int	mt_u32bit_t;
#endif /* MT_MACHINE_BITS */
/*
 * The following value is a fundamental parameter of the algorithm.
 * It was found experimentally using methods described in Matsumoto
 * and Nishimura's paper.  It is exceedingly magic; don't change it.
 */
#define MT_STATE_SIZE	624		/* Size of the MT state vector */

/*
 * Internal state for an MT RNG.  The user can keep multiple mt_state
 * structures around as a way of generating multiple streams of random
 * numbers.
 *
 * In Matsumoto and Nishimura's original paper, the state vector was
 * processed in a forward direction.  I have reversed the state vector
 * in this implementation.  The reason for the reversal is that it
 * allows the critical path to use a test against zero instead of a
 * test against 624 to detect the need to refresh the state.  on most
 * machines, testing against zero is slightly faster.  It also means
 * that a state that has been set to all zeros will be correctly
 * detected as needing initialization; this means that setting a state
 * vector to zero (either with memset or by statically allocating it)
 * will cause the RNG to operate properly.
 */
/**
   Contains state of a random stream.
 */
typedef struct{
	mt_u32bit_t		statevec[MT_STATE_SIZE]; /**< Vector holding current state */
	int			stateptr;	/**< Next state entry to be used */
	int			initialized;	/**< NZ if state was initialized */
}mt_state;
#define MT_TEMPERING_MASK_B 0x9d2c5680
#define MT_TEMPERING_MASK_C 0xefc60000
#define MT_TEMPERING_SHIFT_U(y) (y >> 11)
#define MT_TEMPERING_SHIFT_S(y) (y << 7)
#define MT_TEMPERING_SHIFT_T(y) (y << 15)
#define MT_TEMPERING_SHIFT_L(y) (y >> 18)

/*
 * Macros to do the tempering.  MT_PRE_TEMPER does all but the last step;
 * it's useful for situations where the final step can be incorporated
 * into a return statement.  MT_FINAL_TEMPER does that final step (not as
 * an assignment).  MT_TEMPER does the entire process.  Note that
 * MT_PRE_TEMPER and MT_TEMPER both modify their arguments.
 */
#define MT_PRE_TEMPER(value)						\
    do									\
	{								\
	value ^= MT_TEMPERING_SHIFT_U(value);				\
	value ^= MT_TEMPERING_SHIFT_S(value) & MT_TEMPERING_MASK_B;	\
	value ^= MT_TEMPERING_SHIFT_T(value) & MT_TEMPERING_MASK_C;	\
	}								\
	while (0)
#define MT_FINAL_TEMPER(value) \
			((value) ^ MT_TEMPERING_SHIFT_L(value))
#define MT_TEMPER(value)						\
    do									\
	{								\
	value ^= MT_TEMPERING_SHIFT_U(value);				\
	value ^= MT_TEMPERING_SHIFT_S(value) & MT_TEMPERING_MASK_B;	\
	value ^= MT_TEMPERING_SHIFT_T(value) & MT_TEMPERING_MASK_C;	\
	value ^= MT_TEMPERING_SHIFT_L(value);				\
	}								\
	while (0)
void mts_refresh(mt_state* state);
#define mt_32_to_double (1./4294967296.)
/* Multiplier to convert long to dbl [0,1)*/
#define mt_64_to_double (1./18446744073709551616.)
/* Mult'r to cvt long long to dbl */
/**<32 bit val */
static inline unsigned long mts_lrand(mt_state* state)		/* State for the PRNG */
{
	unsigned long random_value;	/* Pseudorandom value generated */

	if(state->stateptr<=0)
		mts_refresh(state);

	random_value=state->statevec[--state->stateptr];
	MT_PRE_TEMPER(random_value);
	return MT_FINAL_TEMPER(random_value);
}
/**<32 bit precision double. [0,1) */
static inline double mts_drand(mt_state* state)		/* State for the PRNG */
{
	unsigned long random_value;	/* Pseudorandom value generated */

	if(state->stateptr<=0)
		mts_refresh(state);

	random_value=state->statevec[--state->stateptr];
	MT_TEMPER(random_value);

	return random_value*mt_32_to_double;
}
/**<64 bit precision double */
static inline double mts_ldrand(mt_state* state)		/* State for the PRNG */
{
#if MT_MACHINE_BITS == 64
	unsigned long long	final_value;	/* Final (integer) value */
#endif /* MT_MACHINE_BITS */
	unsigned long random_value_1;	/* 1st pseudorandom value generated */
	unsigned long random_value_2;	/* 2nd pseudorandom value generated */

	/*
	 * For maximum speed, we'll handle the two overflow cases
	 * together.  That will save us one test in the common case, at
	 * the expense of an extra one in the overflow case.
	 */
	if(--state->stateptr<=0){
		if(state->stateptr<0){
			mts_refresh(state);
			random_value_1=state->statevec[--state->stateptr];
		} else{
			random_value_1=state->statevec[state->stateptr];
			mts_refresh(state);
		}
	} else
		random_value_1=state->statevec[--state->stateptr];

	MT_TEMPER(random_value_1);

	random_value_2=state->statevec[--state->stateptr];
	MT_TEMPER(random_value_2);

#if MT_MACHINE_BITS == 64
	final_value=((unsigned long long) random_value_1<<32)
		|(unsigned long long) random_value_2;
	return final_value*mt_64_to_double;
#else /* MT_MACHINE_BITS */
	return random_value_1*mt_32_to_double+random_value_2*mt_64_to_double;
#endif /* MT_MACHINE_BITS */
}

void mts_seed32(mt_state* state, unsigned long	seed);
typedef mt_state  rand_t;
#define seed_rand mts_seed32
#define lrand     mts_lrand /*rand integer */
#define randu     mts_drand /*rand double [0,1) */
#define lrandu    mts_ldrand/*rand double with 64 bit precision within [0,1) */
double randn(rand_t* rstat);/*normal distribution with 1*/
long   randp(rand_t* rstat, double xm);
#endif

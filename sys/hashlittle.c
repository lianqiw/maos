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
#include <errno.h>
#include <unistd.h>
#include <string.h>
#include "hashlittle.h"

/*
 * My best guess at if you are big-endian or little-endian.  This may
 * need adjustment.
 */
#if (defined(__BYTE_ORDER) && defined(__LITTLE_ENDIAN) && \
     __BYTE_ORDER == __LITTLE_ENDIAN) || \
    (defined(i386) || defined(__i386__) || defined(__i486__) || \
     defined(__i586__) || defined(__i686__) || defined(vax) || defined(MIPSEL))
# define HASH_LITTLE_ENDIAN 1
# define HASH_BIG_ENDIAN 0
#elif (defined(__BYTE_ORDER) && defined(__BIG_ENDIAN) && \
       __BYTE_ORDER == __BIG_ENDIAN) || \
      (defined(sparc) || defined(POWERPC) || defined(mc68000) || defined(sel))
# define HASH_LITTLE_ENDIAN 0
# define HASH_BIG_ENDIAN 1
#else
# define HASH_LITTLE_ENDIAN 0
# define HASH_BIG_ENDIAN 0
#endif

#define hashsize(n) ((uint32_t)1<<(n))
#define hashmask(n) (hashsize(n)-1)
#define rot(x,k) (((x)<<(k)) | ((x)>>(32-(k))))

#define mix(a,b,c) \
{ \
  a -= c;  a ^= rot(c, 4);  c += b; \
  b -= a;  b ^= rot(a, 6);  a += c; \
  c -= b;  c ^= rot(b, 8);  b += a; \
  a -= c;  a ^= rot(c,16);  c += b; \
  b -= a;  b ^= rot(a,19);  a += c; \
  c -= b;  c ^= rot(b, 4);  b += a; \
}
#define final(a,b,c) \
{ \
  c ^= b; c -= rot(b,14); \
  a ^= c; a -= rot(c,11); \
  b ^= a; b -= rot(a,25); \
  c ^= b; c -= rot(b,16); \
  a ^= c; a -= rot(c,4);  \
  b ^= a; b -= rot(a,14); \
  c ^= b; c -= rot(b,24); \
}
/**
   hashlittle.c contains a hash function used to hash a chunk of data and
   return a hash number that can be used to identify the data. Apply with care
   to double floating numbers because of round off error. */
/**
   Hash the data and return a hash number. For long data only a subset of the
   words (stride length/10000) is sampled.

   @param[in] key The data.
   @param[in] length The length of the data in bytes.
   @param[in] initval The initial hash number.
   @return The hash number.
*/
uint32_t hashlittle(const void* key, size_t length, uint32_t initval){

	uint32_t a, b, c;                                          /* internal state */
	size_t dstep=(length/10000);/*jump */
	dstep=dstep<1?1:dstep;

	size_t dstep12=dstep*12;

	/* Set up the internal state */
	a=b=c=0xdeadbeef+((uint32_t)length)+initval;


	const uint8_t* k=(const uint8_t*)key;
	uint32_t w[3];/*memcpy avoids misaligned reads */

	/*------ all but last block: affect 32 bits of (a,b,c) */
	while(length>dstep12){
		memcpy(w, k, sizeof(w));
		a+=w[0];
		b+=w[1];
		c+=w[2];
		mix(a, b, c);
		length-=dstep12;
		k+=dstep12;
	}

	/*----------------------------- handle the last (probably partial) block */
	if(length==0) return c; /* zero length strings require no mixing */
	if(length<=12){/*zero padding gives the same result as masking the words on little endian machines.*/
		w[0]=w[1]=w[2]=0;
		memcpy(w, k, length);
		a+=w[0];
		b+=w[1];
		c+=w[2];
	}
	final(a, b, c);
	return c;
}

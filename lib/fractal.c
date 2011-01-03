#include "fractal.h"

/**
   \file fractal.c

   Implementation of the fractal operation for atmospheric turbulence screen
   generation and reconstruction.

*/


/**
   Apply forward fractal operator recursively.

   To confirm that fractal_trans is right, compare u'Av with v'A'u, where ' is transpose. 
   u,v are two vectors (reshape to 2-d arrays when applying A, or A')
*/
#define FRACTAL fractal
#define INVERSE 0
#define TRANSPOSE 0
#include "fractal_do.c"
#undef FRACTAL
#undef INVERSE
#undef TRANSPOSE

#define FRACTAL fractal_inv
#define INVERSE 1
#define TRANSPOSE 0
#include "fractal_do.c"
#undef FRACTAL
#undef INVERSE
#undef TRANSPOSE

#define FRACTAL fractal_trans
#define INVERSE 0
#define TRANSPOSE 1
#include "fractal_do.c"
#undef FRACTAL
#undef INVERSE
#undef TRANSPOSE


#define FRACTAL fractal_inv_trans
#define INVERSE 1
#define TRANSPOSE 1
#include "fractal_do.c"
#undef FRACTAL
#undef INVERSE
#undef TRANSPOSE


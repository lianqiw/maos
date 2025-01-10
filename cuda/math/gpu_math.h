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

#ifndef AOS_CUDA_GPU_MATH_H
#define AOS_CUDA_GPU_MATH_H
#define AOS_CUDA_H
/**
 * \file gpu_math.h
 * Routines that can be used by CPU.
 * */
#if defined(__cplusplus) && !USE_CPP
extern "C"{
#endif
#include "../../lib/aos.h"
void gpu_dgemm_test();
void gpu_dsvd_test();
void gpu_ext_assign();
#if defined(__cplusplus) && !USE_CPP
}
#endif
#endif

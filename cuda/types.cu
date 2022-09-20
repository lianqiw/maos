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

#include "types.h"

//Only specialize valid convertors.
//cpu code defines real, comp, dmat, cmat for double or signle precision accourding to CPU_SINGLE
template <>
NumArray<real, Cpu>::NumArray(const dmat *A){
	if(PN(A)){
		NumArray<real, Cpu> temp(NX(A), NY(A), (real *)P(A), 0);//weak reference
		*this=temp;//move data over.
	}
}
template <>
NumArray<comp, Cpu>::NumArray(const cmat *A){
	if(PN(A)){
		NumArray<comp, Cpu> temp(NX(A), NY(A), (comp*)P(A), 0);//weak reference
		*this=temp;//move data over.
	}
}
//Single precision
template <>
NumArray<float, Cpu>::NumArray(const smat *A){
	if(PN(A)){
		NumArray<float, Cpu> temp(NX(A), NY(A), (float *)P(A), 0);//weak reference
		*this=temp;//move data over.
	}
}
template <>
NumArray<fcomplex, Cpu>::NumArray(const zmat *A){
	if(PN(A)){
		NumArray<fcomplex, Cpu> temp(NX(A), NY(A), (fcomplex *)P(A), 0);//weak reference
		*this=temp;//move data over.
	}
}
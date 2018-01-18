/*
  Copyright 2009-2018 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "cucmat.h"
#include "utils.h"

void cucscale(cucmat &in, Real alpha, cudaStream_t stream){
    if(!in) return;
    if(alpha==0){
	cuzero(in, stream);
    }else if(Z(fabs)(alpha-1.f)>EPS){
	int n=in.Nx()*in.Ny();
	scale_do<<<DIM(n,256), 0, stream>>>(in.P(), n, alpha); 
    }
}

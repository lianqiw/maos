#ifndef AOS_MAOS_MVM_CLIENT_H
#define AOS_MAOS_MVM_CLIENT_H
/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "maos.h"
#define USE_SHORT 1
#if USE_SHORT
#define GSCALE 2062650000. //convert grad to int
#define GSCALE1 4.848132257047973e-10 //convert back
#define ASCALE 1.e9 //convert DM commands to int
#define ASCALE1 1.e-9 //convert int to DM commands
typedef short GTYPE; //type used for gradients
typedef short ATYPE; //type used for actuator commands
#else
typedef float GTYPE;
typedef float ATYPE;
#define GSCALE 1
#define GSCALE1 1
#define ASCALE 1
#define ASCALE1 1
#endif
#define N_CMD 4
enum{
    GPU_MVM_ZERO,
    GPU_MVM_M,
    GPU_MVM_G,
    GPU_MVM_A,
};

void mvm_client_init(const char *host, int port);
void mvm_client_send_m(const PARMS_T *parms, dcell *mvm);
void mvm_client_recon(const PARMS_T *parms, dcell *dm, dcell *grad);
void mvm_client_close(void);
#endif

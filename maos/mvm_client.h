/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "common.h"

/*
  \file mvm_client.h

  For testing MVM across two machines.
*/

#ifndef MVM_CLIENT_H
#define MVM_CLIENT_H
#define USE_SHORT 0
#if USE_SHORT
#define GSCALE 2062650000. //convert grad to int
#define GSCALE1 4.848132257047973e-10 //convert back
#define ASCALE 1.e9 //convert DM commands to int
#define ASCALE1 1.e-9 //convert int to DM commands
typedef short GReal; //type used for gradients
typedef short AReal; //type used for actuator commands
#else
#if CUDA_DOUBLE
typedef real GReal;
typedef real AReal;
#else
typedef float GReal;
typedef float AReal;
#endif
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

void mvm_client_init(const char *host, int port, dmat *mvm, int ngpu);
void mvm_client_send_m(dmat *mvm, int ngpu);
void mvm_client_recon(int mvmsize, dcell *dm, dcell *grad);
void mvm_client_close(void);
#endif

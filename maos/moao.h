#ifndef MAOS_MOAO_H
#define MAOS_MOAO_H
#include "types.h"
void free_recon_moao(RECON_T *recon, const PARMS_T *parms);
void setup_recon_moao(RECON_T *recon, const PARMS_T *parms);
void moao_recon(SIM_T *simu);
#endif

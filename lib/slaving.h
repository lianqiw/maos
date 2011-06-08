#ifndef AOS_LIB_SLAVING_H
#define AOS_LIB_SLAVING_H
#include "type.h"
#include "imat.h"
spcell *slaving(loc_t **aloc, spcell *HA, dmat *W1, dcell *NW, icell *actstuck, icell *actfloat, double thres, double scl);
void act_stuck(loc_t **aloc, spcell *HA, icell *stuck);
void act_float(loc_t **aloc, spcell **HA, icell *actfloat);
#endif

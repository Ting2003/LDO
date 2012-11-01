#ifndef _COMMON_H
#define _COMMON_H
#include "ldo.h"
// stores common functions
bool node_in_wspace(double x, double y, WSPACE *wspace);
bool ldo_in_wspace(LDO *ldo, WSPACE *wspace);
#endif

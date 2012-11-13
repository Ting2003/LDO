#ifndef _COMMON_H
#define _COMMON_H
#include "ldo.h"
#include "math.h"
#include <queue>
// stores common functions
bool node_in_wspace(double x, double y, MODULE *module);
bool ldo_in_wspace(LDO *ldo, MODULE *module);
bool ldo_in_wspace_trial(double ref_dist, double ref_x, double ref_y, double &x2, double &y2, LDO &ldo, vector<MODULE*> wspacelist);
#endif

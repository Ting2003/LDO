#ifndef _CHIP_H
#define _CHIP_H

#include "ldo.h"
#include "circuit.h"
// chip includes circuit of power grid and LDO device
class Chip{
private:
public:
	Chip();
	friend class Parser;
	vector<Circuit *>cktlist;
	vector<LDO*> ldolist;
};

#endif

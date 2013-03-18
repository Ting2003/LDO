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
	// stores the circuit from different domain
	// including a global and several local
	vector<Circuit *>cktlist;
	vector<LDO*> ldolist;
	vector<MODULE*> wspacelist;

	//void extract_rec();
	//void extract_rec_single(size_t k);
};

#endif

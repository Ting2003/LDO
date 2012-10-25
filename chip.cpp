#include "chip.h"

Chip::Chip(){
	cktlist.clear();
	ldolist.clear();
}

// add extra node to LDO node A
void Chip::add_net_LDO(){
	LDO *ldo_ptr;
	Node *na, *nb;
	for(size_t i=0;i<ldolist.size();i++){
		ldo_ptr = ldolist[i];	
	}		
}

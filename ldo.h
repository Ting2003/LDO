#ifndef _LDO_H
#define _LDO_H
#include "node.h"

class LDO{
private:
public:
	LDO();	
	Node *A, *B, *C, *D, *E;
	Node *in; // the input pin
	Node *out; // the output pin
	double width; // width of the ldo
	double height;// height of the ldo
};

#endif

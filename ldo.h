#ifndef _LDO_H
#define _LDO_H
#include "node.h"

class LDO{
private:
public:
	LDO();	
	Node *A;
	Node *in; // the input pin
	Node *out; // the output pin
	double width; // width of the ldo
	double height;// height of the ldo
	int degree;// 0, 90, 180, 270
};

#endif

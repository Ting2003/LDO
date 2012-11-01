#ifndef _LDO_H
#define _LDO_H
#include "node.h"

class LDO{
private:
public:
	LDO();
	string name;	
	Node*A;
	// 4 corners
	vector<Point*> node;
	int width;
	int height;
};

class WSPACE{
private:
public:
	WSPACE();
	string name;
	vector<Point*> node;
};

#endif

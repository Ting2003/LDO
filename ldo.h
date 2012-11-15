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
	// first node is the pin for VDD
	vector<Point*> node;
	int width;
	int height;
	// stores the output voltage
	double voltage;
	// current extracted from this LDO
	double current; 
};

// current block in polygon format
class MODULE{
private:
public:
	MODULE();
	string name;
	vector<Point*> node;
	// vector<int> LDO_id; // index of LDO within wspace
};

#endif

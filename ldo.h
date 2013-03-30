#ifndef _LDO_H
#define _LDO_H
#include "node.h"

class LDO{
private:
public:
	LDO();
	string name;	
	Node * A; // output node
	Node * nd_in; // input node
	// 4 corners
	// first node is the pin for VDD
	//vector<Point*> node;
	int width;
	int height;
	// stores the input voltage
	//double vin;
	// stores the output voltage
	double voltage;
	double voltage_old; // stores the old vol for A
	// current extracted from this LDO
	double current;
	double current_old;
};

// current block in polygon format
class MODULE{
private:
public:
	MODULE();
	string name;
	//Point * node;
	int width;
	int height;
	vector<Point*> node;
	// vector<int> LDO_id; // index of LDO within wspace
};

#endif

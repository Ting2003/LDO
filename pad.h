#ifndef __PAD_H__
#define __PAD_H__

#include <string>
#include <algorithm>
#include "global.h"
#include "point.h"
#include "node.h"
#include <vector>
#include <iostream>
#include <map>
//#include "block.h"
using namespace std;

class Pad{
public:
	Pad();
	~Pad();
	
	Node * node;
	// stores the node for LDO geo output node
	Node * nd_out_LDO;
	// record the controlled nodes and weighting
	map<Node*, double> control_nodes;
	// stores the neighboring pads for a pad
	vector<Pad*> nbrs;
	// sorted ref drop values of control nodes
	double ref_vol; // stores the middle voltage value
	vector<double> drop_vec;
	double newx;
	double newy;
	bool visit_flag;
	bool fix_flag;
	double data; // stores the maximum diff
	double ratio;
	bool violate_flag;
	// extracted current of each time step
	vector<double> current;
};
#endif

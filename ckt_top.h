// ----------------------------------------------------------------//
// Filename : ckt_top.h
// Author : Ting Yu <tingyu1@illinois.edu>
//
// declaration of sub_ckt class
// ----------------------------------------------------------------//

#ifndef __SUB_CKT_H__
#define __SUB_CKT_H__
#include <fstream>
#include "triplet.h"
#include "global.h"
#include "vec.h"
#include "net.h"
#include "util.h"
#include "cholmod.h"
#include "circuit.h"
using namespace std;

class CKT_TOP{
public:
	// stores the global and local circuit of each type circuit
	Circuit* ckt1;
	Circuit* ckt2;
        CKT_TOP(string name);
	~CKT_TOP();

	string name;
	vector<Net *> boundary_net;
	// stores interface global nodes
	vector<Node*> node_global;
	// stores all current values of interface global nodes
	vector<double> cur_global;
	// stores interface local nodes
	vector<Node*> node_local;
	// stores all voltage values of interface local nodes
	vector<double> vol_local;
};
#endif

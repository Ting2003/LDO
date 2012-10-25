#ifndef __TRANSIENT__H__
#define __TRANSIENT__H__

#include <iostream>
#include <string.h>
#include <vector>
#include "node.h"
using namespace std;

// define class for print out
class Node_TR_PRINT{
private:
	bool LDO_flag;
public:
	Node_TR_PRINT();
	~Node_TR_PRINT();
	string name;
	Node * node;
	int flag;
	vector<double> value;
	
	void enable_LDO();
	void disable_LDO();
	bool isLDO();
};

class Tran{
  public:
	Tran();
	~Tran();
	double step_t;
	double tot_t;
	int length;
	vector<Node_TR_PRINT> nodes;
	int isTran;
	void print_tr_nodes();
};

#endif


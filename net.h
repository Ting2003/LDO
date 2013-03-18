#ifndef __NET_H__
#define __NET_H__

#include <string>
#include "global.h"
#include "point.h"
using namespace std;

class Node;

class Net{
public:
	//Net(NET_TYPE type, string name, double value, Node * a, Node * b);
	Net(NET_TYPE type, double value, Node * a, Node * b);
	friend ostream & operator << (ostream & os, const Net & net);

	NET_TYPE type;	// identify the type
	//string name;
	double value; // dc value
	Node* ab[2];	// two connection node
	int id;		// index of each net
	bool flag_global; // see if the net is global
	// store the pulse parameter for current net
	double V1, V2, TD, Tr, Tf, PW, Period;

	void enableGlobal();
	bool is_global() const;
};

inline void Net::enableGlobal() {flag_global = true;}
inline bool Net::is_global() const{return flag_global;}
#endif

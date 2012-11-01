#ifndef _LDO_H
#define _LDO_H
#include "node.h"

class LDO{
private:
public:
	LDO();
	string name;	
	Node *A;
	Node *in; // the input pin
	Node *out; // the output pin
	double width; // width of the ldo
	double height;// height of the ldo
	int degree;// 0, 90, 180, 270
};

// white space class
// each white space is a rectangle
class RECTANGLE{
private:
public:
	RECTANGLE();
	// record 4 boundary lines
	string name;
	int xl, xr, yb, yt;
};

class WSPACE{
private:
public:
	WSPACE();
	string name;
	vector<Point*> node;
	vector<RECTANGLE*> rectangle;
};

#endif

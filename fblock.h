// ----------------------------------------------------------------//
// Filename : fblock.h
// Author : Ting Yu <tingyu1@illinois.edu>
//
// declaration of physical block class
// ----------------------------------------------------------------//

#ifndef __FBLOCK_H__
#define __FBLOCK_H__
#include <fstream>
#include "triplet.h"
#include "global.h"
#include "point.h"
#include "vec.h"
#include "util.h"

using namespace std;

// define block in floorplanning stage
class FBlock{
public:
	FBlock();
private:
	string name;
	// 4 corners in counter clockwise
	Point A, B, C, D;
	Point E; // power pin
	// rotation degree in floorplanning
	// degree = {0, 90, 180, 270}
	int degree; 
};

#endif


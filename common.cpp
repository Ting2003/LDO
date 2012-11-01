#include "common.h"
// judge whether a point is within white space or not
bool node_in_wspace(double x0, double y0, WSPACE *wspace){
	bool flag = false;
	int degree = 0;
	double x_new, y_new;
	double x_old, y_old;
	x_old = wspace->node[0]->x;	
	y_old = wspace->node[0]->y;

	// case A. node is exactly a corner of wspace
	for(size_t i=0;i<wspace->node.size();i++){
		if(x0 == wspace->node[i]->x && 
		   y0 == wspace->node[i]->y)
			return true;
	}
	// case B. node is on a boundary of wspace
	for(size_t i=1;i<wspace->node.size();i++){
		x_new = wspace->node[i]->x;	
		y_new = wspace->node[i]->y;
		// on horizontal line
		if(y_old == y_new){
			double dx0 = x0 - x_old;
			double dx1 = x0 - x_new;
			if(dx0 * dx1 < 0)
				return true;
		}// on vertical line
		else if(x_old == x_new){
			double dy0 = y0 - y_old;
			double dy1 = y0 - y_new;
			if(dy0 * dy1 < 0)
				return true;
		}
		
	}
	// case C. node is not A. and B.	
	for(size_t i=1;i<wspace->node.size();i++){
		x_new = wspace->node[i]->x;	
		y_new = wspace->node[i]->y;
		// 1 
		if(x_old > x0 && y_old > y0){
			// 1 -->2
			if(x_new < x0 && y_new > y0)	
				degree += 90;	
			
			// 1 -->4
			else if(x_new > x0 && y_new < y0)
				degree -= 90;
		}// 2
		else if(x_old < x0 && y_old > y0){
			// 2 --> 3
			if(x_new <x0 && y_new <y0)
				degree += 90;
			// 2 -->1
			if(x_new > x0 && y_new > y0)
				degree -= 90;	
		}// 3
		else if(x_old <x0 && y_old <x0){
			// 3 -->4
			if(x_new > x0 && y_new < y0)
				degree += 90;
			// 3 -->2
			if(x_new < x0 && y_new > y0)
				degree -= 90;
		}// 4
		else if(x_old >x0 && y_old <y0){
			// 4 -->1
			if(x_new > x0 && y_new >y0)
				degree += 90;
			// 4 -->3
			if(x_new <x0 && y_new <y0)
				degree -= 90;
		}
		x_old = x_new;
		y_old = y_new;
		
	}
	// if(wspace->name == "W11")
	// clog<<"degree: "<<degree<<endl;
	if(degree == 360)
		flag = true;
	else if(degree == 0)
		flag = false;
	return flag;
}

bool ldo_in_wspace(LDO *ldo, WSPACE *wspace){
	Point *pt;
	double x, y;
	bool node_flag;
	bool ldo_flag = true;
	for(size_t k=0;k<ldo->node.size();k++){
		x = ldo->node[k]->x;
		y = ldo->node[k]->y;
		node_flag = node_in_wspace(x, y, wspace);
		//if(wspace->name == "W11")
		//clog<<"wspace, x, y, flag: "<<wspace->name<<" "<<x<<" "<<y<<" "<<node_flag<<endl;
		if(!node_flag){
			ldo_flag = false;
			break;
		}
	}
	return ldo_flag;
}

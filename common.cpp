#include "common.h"
// judge whether a point is within white space or not
bool node_in_ldo(double x0, double y0, LDO *wspace){
	bool flag = false;
	int degree = 0;
	double x_new, y_new;
	double x_old, y_old;

	// case A. node is exactly a corner of wspace
	for(size_t i=0;i<wspace->node.size();i++){
		if(x0 == wspace->node[i]->x && 
		   y0 == wspace->node[i]->y){
			return true;
		}
	}
	// case B. node is on a boundary of wspace
	for(size_t i=1;i<wspace->node.size();i++){
		x_old = wspace->node[i-1]->x;	
		y_old = wspace->node[i-1]->y;
		x_new = wspace->node[i]->x;	
		y_new = wspace->node[i]->y;
		//clog<<"x_old, y_old, x_new, y_new: "<<
			//x_old<<" "<<y_old<<" "<<x_new<<"	"<<y_new<<endl;
		// on horizontal line
		if(y_old == y_new && y0 == y_old){
			double dx0 = x0 - x_old;
			double dx1 = x0 - x_new;
			if(dx0 * dx1 < 0){
				//clog<<"on edge. "<<endl;
				return true;
			}
		}// on vertical line
		else if(x_old == x_new && x0 == x_old){
			double dy0 = y0 - y_old;
			double dy1 = y0 - y_new;
			if(dy0 * dy1 < 0){
				//clog<<"on edge. "<<endl;
				return true;
			}
		}
		
	}
	// case C. node is not A. and B.	
	for(size_t i=1;i<wspace->node.size();i++){
		x_old = wspace->node[i-1]->x;	
		y_old = wspace->node[i-1]->y;
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
	// clog<<"wspace, degree: "<<wspace->name<<" "<<degree<<endl;
	if(degree == 360)
		flag = true;
	else if(degree == 0)
		flag = false;
	return flag;
}

// judge whether a point is within white space or not
bool node_in_wspace(double x0, double y0, MODULE *wspace){
	bool flag = false;
	int degree = 0;
	double x_new, y_new;
	double x_old, y_old;

	// case A. node is exactly a corner of wspace
	for(size_t i=0;i<wspace->node.size();i++){
		if(x0 == wspace->node[i]->x && 
		   y0 == wspace->node[i]->y){
			return true;
		}
	}
	// case B. node is on a boundary of wspace
	for(size_t i=1;i<wspace->node.size();i++){
		x_old = wspace->node[i-1]->x;	
		y_old = wspace->node[i-1]->y;
		x_new = wspace->node[i]->x;	
		y_new = wspace->node[i]->y;
		//clog<<"x_old, y_old, x_new, y_new: "<<
			//x_old<<" "<<y_old<<" "<<x_new<<"	"<<y_new<<endl;
		// on horizontal line
		if(y_old == y_new && y0 == y_old){
			double dx0 = x0 - x_old;
			double dx1 = x0 - x_new;
			if(dx0 * dx1 < 0){
				//clog<<"on edge. "<<endl;
				return true;
			}
		}// on vertical line
		else if(x_old == x_new && x0 == x_old){
			double dy0 = y0 - y_old;
			double dy1 = y0 - y_new;
			if(dy0 * dy1 < 0){
				//clog<<"on edge. "<<endl;
				return true;
			}
		}
		
	}
	
	// case C. node is not A. and B.	
	for(size_t i=1;i<wspace->node.size();i++){
		x_old = wspace->node[i-1]->x;	
		y_old = wspace->node[i-1]->y;
		x_new = wspace->node[i]->x;	
		y_new = wspace->node[i]->y;
		/*if(wspace->name =="W6"){
			clog<<"x0, y0: "<<x0<<" "<<y0<<endl;
			clog<<"x_old, y_old, x_new, y_new: "<<x_old<<" "<<y_old<<" "<<x_new<<" "<<y_new<<endl;
		}*/
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
		else if(x_old <x0 && y_old <y0){
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
	 //if(wspace->name == "W6")
	 //clog<<"wspace, degree: "<<wspace->name<<" "<<degree<<endl;
	if(degree == 360)
		flag = true;
	else if(degree == 0)
		flag = false;
	return flag;
}

bool ldo_in_wspace(LDO *ldo, MODULE *wspace){
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

/*bool ldo_in_wspace_trial(double ref_dist, double ref_x, double ref_y, double x0, double y0, LDO &ldo, MODULE *wspace){
	bool in_space_flag = false;
	double x1, y1;
	double x, y;
	double dx, dy;
	double dist;
		
	int width = ldo.width;
	int height = ldo.height;

	// try all 8 positions first to find a place
	vector<int> vec_y; // diagonal corner coordinate
	vector<int> vec_x;
	
	if(width == height){
		vec_y.resize(4);
		vec_x.resize(4);
	}else{
		vec_y.resize(8);
		vec_x.resize(8);
	}
	// case 1
	vec_x[0] = x0 + width-1; vec_y[0] = y0 + height-1;
	vec_x[1] = x0 - width+1; vec_y[1] = y0 + height-1;
	vec_x[2] = x0 - width+1; vec_y[2] = y0 - height+1;
	vec_x[3] = x0 + width-1; vec_y[3] = y0 - height+1;

	if(width != height){
		// the other 4
		vec_x[4] = x0 + height-1; vec_y[4] = y0 + width-1;
		vec_x[5] = x0 - height+1; vec_y[5] = y0 + width-1;
		vec_x[6] = x0 - height+1; vec_y[6] = y0 - width+1;
		vec_x[7] = x0 + height-1; vec_y[7] = y0 - width+1;
	}

	bool return_flag = false;
	double x2_new, y2_new;
		
	// scan all the 8 different shapes
	for(int i=0;i<vec_x.size();i++){
		x2_new = vec_x[i];
		y2_new = vec_y[i];
		return_flag = node_in_wspace(vec_x[i], vec_y[i], wspace);
		// if find ont, process overlap
		if(return_flag){
			break;
		}
	}	
	vec_x.clear();
	vec_y.clear();
	return return_flag;
}*/

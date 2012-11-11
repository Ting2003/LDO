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

// simple version utilized in moving process
// at least one corner should be OK for the ldo
bool ldo_in_wspace_trial(double ref_dist, double ref_x, double ref_y, double &x2, double &y2, LDO &ldo, WSPACE *wspace){
	bool in_space_flag = false;
	double x1, y1;
	double x, y;
	double dx, dy;
	double dist;
		
	int width = ldo.width;
	int height = ldo.height;
	int x0 = ldo.node[0]->x;
	int y0 = ldo.node[0]->y;
	//clog<<"x0, y0: "<<x0<<" "<<y0<<endl;

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
	bool flag_1, flag_2;
	double x1_new, y1_new, x2_new, y2_new;
	queue<int> q;
	int id_cur;
	int id_nbr;
	int step[2] = {1,-1};
	vector<bool> process;
	process.resize(wspace->node.size(), false);

	// first see if this node is already within wspace
	bool flag = node_in_wspace(ref_x, ref_y, wspace);
	clog<<"target in wspace: "<<flag<<endl;
	// directly set LDO here
	if(flag ==  true){
			
		//return;
	}	
	// need to move, put the white space as a queue
	// radius search for position around optimum place
	// first push back optimial node
	int min_id;
	for(size_t i=0;i<wspace->node.size();i++){
		x = wspace->node[i]->x;
		y = wspace->node[i]->y;
		if(x == x0 && y == y0){
			min_id = i;
			break;
		}
	}

	// scan all the corners
	q.push(min_id);
	process[min_id] = true;
	while(!q.empty()&&!return_flag){
		id_cur = q.front();
		// clog<<"id_cur: "<<id_cur<<endl;
		id_nbr = id_cur;
		x = wspace->node[id_cur]->x;
		y = wspace->node[id_cur]->y;
		//clog<<"x, y: "<<x<<" "<<y<<endl;
		
		dx = x - ref_x;
		dy = y - ref_y;
		dist = sqrt(dx*dx + dy*dy);
		// avaliable
		if(dist < ref_dist){
			for(int i=0;i<vec_x.size();i++){
				x2_new = vec_x[i];
				y2_new = vec_y[i];
			//clog<<"x, y, xnew, ynew: "<<x<<" "<<y<<" "<<x2_new<<" "<<y2_new<<endl;
				flag_2 = node_in_wspace(vec_x[i], vec_y[i], wspace);
			//clog<<"flag: "<<flag_2<<endl;
			// not satisfied, continue
				if(flag_2){
					// satisfied
					ldo.node[0]->x = x;
					ldo.node[0]->y = y;
					//x1 = x; y1 = y; 
					x2 = vec_x[i]; 
					y2 = vec_y[i];
				
				 clog<<"x, y, xnew, ynew, name: "<<x<<" "<<y<<" "<<x2<<" "<<y2<<" "<<wspace->name<<endl;
					return_flag = true;
					break;
				}
			}
		}	
		for(size_t i=0;i<2;i++){
			id_nbr = id_cur + step[i];
			if(id_nbr < wspace->node.size()
			  && id_nbr >= 0 && 
			  process[id_nbr] == false){
				q.push(id_nbr);
				//clog<<"push: "<<id_nbr<<endl;
				process[id_nbr] = true;
			}
		}
		q.pop();
	}
	// clear queue
	while(!q.empty())
		q.pop();
	vec_x.clear();
	vec_y.clear();
	return return_flag;
}

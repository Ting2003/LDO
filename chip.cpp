#include "chip.h"
#include "common.h"

Chip::Chip(){
	cktlist.clear();
	ldolist.clear();
}

void Chip::extract_rec(){
	for(size_t i=0;i<wspacelist.size();i++){
		extract_rec_single(i);
	}
}

// extract rectangles from whitespace
void Chip::extract_rec_single(size_t k){
	RECTANGLE * rec;
	WSPACE *wspace = wspacelist[k];
	vector<int> x;
	vector<int> y;
	int xl, xr, yb, yt;
	double center_x, center_y;
	bool in_wspace = false;
	stringstream sstream;
	int count_rec = 0;
	// collect all the x and ys
	for(size_t j=0;j<wspace->node.size();j++){
		x.push_back(wspace->node[j]->x);
		y.push_back(wspace->node[j]->y);	
	}
	unique(x.begin(), x.end());
	unique(y.begin(), y.end());
	// form rectangle
	for(size_t i=0;i<x.size()-1;i++){
		xl = x[i]; xr = x[i+1];
		center_x = (xl+xr)*0.5;
		for(size_t j=0;j<y.size()-1;j++){
			yb = y[j]; yt = y[j+1];
			center_y = (yb+yt)*0.5;
			// judge wheter this node is outside of this whitespace
			in_wspace = is_in_wspace(center_x, center_y, wspace);
			if(in_wspace == false)
				continue;
			rec = new RECTANGLE();
			rec->xl = xl; rec->xr = xr;	
			rec->yb = yb; rec->yt = yt;
			count_rec++;
			sstream.str("");
			sstream<<"rec"<<count_rec;
			rec->name = sstream.str();
			wspace->rectangle.push_back(rec);
		}
	}
	x.clear();
	y.clear();
}


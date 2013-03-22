#include "chip.h"

Chip::Chip(){
	cktlist.clear();
	ldolist.clear();
}

/*void Chip::extract_rec(){
	// clog<<"wspacelist.size(): "<<wspacelist.size()<<endl;
	for(size_t i=0;i<wspacelist.size();i++){
		extract_rec_single(i);
	}
}

// extract rectangles from whitespace
void Chip::extract_rec_single(size_t k){
	// clog<<"entering extract rec single. "<<endl;
	RECTANGLE * rec;
	MODULE *wspace = wspacelist[k];
	vector<int> x;
	vector<int> y;
	int xl, xr, yb, yt;
	double center_x, center_y;
	bool in_wspace = false;
	stringstream sstream;
	int count_rec = 0;
	vector<int>::iterator it;
	// collect all the x and ys
	for(size_t j=0;j<wspace->node.size();j++){
		// clog<<"x, y: "<<wspace->node[j]->x<<" "<<wspace->node[j]->y<<endl;
		x.push_back(wspace->node[j]->x);
		y.push_back(wspace->node[j]->y);	
	}
	sort(x.begin(), x.end());
	it = unique(x.begin(), x.end());	
	x.resize(it - x.begin());
	
	sort(y.begin(), y.end());
	it = unique(y.begin(), y.end());
	y.resize(it - y.begin());
	
	// form rectangle
	for(size_t i=0;i<x.size()-1;i++){
		xl = x[i]; xr = x[i+1];
		center_x = (xl+xr)*0.5;
		for(size_t j=0;j<y.size()-1;j++){
			yb = y[j]; yt = y[j+1];
			center_y = (yb+yt)*0.5;
			// clog<<"cx, cy: "<<center_x<<" "<<center_y<<endl;
			// judge wheter this node is outside of this whitespace
			in_wspace = is_in_wspace(center_x, center_y, wspace);
			// clog<<"is in space: "<<center_x<<" "<<center_y<<" "<<in_wspace<<endl;
			if(in_wspace == false)
				continue;
			rec = new RECTANGLE();
			rec->xl = xl; rec->xr = xr;	
			rec->yb = yb; rec->yt = yt;
			//clog<<"rec.xl, xr, yb, yt: "<<
				//xl<<" "<<xr<<" "<<yb<<" "<<yt<<endl;
			count_rec++;
			sstream.str("");
			sstream<<"ws"<<k<<"::rec"<<count_rec;
			rec->name = sstream.str();
			
			//clog<<"rec.name: "<<rec->name<<endl;
			wspace->rectangle.push_back(rec);
		}
	}
	x.clear();
	y.clear();
}*/


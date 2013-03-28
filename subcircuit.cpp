// ----------------------------------------------------------------//
// Filename : SubCircuit.cpp
// Author : Zigang Xiao <zxiao2@illinois.edu>
//          Ting Yu <tingyu1@illinois.edu>
//
// implementation file of SubCircuit.h
// ----------------------------------------------------------------//
// - Zigang Xiao - Sun Jan 30 18:35:56 CST 2011
//   * Add UMFPACK support
// - Zigang Xiao - Tue Jan 25 17:19:21 CST 2011
//   * added framework of PCG
// - Zigang Xiao - Tue Jan 18 21:46:13 CST 2011
//   * added solve() and related function
// - Zigang Xiao - Sun Jan 16 16:04:03 CST 2011
//   * added this log

#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <ctime>
#include <utility>
#include <cassert>
#include <vector>
#include "subcircuit.h"
#include "util.h"
#include "algebra.h"
#include "node.h"
using namespace std;

//////////////////////////////////////////////////////////////////////////
// Constructor and utility functions goes here

vector<LAYER_DIR> SubCircuit::layer_dir(MAX_LAYER);

// constructor of SubCircuit class, name is optional
SubCircuit::SubCircuit(string _name):name(_name),
	circuit_type(UNKNOWN), VDD(0.0){
	// add ground node
	Node * gnd = new Node(string("0"), Point(-1,-1,-1));
	gnd->rep = gnd;
	this->add_node(gnd);

	for(int i=0;i<MAX_LAYER;i++)
		layer_dir[i]=NA;
	id_map = NULL;
	Lx = NULL;
	Li = NULL;
	Lp = NULL;
	Lnz = NULL;
}

// Trick: do not release memory to increase runtime
SubCircuit::~SubCircuit(){
	//delete [] id_map;
	Lx = NULL;
	Li = NULL;
	Lp = NULL;
	Lnz = NULL;
}

void SubCircuit::check_sys() const{
	clog<<"**** CHECKING SYSTEM ENVIRONMENT ****"<<endl;
	clog<<"* int size     = "<< sizeof(int)<<endl;
	clog<<"* long size    = "<< sizeof(long)<<endl;
	clog<<"* size_t size  = "<< sizeof(size_t)<<endl;
	clog<<"* UF_long size = "<< sizeof(UF_long)<<endl;
	clog<<"* Max nodelist = "<<(size_t)nodelist.max_size()<<endl;
	clog<<"****            END              ****"<<endl<<endl;
}

// functor to be used in STL sort function
// order: y > x > z > flag 
// input: two node a, b
// return true if a < b, false o/w
// note that ground note are put to last
bool compare_node_ptr(const Node * a, const Node * b){
	if( a->is_ground() ) return false;
	if (b->is_ground() ) return true;

	if( a->pt.y == b->pt.y ){
		if( a->pt.x == b->pt.x ){
			if( a->pt.z == b->pt.z ){
				return (a->isS() > b->isS());
			}
			else{
				return (a->pt.z > b->pt.z);// top down
			}
		}
		else
			return ( a->pt.x < b->pt.x );
	}
	else
		return (a->pt.y < b->pt.y);
}

// sort the nodes according to their coordinate 
void SubCircuit::sort_nodes(){
   size_t i=0;
   sort(nodelist.begin(), nodelist.end(), compare_node_ptr);
   // update node id mapping, 
   // NOTE: ground node will be the last
   // if nodelist_size < 1e6, use serial
   for(i=0;i<nodelist.size();i++){
	   Node * p = nodelist[i];
	   node_id[p] = i;
   }
}

string SubCircuit::get_name() const{return this->name;}

ostream & operator << (ostream & os, const NodePtrVector & nodelist){
	for(size_t i=0;i<nodelist.size();i++)
		os<<*nodelist[i]<<endl;
	return os;
}

ostream & operator << (ostream & os, const NetPtrVector & nets){
	for(size_t i=0;i<nets.size();i++)
		os<<*nets[i]<<endl;
	return os;
}

ostream & operator << (ostream & os, const SubCircuit & ckt){
	os<<"SubCircuit ["<<ckt.name<<"] info:"<<endl;

	os<<"==== Nodes ===="<<endl;
	os<<ckt.nodelist;

	os<<"==== Reps  ===="<<endl;
	os<<ckt.replist;

	return os;
}

void SubCircuit::print(){
	// uncomment this if want to output to a file
	//freopen("output.txt","w",stdout);

	// don't output ground node
	for(size_t i=0;i<nodelist.size()-1;i++){
		printf("%s  %.5e\n", nodelist[i]->name.c_str(), 
				nodelist[i]->value);
	}
}

///////////////////////////////////////////////////////////////////////////////
// Computation Functions

// initialization before solving the SubCircuit
// 1. sort the nodes
// 2. set node representatives
// 3. find node in which block, update count
// 4. get representative lists
void SubCircuit::solve_init(bool flag){
	replist.clear();
        sort_nodes();
	size_t size = nodelist.size()-1; // exclude ground node!!
	Node *p=NULL;
	for(size_t i=0,nr=0;i<size;i++){
		p=nodelist[i];

		// set the representative
		Net * net = p->nbr[TOP];
		if(flag == true){ // global grid
			if( p->isS()== Y) 
				VDD = p->get_value();
		}
		else	VDD = VDD_G;

		// test short SubCircuit
		if( p->isS() !=Y && // Y must be representative 
		    net != NULL &&
		    fzero(net->value) ){
			// TODO: ensure ab[1] is not p itself
			assert( net->ab[1] != p );
			p->rep = net->ab[1]->rep;
		} // else the representative is itself

		// push the representatives into list
		if( p->rep == p ) {
			replist.push_back(p);
			//rep_id[p] = nr; // save the id
			p->rid = nr;
			++nr;
		}
		else
			p->rid = p->rep->rid;
		//clog<<"p->rep: "<<*p->rep<<endl;
	}
	
	//clog<<"before build pad set. "<<endl;
	// build_pad_set();

	//clog<<"nodelist.size: "<<nodelist.size()<<endl;
	//clog<<"replist.size: "<<replist.size()<<endl;
}

// count number of nodes that can be merged
void SubCircuit::count_merge_nodes(){
	size_t size = replist.size();
	size_t count = 0;
	for(size_t i=0;i<size;i++){
		Node * p = replist[i];
		if(p->nbr[TOP] ==NULL && p->nbr[BOTTOM] ==NULL){
			if((p->nbr[EAST] !=NULL && p->nbr[WEST] !=NULL) 
			&& (p->nbr[NORTH] ==NULL && p->nbr[SOUTH]==NULL))
				count++;
			else if((p->nbr[NORTH] !=NULL && p->nbr[SOUTH] !=NULL)
					&&(p->nbr[EAST]==NULL && p->nbr[WEST] ==NULL))
				count++;
		}
	}
	//clog<<"number of nodes can be merged is: "<<count<<endl;
}

// change the wspace into a local variable
// here is simply copy
void SubCircuit::build_wspacelist(vector<MODULE*> wspace_vec){
	wspacelist.clear();
	wspacelist = wspace_vec;
}

void SubCircuit::build_ldolist(vector<LDO*> ldo_vec){
	Node *nd;
	ldolist.clear();
	for(size_t i=0;i<ldo_vec.size();i++){
		nd = ldo_vec[i]->A;
		if(has_node(nd->name))
			ldolist.push_back(ldo_vec[i]);
	}
	// clog<<"ldolist.size: "<<ldolist.size()<<endl;
	//clog<<"gx, gy: "<<gx<<" "<<gy<<endl; 
}
#if 0
// solve the node voltages using direct LU
void SubCircuit::solve_LU(Tran &tran, bool local_flag){
        solve_init(local_flag);
	clog<<"subckt enter solve LU core. "<<endl;
	if(local_flag == true)
		// build local VDD nets
		build_local_nets();
	else    // build global current nets
		build_global_nets();
	/*for(int type =0 ;type <NUM_NET_TYPE;type++){
		NetList &ns = net_set[type];
		NetList::iterator it;
		for(it = ns.begin();it!=ns.end();++it)
			cout<<*(*it)<<endl;
	}*/
	solve_LU_core(tran, local_flag);
}
#endif

// given vector x that obtained from LU, set the value to the corresponding
// node in nodelist
double SubCircuit::get_voltages_from_LU_sol(double * x){
   size_t i;
   double max_diff = 0;
   double diff = 0;
   for(i=0;i<nodelist.size()-1;i++){
      Node * node = nodelist[i];
      size_t id = node->rep->rid;  // get rep's id in Vec
      double v = x[id];		// get its rep's value
      //cout<<"node_vol, x, diff: "<<node->value<<" "<<v<<" "<<diff<<endl;
      diff = abs(node->value - v);
      if(diff  > max_diff)
	      max_diff = diff;
      node->value = v;
   }
   return max_diff;
}

// copy node voltages from the SubCircuit to a Vec
// from = true then copy SubCircuit to x
// else copy from x to SubCircuit
void SubCircuit::copy_node_voltages(double * x, size_t &size, bool from){
	if( from == true ){
		for(size_t i=0;i<size;i++)
			x[i] = nodelist[i]->value;
	}
	else{
		for(size_t i=0;i<size;i++){
			Node * node = nodelist[i];
			size_t id = node_id[node->rep];
			double v = x[id];
			node->value = v;
		}
	}
}

// stamp the net in each set, 
// *NOTE* at the same time insert the net into boundary netlist
void SubCircuit::stamp_by_set(Matrix & A){
   	A.clear();
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		switch(type){
		case RESISTOR:
			for(size_t i=0;i<ns.size();i++){
				assert( fzero(ns[i]->value) == false );
				stamp_resistor(A, ns[i]);
			}
			break;
		case CURRENT:
			break;
		case VOLTAGE:
			for(size_t i=0;i<ns.size();i++){
				if( fzero(ns[i]->value)  && 
				    !ns[i]->ab[0]->is_ground() &&
				    !ns[i]->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_VDD(A, ns[i]);
			}
			break;
		case CAPACITANCE:
			//for(size_t i=0;i<ns.size();i++)
				//stamp_capacitance_dc(A, ns[i]);
			break;
		case INDUCTANCE:
			for(size_t i=0;i<ns.size();i++){
				stamp_inductance_dc(A, ns[i]);	
			}
			break;
		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}
}

// stamp transient current values into rhs
void SubCircuit::stamp_current_tr(double *b, double &time){
	NetPtrVector & ns = net_set[CURRENT];
	for(size_t i=0;i<ns.size();i++)
		stamp_current_tr_net(b, ns[i], time);
}

// stamp transient current values into rhs
void SubCircuit::stamp_current_tr_1(double *bp, double *b, double &time){
	NetPtrVector & ns = net_set[CURRENT];
	//if(ns.size()<THRESHOLD){
		for(size_t i=0;i<ns.size();i++)
			stamp_current_tr_net_1(bp, b, ns[i], time);
	//}
#if 0
	// else work in parallel
	else{
		size_t i=0;
#pragma omp parallel for private(i)
		for(i=0;i<ns.size();i++)
			stamp_current_tr_net_1(bp, b, ns[i], time);
	}
#endif
}

// stamp the transient matrix
void SubCircuit::stamp_by_set_tr(Matrix & A, Tran &tran){
   	A.clear();
	
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		switch(type){
		case RESISTOR:
			for(size_t i=0;i<ns.size();i++){
				assert( fzero(ns[i]->value) == false );
				stamp_resistor_tr(A, ns[i]);
			}
			break;
		case CURRENT:
			//for(size_t i=0;i<ns.size();i++)
				//stamp_current(b, ns[i]);
			break;
		case VOLTAGE:
			for(size_t i=0;i<ns.size();i++){
				if( fzero(ns[i]->value)  && 
				    !ns[i]->ab[0]->is_ground() &&
				    !ns[i]->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_VDD(A, ns[i]);
			}
			break;
		case CAPACITANCE:
			for(size_t i=0;i<ns.size();i++)
				stamp_capacitance_tr(A, ns[i], tran);
			break;
		case INDUCTANCE:
			for(size_t i=0;i<ns.size();i++){
				stamp_inductance_tr(A, ns[i], tran);	
			}
			break;
		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}
}

// update rhs by transient nets
void SubCircuit::modify_rhs_tr_0(Tran &tran){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		if(type ==CAPACITANCE){	
			for(size_t i=0;i<ns.size();i++)
				modify_rhs_c_tr_0(ns[i], tran);
		}
		else if(type == INDUCTANCE){
			for(size_t i=0;i<ns.size();i++){
				modify_rhs_l_tr_0(ns[i], tran);	
			}
		}
	}
}

// update rhs by transient nets
void SubCircuit::modify_rhs_tr(double * b, double *x){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		if(type ==CAPACITANCE){
			for(size_t i=0;i<ns.size();i++)
				modify_rhs_c_tr(ns[i], b, x);
		}
		else if(type == INDUCTANCE){
			for(size_t i=0;i<ns.size();i++){
				modify_rhs_l_tr(ns[i], b, x);
			}
		}
	}
}

void SubCircuit::stamp_resistor(Matrix & A, Net * net){
	//cout<<"net: "<<*net<<endl;
	double G;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	G = 1./net->value;
	//cout<<"net: "<<*net<<endl;
        if( !nk->is_ground()&& nk->isS()!=Y && 
          (nk->nbr[TOP]== NULL 
           || nk->nbr[TOP]->type != INDUCTANCE)) {
           A.push_back(k,k, G);

           // cout<<"("<<k<<" "<<k<<" "<<G<<")"<<endl;
           if(!nl->is_ground() && nl->isS()!=Y
			   &&(nl->nbr[TOP]==NULL || 
                 nl->nbr[TOP]->type != INDUCTANCE)&&(k > l)){
                    A.push_back(k,l,-G);

                  // cout<<"("<<k<<" "<<l<<" "<<-G<<")"<<endl;
           }
        }

	if( !nl->is_ground() && nl->isS() !=Y && 
			(nl->nbr[TOP] ==NULL 
		||nl->nbr[TOP]->type != INDUCTANCE)) {
		A.push_back(l,l, G);
                // cout<<"("<<l<<" "<<l<<" "<<G<<")"<<endl;

		if(!nk->is_ground()&& nk->isS()!=Y
			&&(nk->nbr[TOP]==NULL ||
		  nk->nbr[TOP]->type != INDUCTANCE) && l > k){
			A.push_back(l,k,-G);
		 //cout<<"("<<l<<" "<<k<<" "<<-G<<")"<<endl;
                }
	}
}

// stamp all resistor net
void SubCircuit::stamp_resistor_tr(Matrix & A, Net * net){
   double G;
   Node * nk = net->ab[0]->rep;
   Node * nl = net->ab[1]->rep;
   size_t k = nk->rid;
   size_t l = nl->rid;
   G = 1./net->value;

   // cout<<"net: "<<*net<<endl;
   if( nk->isS()!=Y && !nk->is_ground()){
           A.push_back(k,k, G);
      // cout<<"("<<k<<" "<<k<<" "<<G<<")"<<endl;
      if(!nl->is_ground() && nl->isS()!=Y && l<k){
            A.push_back(k,l,-G);
            // cout<<"("<<k<<" "<<l<<" "<<-G<<")"<<endl;
      }
   }

   if( nl->isS() !=Y && !nl->is_ground()){
      A.push_back(l,l, G);
      // cout<<"("<<l<<" "<<l<<" "<<G<<")"<<endl;
      if(!nk->is_ground() && nk->isS()!=Y && k<l){
            A.push_back(l,k,-G);
            // cout<<"("<<l<<" "<<k<<" "<<-G<<")"<<endl;
      }
   }
}

void SubCircuit::stamp_induc_rhs_dc(double *b, Net * net){
	//clog<<"net: "<<*net<<endl;
	double G;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	G = 1./net->value;
	if( nk->isS()!=Y && !nk->is_ground()) {
		// general stamping
		if(!nl->is_ground())
			b[k] = b[l];
	}

	if( nl->isS() !=Y && !nl->is_ground()) {
		if(!nk->is_ground())
		// general stamping
		b[l] = b[k];
	}
}

void SubCircuit::stamp_inductance_dc(Matrix & A, Net * net){
	//clog<<"net: "<<*net<<endl;
	double G;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	G = 1./net->value;
	if( nk->isS()!=Y && !nk->is_ground()) {
		A.push_back(k,k, 1);
		// general stamping
	}

	if( nl->isS() !=Y && !nl->is_ground()) {
		A.push_back(l,l, 1);
	}
}

// stamp dc cap, infinitive resistance
void SubCircuit::stamp_capacitance_dc(Matrix & A, Net * net){
	double G;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	G = 1./net->value;
	if( nk->isS()!=Y && !nk->is_ground()) {
		A.push_back(k,k, 0);
		if(!nl->is_ground())
			A.push_back(k,l,0);
	}

	if( nl->isS() !=Y && !nl->is_ground()) {
		A.push_back(l,l, 0);
		if(!nk->is_ground())
			A.push_back(l,k,0);
	}
}

// stamp inductance Geq = delta_t/(2L)
void SubCircuit::stamp_inductance_tr(Matrix & A, Net * net, Tran &tran){
	// cout<<"induc net: "<<*net<<endl;
	double Geq = 0;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	// Geq = delta_t / (2*L)
	Geq = tran.step_t / (2*net->value);
	//net->value = Geq;

	if( nk->isS()!=Y  && !nk->is_ground()) {
		// -1 is to clear formal inserted 1 at (k,k)
		A.push_back(k,k, Geq);
		// cout<<"("<<k<<" "<<k<<" "<<Geq<<")"<<endl;
		//clog<<nl->isS()<<endl;
		if(!nl->is_ground()&& nl->isS()!=Y && k>l){
			A.push_back(k,l,-Geq);
		        // cout<<"("<<k<<" "<<l<<" "<<-Geq<<")"<<endl;
		}
	}

	if( nl->isS() !=Y && !nl->is_ground()) {
		// -1 is to clear formal inserted 1 at (l,l)
		A.push_back(l,l, Geq);
		 cout<<"("<<l<<" "<<l<<" "<<Geq<<")"<<endl;
		if(!nk->is_ground() && nk->isS()!=Y && l>k){
			A.push_back(l,k,-Geq);
			 cout<<"("<<l<<" "<<k<<" "<<-Geq<<")"<<endl;
		}
	}
}

// stamp capacitance Geq = 2C/delta_t
void SubCircuit::stamp_capacitance_tr(Matrix &A, Net *net, Tran &tran){
	// cout<<"cap net: "<<*net<<endl;
	double Geq = 0;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	// Geq = 2*C / delta_t
	Geq = (2*net->value) / tran.step_t;
	//net->value = Geq;
	// clog<<"C delta_t Geq: "<<net->value<<" "<<tran.step_t<<" "<<Geq<<endl;
	// Ieq = i(t) + 2*C / delta_t * v(t)

	if( nk->isS()!=Y  && !nk->is_ground()) {
		A.push_back(k,k, Geq);
		// cout<<"("<<k<<" "<<k<<" "<<Geq<<")"<<endl;
		if(!nl->is_ground()&& k > l){
			A.push_back(k,l,-Geq);
			// cout<<"("<<k<<" "<<l<<" "<<-Geq<<")"<<endl;
		}
	}

	if( nl->isS() !=Y && !nl->is_ground()) {
		A.push_back(l,l, Geq);
		// cout<<"("<<l<<" "<<l<<" "<<Geq<<")"<<endl;
		if(!nk->is_ground()&& l > k){
			A.push_back(l,k,-Geq);
			// cout<<"("<<l<<" "<<k<<" "<<-Geq<<")"<<endl;
		}
	}
}

// add Ieq into rhs
// Ieq = i(t) + 2*C / delta_t *v(t)
void SubCircuit::modify_rhs_c_tr_0(Net *net, Tran &tran){
	double i_t = 0;
	double temp = 0;
	double Ieq = 0;
	// clog<<"c net: "<<*net<<endl;
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
        // nk point to Z node
        if(nk->isS() != Z){
		swap<Node *>(nk, nl);
		swap<Node*>(net->ab[0], net->ab[1]);
	}
	//clog<<"nk, nl: "<<*nk<<" "<<*nl<<endl;
	//size_t k = nk->rid;
	//size_t l = nl->rid;

	Net *r = nk->nbr[TOP];
	Node *a = r->ab[0]->rep;
	Node *b = r->ab[1]->rep;
	// a point to Z node
	if(a->isS()!=Z) {
		swap<Node *>(a, b);
		swap<Node*>(r->ab[0], r->ab[1]);
	}
	//clog<<"a, b: "<<*a<<" "<<*b<<endl;

	//size_t id_a = a->rid;
	//size_t id_b = b->rid;
	i_t = (b->value - a->value) / r->value;
	// i_t = (x[id_b] - x[id_a]) / r->value;
	       
        // push 2 nodes into node_set_x
        //clog<<*nk<<" "<<k<<endl;
	if(nk->is_ground())
	// temp = 2*net->value/tran.step_t*(0-x[l]);
		temp = 2*net->value/tran.step_t *(-nl->value);
        else if(nl->is_ground()){
         //temp = 2*net->value/tran.step_t *(x[k]);
		temp = 2*net->value / tran.step_t *nk->value;
        }
        else
        //temp = 2*net->value/tran.step_t *(x[k] - x[l]);
        	temp = 2*net->value/tran.step_t *
			(nk->value - nl->value);
		
	Ieq  = (i_t + temp);
	net->Ieq = Ieq;
	//clog<< "Ieq is: "<<Ieq<<endl;
	//clog<<"Geq is: "<<2*net->value / tran.step_t<<endl;
	/*if(!nk->is_ground()&& nk->isS()!=Y){
		 rhs[k] += Ieq;	// for VDD SubCircuit
		// clog<<*nk<<" rhs +: "<<rhs[k]<<endl;
	}
	if(!nl->is_ground()&& nl->isS()!=Y){
		 rhs[l] += -Ieq; 
		//clog<<*nl<<" rhs +: "<<rhs[l]<<endl;
	}*/
}

// add Ieq into rhs
// Ieq = i(t) + 2*C / delta_t *v(t)
void SubCircuit::modify_rhs_c_tr(Net *net, double * rhs, double *x){
	double temp = 0;
	//clog<<"c net: "<<*net<<endl;
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
        
	// nk point to Z node
	size_t k = nk->rid;
	size_t l = nl->rid;

	Net *r = nk->nbr[TOP];
	Node *a = r->ab[0]->rep;
	Node *b = r->ab[1]->rep;
	// a point to Z node

	size_t id_a = a->rid;
	size_t id_b = b->rid;
	double i_t = (x[id_b] - x[id_a]) / r->value;
	
	if(nk->is_ground())
	 temp = net->value *(-x[l]);
        else if(nl->is_ground())
	 temp = net->value *x[k];
        else
	 temp = net->value *(x[k]-x[l]);
	
	double Ieq  = i_t + temp;
	if(!nk->is_ground()&& nk->isS()!=Y){
		 rhs[k] += Ieq;	// for VDD SubCircuit
	}
	if(!nl->is_ground()&& nl->isS()!=Y){
		 rhs[l] -= Ieq; 
	}
}

// set eq resistance of induc
void SubCircuit::set_eq_induc(Tran &tran){
	NetPtrVector &ns = net_set[INDUCTANCE];
	for(size_t i=0;i<ns.size();i++)
		ns[i]->value = tran.step_t /(2*ns[i]->value);
}

// set eq resistance of capac
void SubCircuit::set_eq_capac(Tran &tran){
	NetPtrVector &ns = net_set[CAPACITANCE];
	for(size_t i=0;i<ns.size();i++)
		ns[i]->value = 2*ns[i]->value/tran.step_t;
}

// add Ieq into rhs
// Ieq = i(t) + delta_t / (2*L) *v(t)
void SubCircuit::modify_rhs_l_tr_0(Net *net, Tran &tran){
	//  clog<<"l net: "<<*net<<endl;
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
	// nk point to X node
	if(nk->isS() !=X){ 
		swap<Node*>(nk, nl);
		swap<Node*>(net->ab[0], net->ab[1]);
	}
	//size_t k = nk->rid;
	//size_t l = nl->rid;
	double Ieq = 0;

	double i_t = 0;
	double temp = 0;
	temp = tran.step_t / (2*net->value) * 
		(nl->value - nk->value);
	// temp = tran.step_t / (2*net->value)*(x[l] - x[k]);
	//clog<<"Geq: "<<tran.step_t / (2*net->value)<<endl;
	// temp = net->value *(x[l] - x[k]);	
	
	//clog<<"delta_t/2L, nl-nk, temp: "<<tran.step_t / (2*net->value)<<" "<<(nl->value-nk->value)<<" "<<temp<<endl;
	
	Net *r = nk->nbr[BOTTOM];
	Node *a = r->ab[0]->rep;
	Node *b = r->ab[1]->rep;
	// a point to X node
	if(a->isS()!=X) {
		swap<Node*>(a, b);
		swap<Node*>(r->ab[0], r->ab[1]);
	}
	//size_t id_a = a->rid;
	//size_t id_b = b->rid;
	//i_t = (x[id_a] - x[id_b]) / r->value;
	i_t = (a->value - b->value) / r->value;
	// clog<<"a->value, b->value, r->value, it: "<<a->value<<" "<<b->value<<" "<<r->value<<" "<<i_t<<endl;
	// clog<<"temp: "<<temp<<endl;
                
	// push inductance nodes into node_set_x
        //clog<<*nk<<" "<<k<<endl;
        //clog<<*b<<" "<<id_b<<endl;
	Ieq  = i_t + temp;
	net->Ieq = Ieq;
	// cout<<"Ieq: "<<Ieq<<endl;
	/*if(nk->isS() !=Y && !nk->is_ground()){
		 rhs[k] += Ieq; // VDD SubCircuit
		//clog<<*nk<<" "<<rhs[k]<<endl;
	}
	if(nl->isS()!=Y && !nl->is_ground()){
		 rhs[l] += -Ieq; // VDD SubCircuit
		//clog<<*nl<<" "<<rhs[l]<<endl;
	}*/
}

// add Ieq into rhs
// Ieq = i(t) + delta_t / (2*L) *v(t)
void SubCircuit::modify_rhs_l_tr(Net *net, double *rhs, double *x){
	//clog<<"l net: "<<*net<<endl;
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
	// nk point to X node
	//if(nk->isS() !=X){ 
		//swap<Node*>(nk, nl);
	//}
	size_t k = nk->rid;
	size_t l = nl->rid;
	double Ieq = 0;

	double i_t = 0;
	double temp = 0;
	//temp = tran.step_t / (2*net->value) * 
		//(nl->value - nk->value);
	//temp = tran.step_t / (2*net->value)*(x[l] - x[k]);
	temp = net->value *(x[l] - x[k]);	
	//if(nk->value != x[k] || nl->value != x[l])
	   //clog<<"k, l, x_k, x_l: "<<nk->value<<" "<<nl->value<<" "<<
	     //x[k]<<" "<<x[l]<<endl;

	//clog<<"delta_t/2L, nl-nk, temp: "<<tran.step_t / (2*net->value)<<" "<<(nl->value-nk->value)<<" "<<temp<<endl;
	
	Net *r = nk->nbr[BOTTOM];
	Node *a = r->ab[0]->rep;
	Node *b = r->ab[1]->rep;
	// a point to X node
	//if(a->isS()!=X) {
		//swap<Node*>(a, b);
	//}
	size_t id_a = a->rid;
	size_t id_b = b->rid;
	i_t = (x[id_a] - x[id_b]) / r->value;
	//i_t = (a->value - b->value) / r->value;
        //if(b->value != x[id_b] || a->value != x[id_a])
	   //clog<<"a, b, x_a, x_b: "<<a->value<<" "<<b->value<<" "<<
	     //x[id_a]<<" "<<x[id_b]<<endl;

	//clog<<"resiste r: "<<*r<<endl;
	//clog<<*a<<" "<<*b<<endl;
	//clog<<"a, b, r, i_t: "<<a->value<<" "<<b->value<<" "<<
		//r->value<<" "<<i_t<<endl;
       
        // push inductance nodes into node_set_x
        //clog<<*nk<<" "<<k<<endl;
        //clog<<*b<<" "<<id_b<<endl;
 #if 0
        if(iter==0){
           pg.node_set_x.push_back(k);
           pg.node_set_x.push_back(id_b);
        }
#endif
	Ieq  = i_t + temp;
	//clog<<"Ieq: "<<Ieq<<endl;
	if(nk->isS() !=Y && !nk->is_ground()){
		 rhs[k] += Ieq; // VDD SubCircuit
		//clog<<*nk<<" "<<rhs[k]<<endl;
	}
	if(nl->isS()!=Y && !nl->is_ground()){
		 rhs[l] += -Ieq; // VDD SubCircuit
		//clog<<*nl<<" "<<rhs[l]<<endl;
	}
}
// stamp a current source
void SubCircuit::stamp_current(double * b, Net * net){
	//cout<<"net: "<<*net<<endl;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;

	if( !nk->is_ground() && nk->isS()!=Y){// && 
		size_t k = nk->rid;
		b[k] += -net->value;
		//cout<<"b: "<<k<<" "<<-net->value<<endl;
	}
	if( !nl->is_ground() && nl->isS() !=Y){// &&
		size_t l = nl->rid;
		b[l] +=  net->value;
		//cout<<"b: "<<l<<" "<<-net->value<<endl;
	}
}

void SubCircuit::stamp_current_tr_net(double * b, Net * net, double &time){
	current_tr(net, time);
	// cout<<"net: "<<*net<<endl;
	//clog<<"current: "<<current<<endl;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	if( !nk->is_ground()&& nk->isS()!=Y) { 
		size_t k = nk->rid;
		//clog<<"node, rid: "<<*nk<<" "<<k<<endl;
		b[k] += -net->value;//current;
		// cout<<"time, k, b: "<<time<<" "<<k<<" "<<b[k]<<endl;
	}
	if( !nl->is_ground() && nl->isS()!=Y) {
		size_t l = nl->rid;
		//clog<<"node, rid: "<<*nl<<" "<<l<<endl;
		b[l] +=  net->value;// current;
		// cout<<"time, l, b: "<<time<<" "<<l<<" "<<b[l]<<endl;
	}
}

void SubCircuit::stamp_current_tr_net_1(double *bp, double * b, Net * net, double &time){
	double diff = 0;
	double current = net->value;
	current_tr(net, time);
	//if(time / 1e-11 == 22)
	//cout<<"time, co, cn: "<<time<<" "<<current<<" "<<net->value<<endl;
	// only stamps when net got a different current
	if(current != net->value){
		diff = net->value - current;
		
		//cout<<"time, old, new, diff: "<<time <<" "<<current<<" "<<net->value<<" "<<diff<<endl;
		//cout<<"net: "<<*net;
		//clog<<"current: "<<current<<endl;
		Node * nk = net->ab[0]->rep;
		Node * nl = net->ab[1]->rep;
		if( !nk->is_ground()&& nk->isS()!=Y) { 
			size_t k = nk->rid;
			//clog<<"node, rid: "<<*nk<<" "<<k<<endl;
			//clog<<"time, k, b bef: "<<time<<" "<<k<<" "<<b[k]<<endl;
			b[k] += -diff;//current;
			bp[k] = b[k];
			//clog<<"time, k, b: "<<time <<" "<<k<<" "<<b[k]<<endl;
		}
		if( !nl->is_ground() && nl->isS()!=Y) {
			size_t l = nl->rid;
			//clog<<"time, l, b bef: "<<time<<" "<<l<<" "<<b[l]<<endl;
			//clog<<"node, rid: "<<*nl<<" "<<l<<endl;
			b[l] +=  diff;// current;
			bp[l] = b[l];
			//clog<<"time, l, b: "<<time<<" "<<l<<" "<<b[l]<<endl;
		}
	}
}

// stamp a voltage source
// no current source at voltage source, save
void SubCircuit::stamp_VDD(Matrix & A, Net * net){
	// find the non-ground node
	// cout<<"vol net: "<<*net<<endl;
	Node * X = net->ab[0];
	if( X->is_ground() ) X = net->ab[1];
	size_t id = X->rep->rid;
	A.push_back(id, id, 1.0);
	//cout<<"push id, id, 1: "<<id<<" "<<id<<" "<<1<<endl;	
}

// decide transient step current values
void SubCircuit::current_tr(Net *net, double &time){
	double slope = 0;
	double Tr = net->Tr;
	double PW = Tr + net->PW;
	double Tf = PW + net->Tf;
	double t_temp = time - net->TD;
	double t = fmod(t_temp, net->Period);
	if(time <= net->TD)
		net->value = net->V1;
	else if(t > 0 && t<= Tr){
		slope = (net->V2 - net->V1) / 
			(net->Tr);
		net->value = net->V1 + t*slope;
	}
	else if(t > Tr && t<= PW)
		net->value = net->V2;
	else if(t>PW && t<=Tf){
		slope = (net->V1-net->V2)/(net->Tf);
		net->value = net->V2 + slope*(t-PW);
	}
	else
		net->value = net->V1;
	//return current;
}

// assign value back to transient nodes
void SubCircuit:: save_tr_nodes(Tran &tran, double *x){
   size_t id=0;
   for(size_t j=0;j<tran.nodes.size();j++){
      if(tran.nodes[j].node != NULL ){
         id = tran.nodes[j].node->rep->rid;
         tran.nodes[j].value.push_back(x[id]);
      }
   }
}

// assign value back to transient nodes
void SubCircuit:: save_ckt_nodes(){
   size_t id=0;
   // clog<<"ckt nodes number: "<<ckt_nodes.size()<<endl;
   for(size_t j=0;j<ckt_nodes.size();j++){
	 // clog<<"nodes: "<<ckt_nodes[j].node->name<<" / ";
         id = ckt_nodes[j].node->rep->rid;
	 // clog<<"value: "<<xp[id]<<endl;
         ckt_nodes[j].value.push_back(xp[id]);
      }
}

// link transient nodes with nodelist
void SubCircuit:: link_ckt_nodes(Tran &tran){
   Node_TR_PRINT nodes_temp;
   for(size_t i=0;i<nodelist.size();i++){
      for(size_t j=0;j<tran.nodes.size();j++){
         if(nodelist[i]->name == tran.nodes[j].name){
	    // record the index in tran.nodes
	    nodes_temp.flag = j;
	    //cout<<"tran.nodes, index: "<<
		//nodelist[i]->name<<" "<<nodes_temp.flag<<endl;
	    nodes_temp.node = nodelist[i];
	    ckt_nodes.push_back(nodes_temp);
            // record the id for ckt_nodes
            break;
         }
      }
   }
}

// link transient nodes with nodelist
void SubCircuit:: link_tr_nodes(Tran &tran){
   for(size_t i=0;i<nodelist.size();i++){
      for(size_t j=0;j<tran.nodes.size();j++){
         if(tran.nodes[j].flag ==1) continue;
         if(nodelist[i]->name == tran.nodes[j].name){
            tran.nodes[j].node = nodelist[i];
            tran.nodes[j].flag =1;
            // record the id for tran.nodes
            break;
         }
      }
   }
}

void SubCircuit:: release_tr_nodes(Tran &tran){
   for(size_t j=0;j<tran.nodes.size();j++){
      if(tran.nodes[j].flag != -1){
         tran.nodes[j].node = NULL;
      }
   }
}

void SubCircuit:: release_ckt_nodes(){
   for(size_t j=0;j<ckt_nodes.size();j++){
         ckt_nodes[j].node = NULL;
   }
}

// make local ldo resistor net symmetric
void SubCircuit::make_A_symmetric_local(double *b){
	int type = RESISTOR;
	NetList & ns = net_set[type];
	NetList::iterator it;
	Node *p, *q;

	for(it=ns.begin();it!=ns.end();it++){
		if( (*it) == NULL ) continue;
			assert( fzero((*it)->value) == false );
		if((*it)->ab[0]->isS()==Y && (*it)->ab[1]->isS()==Y)
			continue;
		// node a points to X node
		if((*it)->ab[0]->isS()==Y){
			p = (*it)->ab[0]; q = (*it)->ab[1];
		}
		else if((*it)->ab[1]->isS()==Y){
			p = (*it)->ab[1]; q = (*it)->ab[0];
		}
		else continue;
		size_t id = q->rep->rid;
		double G = 1.0 / (*it)->value;
		// clog<<"net: "<<*(*it)<<endl;
		// clog<<"id, p_val, G, b: "<<id<<" "<<p->value<<" "<<G<<" "<<b[id]<<endl;
		b[id] += p->value * G;
	}	
}

void SubCircuit::make_A_symmetric(double *b){
	int type = RESISTOR;
	NetList & ns = net_set[type];
	NetList::iterator it;
	Node *p=NULL, *q=NULL, *r =NULL;

	for(it=ns.begin();it!=ns.end();it++){
           if( (*it) == NULL ) continue;
           assert( fzero((*it)->value) == false );
           if(!((*it)->ab[0]->rep->isS()==X || (*it)->ab[1]->rep->isS()==X)) continue;
           // node p points to X node
           if((*it)->ab[0]->rep->isS()==X && ((*it)->ab[0]->rep->nbr[TOP]!=NULL && 
                (*it)->ab[0]->rep->nbr[TOP]->type==INDUCTANCE)){
              p = (*it)->ab[0]->rep; q = (*it)->ab[1]->rep;
           } 
           else if((*it)->ab[1]->rep->isS()==X && ((*it)->ab[1]->rep->nbr[TOP]!=NULL && 
                (*it)->ab[1]->rep->nbr[TOP]->type==INDUCTANCE)){
              p = (*it)->ab[1]->rep; q = (*it)->ab[0]->rep;
           }           
           r = p->nbr[TOP]->ab[0]->rep;
           if(r->isS()!=Y) 
              r = p->nbr[TOP]->ab[1]->rep;

           size_t id = q->rid;
           double G = 1.0 / (*it)->value;
           
           b[id] += r->value * G;
        }
}

// only inductance is connected to voltage source
void SubCircuit::make_A_symmetric_tr(double *b, Tran &tran){
	int type = INDUCTANCE;
	NetList & ns = net_set[type];
	NetList::iterator it;
	Node *p, *q;

	for(it=ns.begin();it!=ns.end();it++){
           if( (*it) == NULL ) continue;
           assert( fzero((*it)->value) == false );
           if(!((*it)->ab[0]->rep->isS()==Y || (*it)->ab[1]->rep->isS()==Y)) continue;
           // clog<<"net: "<<*(*it)<<endl;
           // node p points to Y node
           if((*it)->ab[0]->rep->isS()==Y){
              p = (*it)->ab[0]->rep; q = (*it)->ab[1]->rep;
           } 
           else if((*it)->ab[1]->rep->isS()==Y ){
              p = (*it)->ab[1]->rep; q = (*it)->ab[0]->rep;
           }           
           // clog<<"p and q: "<<*p<<" "<<*q<<endl;

           size_t id = q->rid;
	   // clog<<"id: "<<q->rid<<endl;
           double G = tran.step_t / ((*it)->value*2);
           
           b[id] += p->value * G;
           // b[id] += x[p->rid] *G;
           // clog<<"id, stamp p->value, G, b: "<<id<<" "<<p->value<<" "<<G<<" "<<b[id]<<endl;
        }
}

bool compare_Node_G(const Node_G *nd_1, const Node_G *nd_2){
   return (nd_1->value < nd_2->value);
 }

void SubCircuit::print_ckt_nodes(Tran &tran){
	double time = 0;
	size_t j=0;
	//clog<<"ckt, nodes: "<<get_name()<<" "<<ckt_nodes.size()<<endl;
	for(size_t i=0;i<ckt_nodes.size();i++){
		time = 0;
		j=0;
		cout<<endl<<"Node: "<<ckt_nodes[i].node->name<<endl<<endl;
		while(time < tran.tot_t){// && iter <1){
			printf(" %.3e  %.6e\n", time, 
				ckt_nodes[i].value[j]);
			j++;
			time += tran.step_t;
		}
		cout<<"END: "<<ckt_nodes[i].node->name<<endl;
	}
}

void SubCircuit::save_ckt_nodes_to_tr(Tran &tran){
	int index = 0;
	double time = 0;
	int j=0;
	double value = 0;
	for(size_t i=0;i<ckt_nodes.size();i++){
		index = ckt_nodes[i].flag;
		//cout<<"ckt_nodes, index: "<<ckt_nodes[i].node->name
			//<<" "<<ckt_nodes[i].flag<<endl;
		//tran.nodes[index].node = ckt_nodes[i];
		j=0;
		for(time = 0; time < tran.tot_t; time+=tran.step_t){
			value = ckt_nodes[i].value[j];
			//cout<<"time, value: "<<time<<" "<<value<<endl;
			tran.nodes[index].value.push_back(value);
			j++;
		}
	}	
}

void SubCircuit::build_pad_set(){
	pad_set.resize(0);
	//push all pad nodes (LDO nodes);
	for(size_t i=0;i<nodelist.size()-1;i++){
		if(nodelist[i]->isS()==Y){
			Pad *pad_ptr = new Pad();
			pad_ptr->node = nodelist[i];
			pad_set.push_back(pad_ptr);
		}
	}
}

// solve grid with ADI method
/*void SubCircuit::solve_partial_ADI(){
	solve_init();
	double diff1 = 1;
	double diff2 = 1;
	int iter = 0;
	// first solve DC
	while((diff1 > 1e-3 || diff2 > 1e-3)
		&& iter <100){
		// update local grid
		diff1 = solve_ADI_DC();
		// restamp b and update global grid
		//diff2 = solve_global();
		iter ++;
	}
}

// DC, can be parallelized with OpenMP
// solve local grid with LDO as voltage sourses
double SubCircuit::solve_ADI_DC(){
	stringstream sstream;
	Node *nd;
	// solve a line along x direction
	for(double i=lx;i<gx;i++){
		sstream.str("");
		sstream<<"n"<<local_layers[0]<<"_"<<i<<"_"<<ly<<endl;
		nd = get_node_pt(map_node_pt_l, 
			sstream.str());
		// solve along vertical direction
		solve_a_line(nd, NORTH);
	}
	// then solve a line along y direction
	for(double i=ly;i<gy;i++){
		sstream.str("");
		sstream<<"n"<<local_layers[0]<<"_"<<lx<<"_"<<i<<endl;
		nd = get_node_pt(map_node_pt_l, 
			sstream.str());
		// solve along horizontal direction
		solve_a_line(nd, EAST);
	}
}

// solve a line along different dir
void solve_a_line(Node *nd, DIRECTION d){
	vector<double> vec_c;
	vector<double> vec_a;
	vector<double> vec_b;
	vector<double> vec_rhs;
	vector<Node *> vec_nd;

	Net *net;
	Node *nd_temp = nd;
	double c0;
	vec_nd.push_back(nd);
	vec_a.push_back(-1);
	if(nd->is_LDO()){
		vec_b.push_back(1);
		vec_rhs.push_back(nd->value);
	}
	// compute c0 and x0
	if(nd_temp->nbr[d]!=NULL){
		net = nd_temp->nbr[d];
		if(net->type == RESISTOR){
			// find conduc from resis net
			c0 = -1.0/net->value;
		}
	}
	if(nd_temp->nbr[BOTTOM]!=NULL){
		net = nd_temp->nbr[BOTTOM];
		if(net->type== CURRENT)
			b0 = -1.0*net->value;
	}
	// compute c0=c0/b0;
	vec_c.push_back(c0/b0);
	// x0 = x0 / b0;
	nd->value /= b0;	
	while(nd_temp->nbr[d]!= NULL){
		net = nd_temp->nbr[d];
		if(nd_temp == net->ab[0])
			nd_temp = net->ab[1];
		else
			nd_temp = net->ab[0];

	c_temp.clear();
}
*/

// build local VDD nets from LDO
void SubCircuit::build_local_nets(){
	// first delete all voltage nets
	int type = VOLTAGE;
	NetPtrVector &ns = net_set[type];
	for(size_t i=0;i<ns.size();i++){
		delete ns[i];
	}
	ns.clear();
	Node *gnd = NULL;
	for(size_t i=0;i<nodelist.size();i++)
		if(nodelist[i]->is_ground()){
			gnd = nodelist[i];
			break;
		}
	// then create new ones
	Node *nd;
	for(size_t i=0;i<ldolist.size();i++){
		nd = ldolist[i]->A;
		nd->value = VDD_G;
		Net *net = new Net(VOLTAGE, VDD_G, nd, gnd);
		add_net(net);
		// update top nbr net
		nd->rep->nbr[TOP] = net;
		// clog<<"add local net: "<<*net<<endl;
	}
	/*for(int type =0 ;type <NUM_NET_TYPE;type++){
		NetList &ns = net_set[type];
		NetList::iterator it;
		for(it = ns.begin();it!=ns.end();++it)
			cout<<*(*it)<<endl;
	}*/
}

// build global current nets from LDO
// only one time step
void SubCircuit::build_global_nets(){
	// first delete all current nets
	int type = CURRENT;
	NetPtrVector &ns = net_set[type];
	for(size_t i=0;i<ns.size();i++)
		delete ns[i];
	ns.clear();
	// then create new ones
	Node *nd;
	Node *gnd = NULL;
	for(size_t i=0;i<nodelist.size();i++)
		if(nodelist[i]->is_ground()){
			gnd = nodelist[i];
			break;
		}

	double current;
	for(size_t i=0;i<ldolist.size();i++){
		nd = ldolist[i]->nd_in;
		current = ldolist[i]->current;
		//current = 0.3;
		Net *net = new Net(CURRENT, current, nd, gnd);
		add_net(net);
		// update top nbr net
		nd->rep->nbr[BOTTOM] = net;
		// clog<<"add global net: "<<*net<<endl;
	}
}

void SubCircuit::configure_init(){
   cm = &c;
   cholmod_start(cm);
   cm->print = 5;
   size_t n= replist.size();
   b = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
   x = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
   bp = static_cast<double *> (b->x);
   map_node.clear();
   for(size_t i=0;i<nodelist.size()-1;i++)
	   map_node[nodelist[i]->name] = nodelist[i];
}

// update length of b and x
void SubCircuit::reconfigure_TR(){
   size_t n= replist.size();
   b = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
   x = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
   bp = static_cast<double *> (b->x);
}

// stmap matrix and rhs, decomp matrix for DC
void SubCircuit::stamp_decomp_matrix_DC(bool local_flag){
   stamp_by_set(A);
   stamp_rhs_DC(local_flag);
   
   A.set_row(replist.size());
   Algebra::CK_decomp(A, L, cm);
   A.clear();
}

// stamp matrix and rhs, decomp matrix for DC
void SubCircuit::stamp_decomp_matrix_TR(Tran &tran, double time, bool local_flag){
   
   stamp_by_set_tr(A, tran);
   stamp_rhs_tr(local_flag, time, tran); 
   Algebra::CK_decomp(A, L, cm);
   // A.merge();
   // clog<<A<<endl;
   A.clear();
}

// solve eq with decomped matrix
double SubCircuit::solve_CK_with_decomp(){
	// for(size_t i=0;i<replist.size();i++)
		// cout<<"i, bp: "<<i<<" "<<bp[i]<<endl; 
	// solve the eq
	x = cholmod_solve(CHOLMOD_A, L, b, cm);
   	xp = static_cast<double *> (x->x);
	// copy solution to nodes
   	double diff = get_voltages_from_LU_sol(xp);
	// cout<<nodelist<<endl;
	return diff;
}

// solve eq with decomped matrix
double SubCircuit::solve_CK_with_decomp_tr(){
	// modify rhs with Ieq
	modify_rhs_Ieq(bp);
	// solve the eq
	x = cholmod_solve(CHOLMOD_A, L, b, cm);
   	xp = static_cast<double *> (x->x);
	//clog<<"after solve tr. "<<endl;
	//for(size_t i=0;i<replist.size();i++){
		 // clog<<"i, bp: "<<i<<" "<<bnewp[i]<<endl;
		// cout<<"i, bnewp, xp: "<<i<<" "<<bnewp[i]<<" "<<xp[i]<<" "<<*replist[i]<<endl;
	 //}
   	// save_ckt_nodes(tran, xp, time);
	// copy solution to nodes
   	double diff = get_voltages_from_LU_sol(xp);
	// clog<<"diff: "<<diff<<endl;
	// cout<<nodelist<<endl;
	return diff;
}

// update current values for all LDOs
void SubCircuit::update_ldo_current(){
	Node *nd;
	Net *net;
	Node *nb;
	double current=0;
	for(size_t i=0;i<ldolist.size();i++){
		nd = ldolist[i]->A;
		for(DIRECTION d = WEST; d!= UNDEFINED; d = (DIRECTION)(d+1)){
			net = nd->nbr[d];
			if(net == NULL) continue;
			// no current net in ldo 
			// locations in input file
			if(net->type != RESISTOR) 
				continue;
			// work on resistance net
			if(net->ab[0] == nd)
				nb = net->ab[1];
			else
				nb = net->ab[0];
			current += (nd->value - nb->value ) / net->value;
			// clog<<"net, nd, nb, current: "<<*net<<" "<<*nd<<" "<<*nb<<current<<endl;
		}
		// copy old current
		//ldolist[i]->current_old = 
			// ldolist[i]->current;
		// assign new current
		ldolist[i]->current = current;
		// also update the current net in ckt_g
		Net *net_g = ldolist[i]->nd_in->nbr[BOTTOM];
		if(net_g->type == CURRENT)
			net_g->value = current;
		// clog<<"ldo old / new current: "<<ldolist[i]->current_old<<" "<<current<<endl;
	}
}

/*// update ldo input voltage
void SubCircuit::update_ldo_vin(){
	for(size_t i=0;i<ldolist.size();i++){
		ldolist[i]->vin = ldolist[i]->nd_in->value;
	}
}*/

// modify rhs with new current value of LDO
void SubCircuit::modify_ldo_rhs(){
	Node *nd;
	size_t rid;
	for(size_t i=0;i<ldolist.size();i++){
		nd = ldolist[i]->nd_in;
		//clog<<"nd: "<<*nd<<endl;
		rid = nd->rid;
		bp[rid] = -ldolist[i]->current;
		// clog<<"new bp: "<<rid<<" "<<bp[rid]<<endl;
	}
}

double SubCircuit::locate_maxIRdrop(){
	max_IRdrop = 0;
	for(size_t i=0;i<replist.size();i++){	
		double IR_drop = VDD_G - replist[i]->value;		
		if(IR_drop > max_IRdrop)
			max_IRdrop = IR_drop;
	}
	return max_IRdrop;
}

// 1. adjust ldo locations
// 2. update the ldo correlated nets
void SubCircuit::relocate_pads(){
	// record origin_pad_set before optimization
	origin_pad_set.resize(pad_set.size());	
	assign_pad_set(pad_set, origin_pad_set);

	vector<double> ref_drop_vec;
	//double min_IR = max_IRdrop;	
	// clog<<"min_IR initial is: "<<min_IR<<endl;
	// for local pad movement
	// find control nodes for each pad
	extract_pads(pad_set);
	// find the tune spot for control nodes	
	update_pad_control_nodes(pad_set);	

	// find new point for all pads	
	update_pad_pos_all(pad_set);

	clog<<"before graph move. "<<endl;
	// move pads according to graph contraints
	graph_move_pads();
	clog<<"after graph move. "<<endl;

	//clog<<"after graph move. "<<endl;
	clear_flags();

	Node *rm_node = origin_pad_set[0]->rep;
	Node *add_node = pad_set[0]->node->rep;
	rebuild_local_nets(rm_node, add_node);
		
	ref_drop_vec.clear();
}

// update the old pad set value
void SubCircuit::assign_pad_set(vector<Pad*> pad_set, vector<Node*>&pad_set_old){
	//clog<<"assign pad set."<<endl;
	// clog<<"size: "<<pad_set_old.size()<<endl;
	// clog<<"pad_set size: "<<pad_set.size()<<endl;
	pad_set_old.resize(pad_set.size());
	for(size_t i=0;i<pad_set_old.size();i++){
		if(pad_set[i]->node == NULL)
			clog<<"NULL node. "<<endl;
		pad_set_old[i] = pad_set[i]->node;
		//if(pad_set[i]->node->get_layer()==local_layers[0])
		// cout<<"pad: "<<i<<" "<<*pad_set_old[i]<<endl;	
	}
}

void SubCircuit::build_pad_graph(vector<Pad*> &pad_set){
	Pad *pad;
	Pad *pad_nbr;
	bool flag_pad = false;
	bool flag_nbr = false;
	// clear content
	for(size_t i=0;i<pad_set.size();i++){
		pad_set[i]->nbrs.clear();
	}
	// find nbr pad nodes
	for(size_t i=0;i<pad_set.size();i++){
		//clog<<"pad: "<<*pad->node<<endl;
		flag_pad = false;
		flag_nbr = false;
		pad = pad_set[i];
		pad_nbr = find_nbr_pad(pad_set, pad);
		// cout<<"pad, nbr: "<<*pad->node<<" "<<*pad_nbr->node<<endl;
		for(size_t j=0;j<pad_nbr->nbrs.size();j++){
			if(pad_nbr->nbrs[j]->node->name== pad->node->name)
				flag_pad = true;
				break;
		}
		for(size_t j=0;j<pad->nbrs.size();j++){
			if(pad->nbrs[j]->node->name== pad_nbr->node->name)
				flag_nbr = true;
				break;
		}
		if(flag_pad == false){
			pad->nbrs.push_back(pad_nbr);
		}
		if(flag_nbr == false)
			pad_nbr->nbrs.push_back(pad);
	}
	/*for(size_t i=0;i<pad_set.size();i++){
		Pad *pad = pad_set[i];
		cout<<"pad: "<<*pad->node<<endl;
		for(size_t j=0;j<pad->nbrs.size();j++){
			cout<<"nbr: "<<*pad->nbrs[j]->node<<endl;
		}
	}*/
}

// use Euclidiean distance to locate nearest nbr pad
Pad * SubCircuit::find_nbr_pad(vector<Pad*> &pad_set, Pad *pad){
	Pad * nbr;
	double distance=-1;
	double min_dist=0;
	bool flag = false;
	size_t min_index=0;
	for(size_t i=0;i<pad_set.size();i++){
		nbr = pad_set[i];
		// need to make sure they are in same layer
		if(nbr->node->pt.z != pad->node->pt.z)
			continue;
		if(nbr->node->name == pad->node->name)
			continue;
		distance = get_distance(nbr->node, pad->node);
		if(flag == false){
			flag = true;
			min_dist = distance;
			min_index = i;
		}else{
			if(distance < min_dist){
				min_dist = distance;	
				min_index = i;
			}	
		}
	}
	return pad_set[min_index];
}

double SubCircuit::get_distance(Node *na, Node *nb){
	double distance = 0;
	double delta_x = 0;
	double delta_y = 0;
	delta_x=(na->pt.x-nb->pt.x);
	delta_y=(na->pt.y-nb->pt.y);
	delta_x *= delta_x;
	delta_y *= delta_y;
	distance = sqrt(delta_x + delta_y);
	return distance;
}

// find control nodes for each pad
void SubCircuit::extract_pads(vector<Pad*> pad_set){
	int pad_number = 1;
	vector<Node*> pair_first;
	vector<double> pair_second;
	pair<Node*, double> pair_nd;
	//int pad_number = 5;
	double distance = 0;
	map<Node*, double>::iterator it;

	clear_pad_control_nodes(pad_set);
	for(size_t i=0;i<replist.size();i++){
		int count = 0;
		pair_first.clear();
		pair_second.clear();
		Node *nd = replist[i];
		// search for closest pads
		for(size_t j=0;j<pad_set.size();j++){
			Node *ptr = pad_set[j]->node;
			// make sure they are in same layer
			if(ptr->get_layer() != nd->get_layer())
				continue;
			distance = get_distance(ptr, nd);

			if(count < pad_number){
				pair_first.push_back(ptr);
				pair_second.push_back(distance);
				count++;
			}else{// substitute the pad node
				double max_dist = 0;
				size_t max_index = 0;
				for(size_t k=0;k<pair_second.size();k++){
					if(pair_second[k]>max_dist){
						max_dist = pair_second[k];
						max_index = k;
					}
				}
				if(distance < max_dist){ 
					pair_first[max_index] = ptr;
					pair_second[max_index] = distance;
				}
			}
		}
		// then map these distance into pads
		for(size_t j=0;j<pair_first.size();j++){
			Node *ptr = pair_first[j];
			//if(nd->name == "n0_0_0")
			//clog<<"ptr: "<<*ptr<<endl;
			for(size_t k=0;k<pad_set.size();k++){
				if(pad_set[k]->node->name == ptr->name){
					// control nodes
					pair_nd.first = nd;
					// distance
					pair_nd.second = pair_second[j];
					pad_set[k]->control_nodes.insert(pair_nd);
					break;
				}
			}
		}
	}
	//print_pad_map();	
	pair_first.clear();
	pair_second.clear();
	//print_pad_map();
}

void SubCircuit::clear_pad_control_nodes(vector<Pad*> &pad_set){
	for(size_t i=0;i<pad_set.size();i++){
		pad_set[i]->control_nodes.clear();
	}
}

// tune 50% nodes with the small IR drops
void SubCircuit::update_pad_control_nodes(vector<Pad*> pad_set){
	for(size_t i=0;i<pad_set.size();i++){
		if(pad_set[i]->control_nodes.size()==0)
			continue;
		//clog<<"before middle value. "<<endl;	
		double middle_value = locate_ref(pad_set, i);
		//clog<<"middle value: "<<middle_value<<endl;
		pad_set[i]->ref_vol = middle_value;
		//clog<<"after assign. "<<endl;
		//cout<<"middle value: "<<middle_value<<endl;
	}
}

bool compare_values(double a, double b){
	return (a<b);
}

// locate the tune spot for the control nodes.
double SubCircuit::locate_ref(vector<Pad*> pad_set, size_t i){
	Pad *pad_ptr;
	map<Node*, double>::iterator it;
	Node *nd;
	double weight = 0;
	//vector<double> drop_vec;
	pad_ptr = pad_set[i];
	// clog<<"get pad_ptr. "<<endl;
	pad_ptr->drop_vec.clear();
	// cout<<"pad: "<<*pad<<endl;
	for(it = pad_ptr->control_nodes.begin();
			it != pad_ptr->control_nodes.end();
			it++){
		nd = it->first;
		// cout<<"control: "<<*nd<<endl;
		// need to be generated with worst_cur
		weight = nd->value;
		if(weight <0)
			weight *=10;
		// clog<<"before assign weight. "<<endl;
		pad_ptr->control_nodes[nd] = weight;
		// clog<<"before push drop vec. "<<endl;
		pad_ptr->drop_vec.push_back(nd->value); 
	}
	sort(pad_ptr->drop_vec.begin(), pad_ptr->drop_vec.end(),
			compare_values);
	pad_ptr->ratio = 2;
	size_t id = pad_ptr->drop_vec.size() / pad_ptr->ratio;
	double middle_value = pad_ptr->drop_vec[id];
	//drop_vec.clear();
	return middle_value;
}

double SubCircuit::calc_avg_ref_drop(vector<double> &ref_drop_vec){
	Node *pad;
	Pad *pad_ptr;
	double max_drop, min_drop;
	//double sum_max = 0;
	//double sum_min = 0;
	double sum_diff = 0;
	size_t count = 0;
	double ref_drop_value = 0;

	for(size_t i=0;i<pad_set.size();i++){
		if(pad_set[i]->control_nodes.size()==0)
			continue;
		count ++;	
		ref_drop_value = ref_drop_vec[i];
		map<Node *, double>::iterator it;
		pad_ptr = pad_set[i];
		pad = pad_ptr->node;
		max_drop = 0;
		min_drop = -1;
		
		for(it = pad_ptr->control_nodes.begin();
		    it != pad_ptr->control_nodes.end();
		    it++){
			if(it->second > ref_drop_value)
				continue;
			  if(it->second > max_drop)
				max_drop = it->second;
			  if(min_drop == -1)
				min_drop = it->second;
			  else if(it->second < min_drop)
				min_drop = it->second;
		}
		pad_ptr->data = max_drop - min_drop;
		sum_diff += pad_ptr->data;
	}
	double avg_drop = sum_diff / count;
	return avg_drop;
}

// expand from (x,y) to nearest node in grid
// fits for non-uniform grid
Node * SubCircuit::pad_projection( 
	Pad *pad, Node *nd){
	//Node *nd;
	stringstream sstream;
	string pt_name;
	Node *nd_new=NULL;
   	Point pt;
	//int gap - 10;

	LDO *ldo = NULL;
	for(size_t i=0;i<ldolist.size();i++){
		if(ldolist[i]->A->name == nd->name){
			ldo = ldolist[i];
		}
	}
	//nd = pad->node;
	pt.z = nd->get_layer();
	pt.x = pad->newx;
	pt.y = pad->newy;
	if(pad->newx == nd->pt.x && 
		pad->newy == nd->pt.y){
		nd_new = nd;	
		return nd_new;
	}

	sstream<<"n"<<pt.z<<"_"<<pt.x<<"_"<<pt.y; 
	name = sstream.str();
	clog<<"pt_name: "<<name<<endl;
	// first see if this node is on grid
	// and if it is occupied by pad or not
	nd_new = get_node(name);
	if(nd_new == NULL){
		// clog<<"null node: "<<endl;
		return nd;
	}
	nd_new = nd_new->rep;
	// if this node is not occupied by pad
	// need to adjust the local pads
	Node *nb = project_local_pad(nd, nd_new, ldo);
	if(nb->name == nd->name)
		return nd;
	else
		return nb;
}

// project local pad, setting ldo into new white spaces near current/device blocks
Node * SubCircuit::project_local_pad(Node *nd, Node *nd_new, LDO *ldo){
	Node *nd_new_ldo;
	int ref_x = nd_new->pt.x;
	int ref_y = nd_new->pt.y;
	int ref_dx, ref_dy;
	double ref_dist;
	// if min_id == -1, stay on old wspace
	LDO ldo_ptr;
	// find the reference distance
	ref_dx = fabs(ref_x - nd->pt.x);
	ref_dy = fabs(ref_y - nd->pt.y);
	ref_dist = sqrt(ref_dx * ref_dx + ref_dy * ref_dy);

	// ldo_ptr points to current LDO
	for(size_t i=0;i<ldolist.size();i++){
		Node *A = ldolist[i]->A;
		if(A->name == nd->name)
			ldo_ptr = *ldolist[i];
	}
	// 1. project to nearest LDO candi nodes
	// update ref_x and ref_y
	project_ldo_node(ref_x, ref_y, ldo_ptr);
	clog<<"projected ref_x and y: "<<ref_x<<" "<<ref_y<<endl;

	// then start to expand to find candi LDO nodes
	nd_new_ldo = expand_ldo_location(ref_dist, 
			ref_x, ref_y, ldo_ptr);

	ldo->A = nd_new_ldo;
	// also need to change nd_in of LDO in ckt_g
	return nd_new_ldo; 
}

double SubCircuit::update_pad_pos_all(vector<Pad*> pad_set){
	double total_dist = 0;
	double dist = 0;
	for(size_t i=0;i<pad_set.size();i++){
		if(pad_set[i]->control_nodes.size()==0)
			continue;

		double ref_drop_value = pad_set[i]->ref_vol;

		dist = update_pad_pos(ref_drop_value, i);
		total_dist += dist;
	}
	return total_dist;
}

// decide pad's new pos with the weights
// need to be tuned
double SubCircuit::update_pad_pos(double ref_drop_value, size_t i){
	Node *pad;
	Pad *pad_ptr;
	Node *nd;
	double weight = 0;
	double pad_newx;
	double pad_newy;
	map<Node *, double>::iterator it;

	double sum_weight = 0;
	double weighted_x =0;
	double weighted_y =0;
	pad_ptr = pad_set[i];
	pad = pad_ptr->node;
	for(it = pad_ptr->control_nodes.begin();
			it != pad_ptr->control_nodes.end();
			it++){
		if(it->second > ref_drop_value)
			continue;

		nd = it->first;
		weight = 1.0/it->second;
		weighted_x += weight * nd->pt.x;
		weighted_y += weight * nd->pt.y;
		sum_weight += weight; 	
	}

	if(sum_weight !=0){
		pad_newx = weighted_x / sum_weight;
		pad_newy = weighted_y / sum_weight;

		round_data(pad_newx);
		round_data(pad_newy);
	
		pad_ptr->newx = pad_newx;
		pad_ptr->newy = pad_newy;
	}else{
		pad_ptr->newx = pad->pt.x;
		pad_ptr->newy = pad->pt.y;
	}

	double dist = sqrt(weighted_x*weighted_x 			 + weighted_y*weighted_y);

	return dist;
}

// round double, depends on whether the frac >=0.5
void SubCircuit::round_data(double &data){
	double fractpart, intpart;
	fractpart = modf(data, &intpart);
	if(fractpart >= 0.5)
		data = ceil(data);
	else
		data = floor(data);
}

// project the wspace reference node to nearest ldo
// candi locations
void SubCircuit::project_ldo_node(int &ref_x, int &ref_y, LDO &ldo){
	int width = ldo.width;
	int height = ldo.height;
	
	ref_x = ref_x / width;
	ref_x *= width;
	//x[1] = x[0] + width;

	ref_y = ref_y / height;
	ref_y *= height;
}

double SubCircuit::calc_avg_ref(vector<Pad*> &pad_set, vector<double> ref_drop_vec){
	//Node *pad;
	//Pad *pad_ptr;
	double sum_ref = 0;
	size_t count = 0;
	double ref_drop_value = 0;
	double avg_ref = 0;

	for(size_t i=0;i<pad_set.size();i++){
		if(pad_set[i]->control_nodes.size()==0)
			continue;

		count ++;	
		ref_drop_value = ref_drop_vec[i];	
		sum_ref += ref_drop_value;
	}
	avg_ref = sum_ref / count;
	return avg_ref;
}

void SubCircuit::graph_move_pads(){
	Node *new_pad;
	int id=0;
	// for(size_t i=0;i<5;i++){
	// do{
		// id = locate_max_drop_pad(pad_set);
		clog<<"id: "<<id<<endl;
		// if(id==-1) break;
		Pad *pad_ptr = pad_set[id];
		Pad *pad_nbr = NULL;
		Node *pad = pad_ptr->node;
		clog<<"old_pad: "<<*pad<<endl;
		new_pad = pad_projection(pad_ptr, pad);
		clog<<" old pad / new pad: "<<*pad<<" "<<*new_pad<<endl;

		pad_ptr->visit_flag = true;
		for(size_t i=0;i<pad_ptr->nbrs.size();i++){
			pad_nbr = pad_ptr->nbrs[i];
			if(pad_nbr->fix_flag == false){
				pad_nbr->fix_flag = true;
				//cout<<"pad_nbr: "<<*pad_nbr->node<<endl;
			}
		}
		// assign later, after creating X node
		pad_ptr->node = new_pad;
		pad_ptr->control_nodes.clear();
	// }while(id != -1);
}

// locate id that has minimum value and not visited or fixed 
int SubCircuit::locate_max_drop_pad(vector<Pad*> pad_set){
	int min_id = -1;
	double min_ref = 0;
	bool flag = false;
	for(size_t i=0;i<pad_set.size();i++){
		if(pad_set[i]->control_nodes.size()==0)
			continue;
		if(pad_set[i]->visit_flag == true ||
			pad_set[i]->fix_flag ==  true)
			continue;
		//clog<<"i, vec: "<<i<<" "<<vec[i]<<endl;
		if(flag == false){
			flag = true;
			min_ref = pad_set[i]->ref_vol;
			min_id = i;
		}
		else if(pad_set[i]->ref_vol < min_ref){
			min_ref = pad_set[i]->ref_vol;
			min_id = i;
		}
	}	
	//clog<<"min_ref, min_pad: "<<min_ref<<" "<<*pad_set[min_id]->node<<endl;
	return min_id;
}

void SubCircuit::clear_flags(){
	Node *nd;
	for(size_t i=0;i<candi_pad_set.size();i++){
		nd = candi_pad_set[i]->node;
		nd->flag_visited = false;
	}
}

// mark the grid with block and LDO occupation
void SubCircuit::mark_geo_occupation(){
	int width=0;
	int height = 0;

	if(ldolist.size()!=0){
		width = ldolist[0]->width;
		height = ldolist[0]->height;
	}

	Node *nd;
	for(size_t i=0;i<nodelist.size()-1;i++){
		nd = nodelist[i];
		if(nd->isS()==Z) continue;
		int x = nd->pt.x;
		int y = nd->pt.y;
		// mark node with geo occupation
		if(!(x % width ==0 && y % height ==0))
			continue;
		// clog<<"start mark node: "<<*nd<<endl;
		bool flag_module = false;
		// check if it is module
		for(size_t j=0;j<wspacelist.size();j++){
			int module_xl = wspacelist[j]->node[0]->x;
			int module_xr = wspacelist[j]->node[0]->x + wspacelist[j]->width;
			int module_yl = wspacelist[j]->node[0]->y;
			int module_yr = wspacelist[j]->node[0]->y + wspacelist[j]->height;
			// clog<<"bxl, bxr, byl, byr: "<<module_xl<<" "<<module_xr<<" "<<module_yl<<" "<<module_yr<<endl;
			if(x>=module_xl && x <= module_xr && y>=module_yl && y <= module_yr){
				// clog<<"mark block. "<<endl<<endl;
				// in module
				nd->assign_geo_flag(SBLOCK);
				flag_module = true;
				break;
			}
		}
		if(flag_module == true)
			continue;
		// check if it is LDO
		bool flag_ldo = false;
		for(size_t j=0;j<ldolist.size();j++){
			int ldo_xl = ldolist[j]->A->pt.x;
			int ldo_xr = ldolist[j]->A->pt.x + width;
			int ldo_yl = ldolist[j]->A->pt.y;
			int ldo_yr = ldolist[j]->A->pt.y + height;
			// clog<<"lxl, lxr, lyl, lyr: "<<ldo_xl<<" "<<ldo_xr<<" "<<ldo_yl<<" "<<ldo_yr<<endl;
			if(x >= ldo_xl && x <= ldo_xr && y >= ldo_yl && y <= ldo_yr){
				// clog<<"mark ldo. "<<endl<<endl;
				// in module
				nd->assign_geo_flag(SLDO);
				Pad *pad_ptr = new Pad();
				pad_ptr->node = nd;
				candi_pad_set.push_back(pad_ptr);
				flag_ldo = true;
				break;
			}
		}
		if(flag_ldo == true)
			continue;
		// else assign blank
		// clog<<"mark blank. "<<endl<<endl;
		nd->assign_geo_flag(SBLANK);
		Pad *pad_ptr = new Pad();
		pad_ptr->node = nd;
		candi_pad_set.push_back(pad_ptr);
	}
	// count is the maximum candidate number for LDO
	MAX_NUM_LDO = candi_pad_set.size();
	/* cout<<"start to output candi_pad set. "<<endl<<endl<<endl;
	for(size_t i=0;i<MAX_NUM_LDO;i++){
		Pad *pad = candi_pad_set[i];
		Node *nd = pad->node;
		if(nd != NULL)
			cout<<"nd: "<<*nd<<endl;
	}*/
}

bool SubCircuit::node_in_ldo_or_block(double x, double y){
	int width = ldolist[0]->width;
	int height = ldolist[0]->height;

	bool flag_module = false;
	// check if it is module
	for(size_t j=0;j<wspacelist.size();j++){
		int module_xl = wspacelist[j]->node[0]->x;
		int module_xr = wspacelist[j]->node[0]->x + wspacelist[j]->width;
		int module_yl = wspacelist[j]->node[0]->y;
		int module_yr = wspacelist[j]->node[0]->y + wspacelist[j]->height;
		if(x>=module_xl && x <= module_xr && y>=module_yl && y <= module_yr){
			// in module
			flag_module = true;
			break;
		}
	}
	if(flag_module == true)
		return true;
	// check if it is LDO
	bool flag_ldo = false;
	for(size_t j=0;j<ldolist.size();j++){
		int ldo_xl = ldolist[j]->A->pt.x;
		int ldo_xr = ldolist[j]->A->pt.x + width;
		int ldo_yl = ldolist[j]->A->pt.y;
		int ldo_yr = ldolist[j]->A->pt.y + height;
		if(x >= ldo_xl && x <= ldo_xr && y >= ldo_yl && y <= ldo_yr){
			// in module
			flag_ldo = true;
			break;
		}
	}
	if(flag_ldo == true)
		return true;
	return false;
}

// expand to find the location of LDO around ref_x, ref_y
// if the location is farther than ref_dist, restore back
Node * SubCircuit::expand_ldo_location(double ref_dist, int ref_x, int ref_y, LDO &ldo_ptr){
	Node *nd;
	double diff_x = 0;
	double diff_y = 0;
	double dist = 0;
	double min_dist = 0;
	Node *min_nd = NULL;
	bool flag = false;
	for(size_t i=0;i<replist.size();i++){
		nd = replist[i];
		if(nd->get_geo_flag() != SBLANK)
			continue;
		if(nd->isS()==Y || nd->isS() ==Z) 
			continue;
		diff_x = nd->pt.x - ref_x;
		diff_y = nd->pt.y - ref_y;
		dist = sqrt(diff_x*diff_x + diff_y*diff_y);
		if(dist < ref_dist){
			if(flag == false){
				flag  = true;
				min_dist = dist;
				min_nd = nd;
			}else if(dist < min_dist){
				min_dist = dist;
				min_nd = nd;
			}
		}
	}
	// if no candidates, return original node
	if(min_nd == NULL)
		return ldo_ptr.A;
	clog<<"min_dist, min_nd: "<<min_dist<<" "<<*min_nd<<" "<<endl;
	// else return candidate node
	return nd;
}

// perform LDO node change
// modify nets and nodes with the change of LDO
void SubCircuit::rebuild_local_nets(Node *rm_node, Node *add_node){
		if(rm_node->name == add_node->name)
			return;
		 // clog<<"rm_nod, add_node: "<<*rm_node<<" "<<*add_node<<endl;	
		// modify node info
		rm_node->disableY();
		add_node->enableY();
		
		if(rm_node->is_LDO()){
			add_node->enableLDO();
			rm_node->disableLDO();
		}
		add_node->value = rm_node->value;
		rm_node->value = 0;
		add_node->rep = add_node;

		// clog<<"rm_node: "<<*rm_node<<" "<<rm_node->is_LDO()<<" "<<rm_node->isS()<<endl;

		// clog<<"add_node: "<<*add_node<<" "<<add_node->is_LDO()<<" "<<add_node->isS()<<endl;
		Net *net = rm_node->nbr[TOP];
		// should print error info
		if(net == NULL) return;

		// modify net info
		if(net->ab[0]->is_ground())
			net->ab[1] = add_node;
		else if(net->ab[1]->is_ground())
			net->ab[0] = add_node;
		ldolist[0]->A = add_node;
		// clog<<"ldolist[0]->A: "<<*ldolist[0]->A<<endl;
		rm_node->nbr[TOP] = NULL;
		add_node->nbr[TOP] = net;
	//}
}
	
// for global circuit
// rebuild all the RLV nets from ldolist
void SubCircuit::rebuild_global_nets(){	
	int type = CURRENT;
	// then delete all the nets in global circuit
	NetPtrVector & ns = net_set[type];
	Net *net = ns[0];
	Node *nd = net->ab[0];
	if(nd->is_ground())
		nd = ns[0]->ab[1];
	nd->nbr[BOTTOM] = NULL;
	nd = ldolist[0]->nd_in;
	nd->nbr[BOTTOM] = net;
	nd->value = ldolist[0]->current;
}

void SubCircuit::extract_node(char * str, Node & nd){
	//static Node gnd(string("0"), Point(-1,-1,-1));
	if( str[0] == '0' ) {
		nd.name="0";
		nd.pt.set(-1,-1,-1);
		return;
	}

	long z, y, x;
	int flag = -1;
	char * chs;
	char * saveptr;
	char l[MAX_BUF];
	strcpy(l, str);
	const char * sep = "_n";
	chs = strtok_r(l, sep, &saveptr); // initialize
	// for transient, 'Y' is the VDD source node
	if( chs[0] == 'X' || chs[0]== 'Y' || chs[0]== 'Z' ){
		flag = chs[0]-'X';
		chs = strtok_r(NULL, sep, &saveptr);
	}
	z = atol(chs);
	chs = strtok_r(NULL, sep, &saveptr);
	x = atol(chs);
	chs = strtok_r(NULL, sep, &saveptr);
	y = atol(chs);

	nd.name.assign(str);
	nd.pt.set(x,y,z);
	nd.flag = flag;
	//return Node(string(str), Point(x,y,z), flag);
}

// Given a net with its two nodes, update the connection information
// for thet two nodes
void SubCircuit::update_node(Net * net){
	// first identify their connection type:
	// 1. horizontal/vertical   2. via/VDD 3. current
	//
	// swap *a and *b so that a is:
	// WEST   for horizontal
	// SOUTH  for vertical
	// BOTTOM for via / XVDD
	// ground node for CURRENT
	Node *a=net->ab[0], *b=net->ab[1];
	//cout<<"setting "<<net->name<<" nd1="<<nd1->name<<" nd2="<<nd2->name<<endl;
	
	if(net->type == CAPACITANCE){
		// make sure a is the Z node, b is ground
		if(a->isS() != Z) swap<Node*>(a,b);
		// only needs single dir nbr net for index
		// TOP for resistance
		// BOTTOM for capacitance
		a->set_nbr(BOTTOM, net);
	}
	else if(net->type == INDUCTANCE){
		// a is Y, b is X
		if(a->isS()== X) swap<Node*>(a,b);
		a->set_nbr(BOTTOM, net);
		b->set_nbr(TOP, net);
	}
	// resistance type net with special nodes
	else if(net->type == RESISTOR && a->isS()!= -1 
		&& b->isS() == -1){
		if(a->isS()== X){
			a->set_nbr(BOTTOM, net);
			b->set_nbr(TOP, net);
		}
		else if(a->isS() == Z){
			// bottom for z node has been taken by 
			// current
			a->set_nbr(TOP, net);
		}
	}
	else if(net->type == RESISTOR && b->isS()!= -1 && a->isS() ==-1){
		if(b->isS()== X){
			b->set_nbr(BOTTOM, net);
			a->set_nbr(TOP, net);
		}
		else if(b->isS() == Z){
			b->set_nbr(TOP, net);
		}
	}
	else if( a->get_layer() == b->get_layer() ){
		// horizontal or vertical resistor in the same layer
		int layer = a->get_layer();
		if( a->pt.y == b->pt.y ){// horizontal
			if(a->pt.x > b->pt.x) swap<Node*>(a,b);
			a->set_nbr(EAST, net);
			b->set_nbr(WEST, net);
			layer_dir[layer] = HR;
		}
		else if( a->pt.x == b->pt.x ){// vertical
			if(a->pt.y > b->pt.y) swap<Node*>(a,b);
			a->set_nbr(NORTH, net);
			b->set_nbr(SOUTH, net);
			layer_dir[layer] = VT;
		}
		else{
			clog<<*net<<endl;
			report_exit("Diagonal net\n");
		}
	}
	else if( //fzero(net->value) && 
		 !a->is_ground() &&
		 !b->is_ground() ){// this is Via (Voltage or Resistor )
		if( a->get_layer() > b->get_layer() ) swap<Node*>(a,b);
		a->set_nbr(TOP, net);
		b->set_nbr(BOTTOM, net);
	}
	else if (net->type == VOLTAGE){// Vdd Voltage
		// one is X node, one is ground node
		// Let a be X node, b be another
		if( a->is_ground() ) swap<Node*>(a,b);
		a->flag = Y;		// set a to be X node
		a->set_nbr(TOP, net);	// X -- VDD -- Ground
		a->set_value(net->value);
	}
	
	else{// if( net->type == CURRENT ){// current source
		// let a be ground node
		if( !a->is_ground() ) swap<Node*>(a,b);
		b->set_nbr(BOTTOM, net);
	}
}

void SubCircuit::create_current_LDO_graph(){
	// 1. build candidate pad graph
	build_pad_graph(candi_pad_set);
	// update the flag for the single LDO
	update_single_pad_flag(pad_set[0]);
	// 2. search for control nodes for all candi
	extract_pads(candi_pad_set);
	// 3. get ref_value as IR drop for each candi
	update_pad_control_nodes(candi_pad_set);
}

void SubCircuit::extract_add_LDO_dc_info(vector<Pad*> & LDO_pad_vec){	
	int iter = 0;
	double THRES = VDD_G * 0.1;
	// for(size_t i=0;i<ldolist.size();i++)
		//clog<<"current ldo: "<<i<<" "<<*ldolist[i]->A<<endl;
	// while not satisfied and still have room,
	// perform optimization
	// while(max_IRdrop > THRES && )
	while(ldolist.size() < MAX_NUM_LDO &&
			iter ++ <1){
		// 4. LDO should go to candi with 
		// maximum IR
		Pad *pad_ptr = locate_candi_pad_maxIR(candi_pad_set);
		// clog<<"new location for LDO: "<<*pad_ptr->node<<endl;
		LDO_pad_vec.push_back(pad_ptr);
		// 5. update the nbr flags for 
		// candi in graph
		update_single_pad_flag(pad_ptr);
	}
}

// return pad candi with max IR (still available candi)
Pad* SubCircuit::locate_candi_pad_maxIR(vector<Pad*> pad_set){
	double min_vol = VDD_G;
	Pad *pad_ptr = NULL;
	bool flag = false;
	// for(size_t i=0;i<pad_set.size();i++)
		// clog<<"i, pad_set_flag: "<<i<<" "<<*pad_set[i]->node<<" "<<pad_set[i]->node->flag_visited<<endl;
	for(size_t i=0;i<pad_set.size();i++){
		if(pad_set[i]->node->flag_visited == true)
			continue;
		// skip the one that already has pad
		if(pad_set[i]->node->isS()==Y)
			 continue;
		if(flag == false){
			min_vol = pad_set[i]->ref_vol;
			pad_ptr = pad_set[i];
			// clog<<"min, pad: "<<min_vol<<" "<<pad_set[i]->ref_vol<<" "<<*pad_ptr->node<<endl;
			flag = true;
		}
		else if(pad_set[i]->ref_vol < min_vol){
			min_vol = pad_set[i]->ref_vol;
			pad_ptr = pad_set[i];
			// clog<<"min, pad: "<<min_vol<<" "<<pad_set[i]->ref_vol<<" "<<*pad_ptr->node<<endl;
		}
	}
	// clog<<"min_vol, pad: "<<min_vol<<" "<<*pad_ptr->node<<endl;
	return pad_ptr; 
}

// mark pad flag_visited in pad graph
void SubCircuit::update_single_pad_flag(Pad* pad){
	//Pad *pad_nbr = NULL;
	Node *nd;
	Node *nd_nbr;

	nd = pad->node;
	//clog<<"pad node: "<<*nd<<endl;
	nd->flag_visited = true;
	// mark nbr pads with visited
	for(size_t j=0;j<pad->nbrs.size();j++){
		nd_nbr = pad->nbrs[j]->node;
		nd_nbr->flag_visited = true;
		// clog<<"mark nbr: "<<*nd_nbr<<endl;
	}
}

// create local current nets for new LDO
void SubCircuit::create_local_LDO_new_nets(vector<Pad*> LDO_pad_vec){
	Node *nd;
	Node *gnd =NULL;
	for(size_t i=0;i<nodelist.size();i++)
		if(nodelist[i]->is_ground()){
			gnd = nodelist[i];
			break;
		}

	for(size_t i=0;i<LDO_pad_vec.size();i++){
		nd = LDO_pad_vec[i]->node;
		nd->enableY();
		nd->enableLDO();
		nd->value = VDD_G;
		Net *net = new Net(VOLTAGE, VDD_G, nd, gnd);
		add_net(net);
		// update top nbr net
		nd->rep->nbr[TOP] = net;
	}
}

// need to create all the RLVI nets
void SubCircuit::create_global_LDO_new_nets(vector<Pad*> LDO_pad_vec){
	Node *nd, *nd_in;
	double current_value = 0;
	// Net *net;
	Node *gnd = NULL;
	for(size_t i=0;i<nodelist.size();i++)
		if(nodelist[i]->is_ground()){
			gnd = nodelist[i];
			break;
		}

	stringstream sstream;
	// Node *gnd = nodelist[nodelist.size()-1];
	NET_TYPE net_type;
	long ldo_z = ldolist[0]->nd_in->pt.z;
	for(size_t i=0;i<LDO_pad_vec.size();i++){
		nd = LDO_pad_vec[i]->node;
		sstream.str("");
		sstream<<"n"<<ldo_z<<"_"<<nd->pt.x<<"_"<<nd->pt.y;
		nd_in = get_node(sstream.str());
		// clog<<"nd:"<<*nd<<" nd_in: "<<*nd_in<<endl;
		// update_node(net_resis);
		// then create current net 
		net_type = CURRENT;
		Net *net_current = new Net(net_type, 
			current_value, nd_in, gnd);
		add_net(net_current);

		// clog<<"new cur net: "<<*net_current<<endl;
		nd_in->nbr[BOTTOM] = net_current;
		//update_node(net_current);
	}
}

// only modify vol net values for LDO
void SubCircuit::modify_local_nets(){
	Node *nd;
	Net *net;
	for(size_t i=0;i<ldolist.size();i++){
		nd = ldolist[i]->A;
		net = nd->nbr[TOP];
		// should report error
		if(net == NULL) continue;
		net->value = nd->value;
		// clog<<"local nets: "<<*net<<endl;
	}
}

// only modify current net values for LDO
void SubCircuit::modify_global_nets(){
	Node *nd;
	Net *net;
	for(size_t i=0;i<ldolist.size();i++){
		nd = ldolist[i]->nd_in;
		net = nd->nbr[BOTTOM];
		// should report error
		if(net == NULL){
			clog<<"null net: "<<*nd<<endl;
			continue;
		}
		net->value = ldolist[i]->current;	
		// clog<<"global nets: "<<*net<<endl;
	}
}

void SubCircuit::modify_rhs_Ieq(double *rhs){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		if(type ==CAPACITANCE){	
			for(size_t i=0;i<ns.size();i++)
				modify_rhs_Ieq_c(ns[i], rhs);
		}
		else if(type == INDUCTANCE){
			for(size_t i=0;i<ns.size();i++){
				modify_rhs_Ieq_l(ns[i], rhs);	
			}
		}
	}
}

void SubCircuit::modify_rhs_Ieq_c(Net *net, double *rhs){
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
        // nk point to Z node
        if(nk->isS() != Z){
		swap<Node *>(nk, nl);
		swap<Node*>(net->ab[0], net->ab[1]);
	}
	//clog<<"nk, nl: "<<*nk<<" "<<*nl<<endl;
	size_t k = nk->rid;
	size_t l = nl->rid;
	
	if(!nk->is_ground()&& nk->isS()!=Y){
		 rhs[k] += net->Ieq;	// for VDD SubCircuit
		// clog<<*nk<<" rhs +: "<<rhs[k]<<endl;
	}
	if(!nl->is_ground()&& nl->isS()!=Y){
		 rhs[l] += -net->Ieq; 
		//clog<<*nl<<" rhs +: "<<rhs[l]<<endl;
	}
}

void SubCircuit::modify_rhs_Ieq_l(Net *net, double *rhs){
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
	// nk point to X node
	if(nk->isS() !=X){ 
		swap<Node*>(nk, nl);
		swap<Node*>(net->ab[0], net->ab[1]);
	}
	size_t k = nk->rid;
	size_t l = nl->rid;

	if(nk->isS() !=Y && !nk->is_ground()){
		 rhs[k] += net->Ieq; // VDD SubCircuit
		//clog<<*nk<<" "<<rhs[k]<<endl;
	}
	if(nl->isS()!=Y && !nl->is_ground()){
		 rhs[l] += -net->Ieq; // VDD SubCircuit
		//clog<<*nl<<" "<<rhs[l]<<endl;
	}
}

// reset b
void SubCircuit::reset_b(){
   size_t n = replist.size();
   b = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
   bp = static_cast<double *> (b->x);
}

void SubCircuit::stamp_rhs_tr(bool local_flag, double time, Tran &tran){
	size_t n = replist.size();
	b = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
   	bp = static_cast<double *> (b->x);

	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		switch(type){
		case CURRENT:
			if(local_flag == false){
				// stamp global current
				for(size_t i=0;i<ns.size();i++)
					stamp_current(bp, ns[i]);
			}
			else{// stamp local current
				for(size_t i=0;i<ns.size();i++)
					stamp_current_tr_net(bp, ns[i], time);
			}
			break;
		case VOLTAGE:
			for(size_t i=0;i<ns.size();i++){
				if( fzero(ns[i]->value)  && 
				    !ns[i]->ab[0]->is_ground() &&
				    !ns[i]->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_rhs_VDD(bp, ns[i]);
			}
			break;
		case RESISTOR:
		case CAPACITANCE:
		case INDUCTANCE:
			break;
		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}
	// make_A_symmetric
	if(local_flag == false)
		make_A_symmetric_tr(bp, tran);
	else
		make_A_symmetric_local(bp);
}

// restamp bp with LDO
void SubCircuit::stamp_rhs_DC(bool local_flag){
	size_t n = replist.size();
   	b = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
   	bp = static_cast<double *> (b->x);

	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		switch(type){
		case CURRENT:
			for(size_t i=0;i<ns.size();i++)
				stamp_current(bp, ns[i]);
			break;
		case VOLTAGE:
			for(size_t i=0;i<ns.size();i++){
				if( fzero(ns[i]->value)  && 
				    !ns[i]->ab[0]->is_ground() &&
				    !ns[i]->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_rhs_VDD(bp, ns[i]);
			}
			break;
		case RESISTOR:
		case CAPACITANCE:
			break;
		case INDUCTANCE:
			for(size_t i=0;i<ns.size();i++){
				stamp_induc_rhs_dc(bp, ns[i]);	
			}
			break;
		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}
	// only make A symmetric for local
	if(local_flag == true)
		make_A_symmetric_local(bp);	
	else
		make_A_symmetric(bp);
}

// stamp ldo VDD net into bp
void SubCircuit::stamp_rhs_VDD(double *bp, Net *net){
	Node * X = net->ab[0];
	if( X->is_ground() ) X = net->ab[1];
	size_t id = X->rep->rid;
	Net * south = X->rep->nbr[SOUTH];
	if( south != NULL &&
	    south->type == CURRENT ){
		bp[id] = net->value;	    // modify it
		// cout<<"b: ="<<id<<" "<<net->value<<endl;
	}
	else{
		bp[id] += net->value;
		// cout<<"b: +"<<id<<" "<<net->value<<endl;
	}
}

void SubCircuit::release_resources(){
	cholmod_free_dense(&b, cm);
	cholmod_free_factor(&L, cm);
	cholmod_free_dense(&x, cm);
	cholmod_finish(&c);
}

// ----------------------------------------------------------------//
// Filename : circuit.cpp
// Author : Zigang Xiao <zxiao2@illinois.edu>
//          Ting Xu <tingyu1@illinois.edu>
//
// implementation file of circuit.h
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
#include <algorithm>
#include <iostream>
#include <ctime>
#include <utility>
#include <cassert>
#include <vector>
#include "circuit.h"
#include "util.h"
#include "algebra.h"
#include "node.h"
using namespace std;

//////////////////////////////////////////////////////////////////////////
// Constructor and utility functions goes here

vector<LAYER_DIR> Circuit::layer_dir(MAX_LAYER);

// constructor of Circuit class, name is optional
Circuit::Circuit(string _name):name(_name),
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
Circuit::~Circuit(){
	pad_set.clear();
	special_nodes.clear();
	worst_cur.clear();
	map_node_pt.clear();
	for(size_t i=0;i<nodelist.size();i++) 
		delete nodelist[i];
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		for(size_t j=0;j<ns.size();j++) delete ns[j];
	}
	//delete [] id_map;
	Lx = NULL;
	Li = NULL;
	Lp = NULL;
	Lnz = NULL;
}

void Circuit::check_sys() const{
	clog<<"**** CHECKING SXSTEM ENVIRONMENT ****"<<endl;
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

bool compare_values(double a, double b){
	return (a<b);
}
// sort the nodes according to their coordinate 
void Circuit::sort_nodes(){
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

string Circuit::get_name() const{return this->name;}

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

ostream & operator << (ostream & os, const Circuit & ckt){
	os<<"Circuit ["<<ckt.name<<"] info:"<<endl;

	os<<"==== Nodes ===="<<endl;
	os<<ckt.nodelist;

	os<<"==== Reps  ===="<<endl;
	os<<ckt.replist;

	return os;
}

void Circuit::print(){
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

// initialization before solving the circuit
// 1. sort the nodes
// 2. set node representatives
// 3. find node in which block, update count
// 4. get representative lists
void Circuit::solve_init(){
        sort_nodes();
	size_t size = nodelist.size()-1; // exclude ground node!!
	Node *p=NULL;
	for(size_t i=0,nr=0;i<size;i++){
		p=nodelist[i];
		p->id = i;
		// set the representative
		Net * net = p->nbr[TOP];
		//if(net!=NULL)
			//cout<<"p, net: "<<*p<<" "<<p->isS()<<" "<<*net<<endl;
		if( p->isS()== X) VDD = p->get_value();
		// test short circuit
		if( p->isS() !=X && // X must be representative 
		    net != NULL &&
		    fzero(net->value) ){
			// TODO: ensure ab[1] is not p itself
			assert( net->ab[1] != p );
			p->rep = net->ab[1]->rep;
		} // else the representative is itself

		//cout<<"p, p->rep: "<<*p<<" "<<*p->rep<<endl;
		// push the representatives into list
		if( p->rep == p ) {
			replist.push_back(p);
			//rep_id[p] = nr; // save the id
			p->rid = nr;
			++nr;
		}
		else
			p->rid = p->rep->rid;
		//cout<<"p->rep, rid: "<<*p->rep<<" "<<p->rid<<endl;
	}
	//cout<<"nodelist.size: "<<nodelist<<endl;
	//clog<<"replist.size: "<<replist<<endl;
	build_pad_set();
	mark_special_nodes();
}

// count number of nodes that can be merged
void Circuit::count_merge_nodes(){
	size_t size = replist.size();
	size_t count = 0;
	for(size_t i=0;i<size;i++){
		Node * p = replist[i];
		if(p->nbr[TOP] ==NULL && p->nbr[BOTTOM] ==NULL){
			if((p->nbr[EAST] !=NULL && p->nbr[WEST] !=NULL) 
			&& (p->nbr[NORTH] ==NULL && p->nbr[SOUTH]==NULL))
				count++;
			else if((p->nbr[NORTH] !=NULL && p->nbr[SOUTH] !=NULL)&&(p->nbr[EAST]==NULL && p->nbr[WEST] ==NULL))
				count++;
		}
	}
	//clog<<"number of nodes can be merged is: "<<count<<endl;
}

void Circuit::solve(Tran &tran){
	solve_LU(tran);
	// cout<<nodelist<<endl;
	double max_IRdrop = locate_maxIRdrop_tr(tran);
			
	clog<<"max IRdrop for tr is: "<<max_IRdrop<<endl;	
	clog<<endl;

}

void Circuit::solve_DC(){
   replist.clear();
   for(size_t i=0;i<nodelist.size()-1;i++){
	nodelist[i]->rep = nodelist[i];
   }
   solve_init();
   size_t n = replist.size();	// replist dosn't contain ground node
   if( n == 0 ) return;		// No node    
   
   /*cout<<"========nodelist: ====="<<endl;
   for(size_t i=0;i<nodelist.size()-1;i++){
	cout<<"id, isS, node, rep: "<<nodelist[i]->rid<<" "<<nodelist[i]->rep->isS()<<" "<<*nodelist[i]<<" "<<*nodelist[i]->rep<<endl;
   }*/
   //cout<<nodelist<<endl;
   cm = &c;
   cholmod_start(cm);
   cm->print = 5;

   b = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
   x = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
   bp = static_cast<double *> (b->x);

   Matrix A;
   stamp_by_set(A, bp);
   make_A_symmetric(bp);
   A.set_row(n);
   //cout<<endl<<A<<endl;
   //cout<<"replist size: "<<n<<endl;
   //for(size_t i=0;i<n;i++)
	//cout<<"i, bp: "<<i<<" "<<bp[i]<<endl;
   Algebra::solve_CK(A, L, x, b, cm);
   A.clear();
   //return;
   xp = static_cast<double *> (x->x);
   for(size_t i=0;i<n;i++)
	replist[i]->value = xp[i];

   for(size_t i=0;i<nodelist.size();i++){
	nodelist[i]->value = nodelist[i]->rep->value;
   }
   cout<<nodelist<<endl;

   double max_IRdrop = locate_maxIRdrop();
			
   clog<<"max IRdrop is: "<<max_IRdrop<<endl;	
   double special_IRdrop = locate_special_maxIRdrop();
   clog<<"special IRdrop is: "<<special_IRdrop<<endl;
   clog<<endl;
}

// stamp the matrix and solve
void Circuit::solve_LU_core(Tran &tran){
   size_t n = replist.size();	// replist dosn't contain ground node
   if( n == 0 ) return;		// No node    
   cm = &c;
   cholmod_start(cm);
   cm->print = 5;

   b = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
   x = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
   bp = static_cast<double *> (b->x);

   Matrix A;
   stamp_by_set(A, bp);
   make_A_symmetric(bp);
   A.set_row(n);
   Algebra::solve_CK(A, L, x, b, cm);
   //return;
   xp = static_cast<double *> (x->x);
   for(size_t i=0;i<n;i++)
	replist[i]->value = xp[i];

   for(size_t i=0;i<nodelist.size();i++){
	nodelist[i]->value = nodelist[i]->rep->value;
   }
   //cout<<nodelist<<endl;

   // A is already being cleared   
   size_t i=0;
   for(i=0;i<replist.size();i++)
	bp[i] = 0;
   link_ckt_nodes(tran);

   bnew = cholmod_zeros(n,1,CHOLMOD_REAL, cm);
   bnewp = static_cast<double *>(bnew->x);
 
   double time = 0;
   //int iter = 0;
   stamp_by_set_tr(A, bp, tran);
   make_A_symmetric(bp);
   stamp_current_tr(bp, time);
 
   Algebra::CK_decomp(A, L, cm);
   Lp = static_cast<int *>(L->p);
   Lx = static_cast<double*> (L->x);
   Li = static_cast<int*>(L->i) ;
   Lnz = static_cast<int *>(L->nz); 
   A.clear();
   for(size_t i=0;i<n;i++){
	bnewp[i] = bp[i];
   }
  
   // set_eq_induc(tran);
   set_eq_capac(tran);
   modify_rhs_tr_0(bnewp, xp);
   map_node.clear();
   x = cholmod_solve(CHOLMOD_A, L, bnew, cm);
   xp = static_cast<double *> (x->x);
   //for(size_t i=0;i<n;i++)
	//replist[i]->value = xp[i];
    
   save_ckt_nodes(xp);
   time += tran.step_t;
   // maintain the old current values
   worst_cur_new = worst_cur;
   double IRdrop_old = locate_maxIRdrop(xp, n);
   double IRdrop_new = 0;
   int iter=1;
   // then start other iterations
   while(time < tran.tot_t){// && iter < 0){
	for(size_t i=0;i<n;i++)
		bnewp[i] = bp[i];
      // only stamps if net current changes
      // set bp into last state
      stamp_current_tr_1(bp, bnewp, time);
     // get the new bnewp
      modify_rhs_tr(bnewp, xp); 
	
      x = cholmod_solve(CHOLMOD_A, L, bnew, cm);
      xp = static_cast<double *> (x->x); 
      // then locate max_IRdrop
      IRdrop_new = locate_maxIRdrop(xp, n);
      iter++;
      // update the saved worst_cur
      if(IRdrop_new > IRdrop_old){
	worst_cur = worst_cur_new;
	IRdrop_old = IRdrop_new;		
      }

      save_ckt_nodes(xp);
      time += tran.step_t;
   }
   
   save_ckt_nodes_to_tr(tran);
   print_ckt_nodes(tran);
   // release_resource();
   release_ckt_nodes();
   cholmod_free_dense(&b, cm);
   cholmod_free_dense(&bnew, cm);
   cholmod_free_factor(&L, cm);
   cholmod_free_dense(&x, cm);
   cholmod_finish(&c);   
}

// solve the node voltages using direct LU
void Circuit::solve_LU(Tran &tran){
        solve_init();
	// initialize worst_cur vector
	worst_cur.resize(replist.size());
	for(size_t i=0;i<replist.size();i++){
		worst_cur[i] = 0;	
	}
	solve_LU_core(tran);
	//solve_LU_core_all(tran);
}

// given vector x that obtained from LU, set the value to the corresponding
// node in nodelist
void Circuit::get_voltages_from_LU_sol(double * x){
   size_t i;
   for(i=0;i<nodelist.size()-1;i++){
      Node * node = nodelist[i];
      size_t id = node->rep->rid;  // get rep's id in Vec
      double v = x[id];		// get its rep's value
      node->value = v;

   }
}

// copy node voltages from the circuit to a Vec
// from = true then copy circuit to x
// else copy from x to circuit
void Circuit::copy_node_voltages(double * x, size_t &size, bool from){
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
void Circuit::stamp_by_set(Matrix & A, double * b){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		switch(type){
		case RESISTOR:
			for(size_t i=0;i<ns.size();i++){
				if(fzero(ns[i]->value))
					continue;
				//cout<<"net: "<<*ns[i]<<endl;
				//assert( fzero(ns[i]->value) == false );
				stamp_resistor(A, ns[i]);
			}
			break;
		case CURRENT:
			for(size_t i=0;i<ns.size();i++){
				//cout<<"net: "<<*ns[i]<<endl;
				stamp_current(b, ns[i]);
			}
			break;
		case VOLTAGE:
			for(size_t i=0;i<ns.size();i++){
				if( fzero(ns[i]->value)  && 
				    !ns[i]->ab[0]->is_ground() &&
				    !ns[i]->ab[1]->is_ground() )
					continue; // it's a 0v via

				//cout<<"net: "<<*ns[i]<<endl;
				stamp_VDD(A, b, ns[i]);
			}
			break;
		case CAPACITANCE:
			//for(size_t i=0;i<ns.size();i++)
				//stamp_capacitance_dc(A, ns[i]);
			break;
		case INDUCTANCE:
			for(size_t i=0;i<ns.size();i++){
				stamp_inductance_dc(A, b, ns[i]);	
			}
			break;
		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}
}

// stamp transient current values into rhs
void Circuit::stamp_current_tr(double *b, double &time){
	NetPtrVector & ns = net_set[CURRENT];
	for(size_t i=0;i<ns.size();i++)
		stamp_current_tr_net(b, ns[i], time);
}

// stamp transient current values into rhs
void Circuit::stamp_current_tr_1(double *bp, double *b, double &time){
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
void Circuit::stamp_by_set_tr(Matrix & A, double *b, Tran &tran){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		switch(type){
		case RESISTOR:
			/*for(size_t i=0;i<ns.size();i++){
				if(fzero(ns[i]->value))
					continue;
				stamp_resistor_tr(A, ns[i]);
			}*/
			break;
		case CURRENT:
			//for(size_t i=0;i<ns.size();i++)
				//stamp_current_tr(b, ns[i]);
			break;
		case VOLTAGE:
			for(size_t i=0;i<ns.size();i++){
				if( fzero(ns[i]->value)  && 
				    !ns[i]->ab[0]->is_ground() &&
				    !ns[i]->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_VDD_tr(b, ns[i]);
			}
			break;
		case CAPACITANCE:
			for(size_t i=0;i<ns.size();i++)
				stamp_capacitance_tr(A, ns[i], tran);
			break;
		case INDUCTANCE:
			//for(size_t i=0;i<ns.size();i++){
				//stamp_inductance_tr(A, ns[i], tran);	
			//}
			break;
		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}
}

// update rhs by transient nets
void Circuit::modify_rhs_tr_0(double * b, double *x){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		if(type ==CAPACITANCE){	
			for(size_t i=0;i<ns.size();i++){
				modify_rhs_c_tr_0(ns[i], b, x);
			}
		}	
	}
}
// update rhs by transient nets
void Circuit::modify_rhs_tr(double * b, double *x){
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

void Circuit::stamp_resistor(Matrix & A, Net * net){
	//clog<<"net: "<<*net<<endl;
	double G;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	G = 1./net->value;
	//cout<<"nk, nl: "<<nk->isS()<<" "<<nl->isS()<<endl;
        if( !nk->is_ground()&& nk->isS()!=X) {
           A.push_back(k,k, G);

           //cout<<"("<<k<<" "<<k<<" "<<G<<")"<<endl;
           if(!nl->is_ground() &&nl->isS()!=X &&(k > l)){
                    A.push_back(k,l,-G);

                  //cout<<"("<<k<<" "<<l<<" "<<-G<<")"<<endl;
           }
        }

	if( !nl->is_ground() && nl->isS() !=X) {
		//cout<<"nb: "<<*nl<<" "<<nl->isS()<<" "<<X<<endl;
		A.push_back(l,l, G);
                //cout<<"("<<l<<" "<<l<<" "<<G<<")"<<endl;
		if(!nk->is_ground()&& nk->isS()!=X && l > k){
			A.push_back(l,k,-G);
		//cout<<"("<<l<<" "<<k<<" "<<-G<<")"<<endl;
                }
	}
}

// only stamp the resistor node connected to inductance
void Circuit::stamp_resistor_tr(Matrix & A, Net * net){
   double G;
   Node * nk = net->ab[0]->rep;
   Node * nl = net->ab[1]->rep;
   size_t k = nk->rid;
   size_t l = nl->rid;
   G = 1./net->value;

   if( nk->isS()!=X && !nk->is_ground()) {
      //cout<<"net: "<<*net<<endl;
      A.push_back(k,k, G);
      //cout<<"("<<k<<" "<<k<<" "<<G<<")"<<endl;
      if(!nl->is_ground() && nl->isS()!=X){
         if(l < k){
            A.push_back(k,l,-G);
            //cout<<"("<<k<<" "<<l<<" "<<-G<<")"<<endl;
         }
         else if(l > k){ 
            A.push_back(l, k, -G);
            //cout<<"("<<l<<" "<<k<<" "<<-G<<")"<<endl;
         }
      }
   }

   if( nl->isS() !=X && !nl->is_ground()) {

      //cout<<"net: "<<*net<<endl;
      A.push_back(l,l, G);
      //cout<<"("<<l<<" "<<l<<" "<<G<<")"<<endl;
      if(!nk->is_ground()&& nk->isS()!=X){
         if(k < l){
            A.push_back(l,k,-G);
            //cout<<"("<<l<<" "<<k<<" "<<-G<<")"<<endl;
         }
         else if(k > l){
            A.push_back(k, l, -G);
            //cout<<"("<<k<<" "<<l<<" "<<-G<<")"<<endl;
         }
      }
   }
}

void Circuit::stamp_inductance_dc(Matrix & A, double *b, Net * net){
	double G;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	G = 1./net->value;
	if( nk->isS()!=X && !nk->is_ground()) {
		A.push_back(k,k, 1);
		// general stamping
		if(!nl->is_ground())
		// A.push_back(k,l,-1);
		// make it symmetric
			b[k] = b[l];
		//clog<<"("<<k<<" "<<k<<" "<<1<<")"<<endl;
		//clog<<"("<<k<<" "<<l<<" "<<-1<<")"<<endl;
	}

	if( nl->isS() !=X && !nl->is_ground()) {
		A.push_back(l,l, 1);
		if(!nk->is_ground())
		// general stamping
		// A.push_back(l,k,-1);
		b[l] = b[k];
		//clog<<"("<<l<<" "<<l<<" "<<1<<")"<<endl;
		//clog<<"("<<l<<" "<<k<<" "<<-1<<")"<<endl;
	}
}

// stamp dc cap, infinitive resistance
void Circuit::stamp_capacitance_dc(Matrix & A, Net * net){
	double G;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	G = 1./net->value;
	if( nk->isS()!=X && !nk->is_ground()) {
		A.push_back(k,k, 0);
		if(!nl->is_ground())
			A.push_back(k,l,0);
	}

	if( nl->isS() !=X && !nl->is_ground()) {
		A.push_back(l,l, 0);
		if(!nk->is_ground())
			A.push_back(l,k,0);
	}
}

// stamp inductance Geq = delta_t/(2L)
void Circuit::stamp_inductance_tr(Matrix & A, Net * net, Tran &tran){
	//clog<<"net: "<<*net<<endl;
	double Geq = 0;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	// Geq = delta_t / (2*L)
	Geq = tran.step_t / (2*net->value);
	//net->value = Geq;

	if( nk->isS()!=X  && !nk->is_ground()) {
		// -1 is to clear formal inserted 1 at (k,k)
		A.push_back(k,k, Geq-1);
		//clog<<"("<<k<<" "<<k<<" "<<Geq-1<<")"<<endl;
		//clog<<nl->isS()<<endl;
		if(!nl->is_ground()&& nl->isS()!=X && k>l){
			A.push_back(k,l,-Geq);
		        //clog<<"("<<k<<" "<<l<<" "<<-Geq<<")"<<endl;
		}
	}

	if( nl->isS() !=X && !nl->is_ground()) {
		// -1 is to clear formal inserted 1 at (l,l)
		A.push_back(l,l, Geq-1);
		//clog<<"("<<l<<" "<<l<<" "<<Geq-1<<")"<<endl;
		if(!nk->is_ground() && nk->isS()!=X && l>k){
			A.push_back(l,k,-Geq);
			//clog<<"("<<l<<" "<<k<<" "<<-Geq<<")"<<endl;
		}
	}
}

// stamp capacitance Geq = 2C/delta_t
void Circuit::stamp_capacitance_tr(Matrix &A, Net *net, Tran &tran){
	//cout<<"net: "<<*net<<endl;
	double Geq = 0;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	// Geq = 2*C / delta_t
	Geq = (2*net->value) / tran.step_t;
	//net->value = Geq;
	//cout<<"C delta_t Geq: "<<net->value<<" "<<tran.step_t<<" "<<Geq<<endl;
	// Ieq = i(t) + 2*C / delta_t * v(t)

	if( nk->isS()!=X  && !nk->is_ground()) {
		A.push_back(k,k, Geq);
		//cout<<"("<<k<<" "<<k<<" "<<Geq<<")"<<endl;
		if(!nl->is_ground()&& k > l){
			A.push_back(k,l,-Geq);
			//cout<<"("<<k<<" "<<l<<" "<<-Geq<<")"<<endl;
		}
	}

	if( nl->isS() !=X && !nl->is_ground()) {
		A.push_back(l,l, Geq);
		//cout<<"("<<l<<" "<<l<<" "<<Geq<<")"<<endl;
		if(!nk->is_ground()&& l > k){
			A.push_back(l,k,-Geq);
			//cout<<"("<<l<<" "<<k<<" "<<-Geq<<")"<<endl;
		}
	}
}

// add Ieq into rhs
// Ieq = i(t) + 2*C / delta_t *v(t)
void Circuit::modify_rhs_c_tr_0(Net *net, double * rhs, double *x){
	double i_t = 0;
	double temp = 0;
	double Ieq = 0;
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;

        // nk point to Z node
        if(nk->isS()!=Z){
		swap<Node *>(nk, nl);
		swap<Node*>(net->ab[0], net->ab[1]);
	}
	size_t k = nk->rid;
	size_t l = nl->rid;

	Net *r = nk->nbr[TOP];
	Node *a = r->ab[0]->rep;
	Node *b = r->ab[1]->rep;
	// a point to Z node
	if(a->isS()!=Z) {
		swap<Node *>(a, b);
		swap<Node*>(r->ab[0], r->ab[1]);
	}

	size_t id_a = a->rid;
	size_t id_b = b->rid;
	//i_t = (b->value - a->value) / r->value;
	i_t = (x[id_b] - x[id_a]) / r->value;
	//#if 0
	
	if(nk->is_ground())
	 temp = net->value *(-x[l]);
        else if(nl->is_ground()){
	 temp = net->value *x[k];
        }
        else
	 temp = net->value *(x[k]-x[l]);
	Ieq  = (i_t + temp);
	if(!nk->is_ground()&& nk->isS()!=X){
		 rhs[k] += Ieq;	// for VDD circuit
		 worst_cur[k] += Ieq;
	}
	if(!nl->is_ground()&& nl->isS()!=X){
		 rhs[l] += -Ieq; 
		 worst_cur[l] += -Ieq;
	}
}

// add Ieq into rhs
// Ieq = i(t) + 2*C / delta_t *v(t)
void Circuit::modify_rhs_c_tr(Net *net, double * rhs, double *x){
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
	if(!nk->is_ground()&& nk->isS()!=X){
		 rhs[k] += Ieq;	// for VDD circuit
		 worst_cur_new[k] += Ieq;
	}
	if(!nl->is_ground()&& nl->isS()!=X){
		 rhs[l] -= Ieq; 
		 worst_cur_new[l] += -Ieq;
	}
}

void Circuit::set_eq_induc(Tran &tran){
	NetPtrVector &ns = net_set[INDUCTANCE];
	for(size_t i=0;i<ns.size();i++)
		ns[i]->value = tran.step_t /(2*ns[i]->value);
}

void Circuit::set_eq_capac(Tran &tran){
	NetPtrVector &ns = net_set[CAPACITANCE];
	for(size_t i=0;i<ns.size();i++)
		ns[i]->value = 2*ns[i]->value/tran.step_t;
}

// add Ieq into rhs
// Ieq = i(t) + delta_t / (2*L) *v(t)
void Circuit::modify_rhs_l_tr_0(Net *net, double *rhs, double *x){
	//clog<<"l net: "<<*net<<endl;
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
	// nk point to X node
	if(nk->isS() !=X){ 
		swap<Node*>(nk, nl);
		swap<Node*>(net->ab[0], net->ab[1]);
	}
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
	if(a->isS()!=X) {
		swap<Node*>(a, b);
		swap<Node*>(r->ab[0], r->ab[1]);
	}
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
//#if 0
//#endif
	Ieq  = i_t + temp;
	//clog<<"Ieq: "<<Ieq<<endl;
	if(nk->isS() !=X && !nk->is_ground()){
		 rhs[k] += Ieq; // VDD circuit
		//clog<<*nk<<" "<<rhs[k]<<endl;
	}
	if(nl->isS()!=X && !nl->is_ground()){
		 rhs[l] += -Ieq; // VDD circuit
		//clog<<*nl<<" "<<rhs[l]<<endl;
	}
}

// add Ieq into rhs
// Ieq = i(t) + delta_t / (2*L) *v(t)
void Circuit::modify_rhs_l_tr(Net *net, double *rhs, double *x){
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
	if(nk->isS() !=X && !nk->is_ground()){
		 rhs[k] += Ieq; // VDD circuit
		//clog<<*nk<<" "<<rhs[k]<<endl;
	}
	if(nl->isS()!=X && !nl->is_ground()){
		 rhs[l] += -Ieq; // VDD circuit
		//clog<<*nl<<" "<<rhs[l]<<endl;
	}
}
// stamp a current source
void Circuit::stamp_current(double * b, Net * net){
	//cout<<"net: "<<*net<<endl;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;

	if( !nk->is_ground() && nk->isS()!=X){// && 
		size_t k = nk->rid;
		b[k] += -net->value;
		//cout<<"b: "<<k<<" "<<-net->value<<endl;
	}
	if( !nl->is_ground() && nl->isS() !=X){// &&
		size_t l = nl->rid;
		b[l] +=  net->value;
		//cout<<"b: "<<l<<" "<<-net->value<<endl;
	}
}

void Circuit::stamp_current_tr_net(double * b, Net * net, double &time){
	current_tr(net, time);
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	if( !nk->is_ground()&& nk->isS()!=X) { 
		size_t k = nk->rid;
		//clog<<"node, rid: "<<*nk<<" "<<k<<endl;
		b[k] += -net->value;//current;
		worst_cur[k] = -net->value;
		//clog<<"time, k, b: "<<time<<" "<<k<<" "<<b[k]<<endl;
	}
	if( !nl->is_ground() && nl->isS()!=X) {
		size_t l = nl->rid;
		//clog<<"node, rid: "<<*nl<<" "<<l<<endl;
		b[l] +=  net->value;// current;
		worst_cur[l] = net->value;
		//clog<<"time, l, b: "<<time<<" "<<l<<" "<<b[l]<<endl;
	}
}

void Circuit::stamp_current_tr_net_1(double *bp, double * b, Net * net, double &time){
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
		if( !nk->is_ground()&& nk->isS()!=X) { 
			size_t k = nk->rid;
			//clog<<"node, rid: "<<*nk<<" "<<k<<endl;
			//clog<<"time, k, b bef: "<<time<<" "<<k<<" "<<b[k]<<endl;
			b[k] += -diff;//current;
			bp[k] = b[k];
			worst_cur_new[k] = -diff;
			//clog<<"time, k, b: "<<time <<" "<<k<<" "<<b[k]<<endl;
		}
		if( !nl->is_ground() && nl->isS()!=X) {
			size_t l = nl->rid;
			//clog<<"time, l, b bef: "<<time<<" "<<l<<" "<<b[l]<<endl;
			//clog<<"node, rid: "<<*nl<<" "<<l<<endl;
			b[l] +=  diff;// current;
			bp[l] = b[l];
			worst_cur_new[l] = diff;
			//clog<<"time, l, b: "<<time<<" "<<l<<" "<<b[l]<<endl;
		}
	}
}

// stamp a voltage source
void Circuit::stamp_VDD(Matrix & A, double * b, Net * net){
	// find the non-ground node
	//clog<<"net: "<<*net<<endl;
	Node * X = net->ab[0];
	if( X->is_ground() ) X = net->ab[1];
	size_t id = X->rep->rid;
	A.push_back(id, id, 1.0);
	//clog<<"push id, id, 1: "<<id<<" "<<id<<" "<<1<<endl;
	Net * south = X->rep->nbr[SOUTH];
	if( south != NULL &&
	    south->type == CURRENT ){
		// this node connects to a VDD and a current
		assert( feqn(1.0, b[id]) ); // the current should be stamped
		b[id] = net->value;	    // modify it
		//clog<<"b: ="<<id<<" "<<net->value<<endl;
	}
	else{
		b[id] += net->value;
		//clog<<"b: +"<<id<<" "<<net->value<<endl;
	}
}

// stamp a voltage source
void Circuit::stamp_VDD_tr(double * b, Net * net){
	// find the non-ground node
	//clog<<"net: "<<*net<<endl;
	Node * X = net->ab[0];
	if( X->is_ground() ) X = net->ab[1];
	size_t id = X->rep->rid;
	//A.push_back(id, id, 1.0);
	Net * south = X->rep->nbr[SOUTH];
	if( south != NULL &&
	    south->type == CURRENT ){
		// this node connects to a VDD and a current
		assert( feqn(1.0, b[id]) ); // the current should be stamped
		b[id] = net->value;	    // modify it
		//clog<<"b: ="<<id<<" "<<net->value<<endl;
	}
	else{
		b[id] += net->value;
		//clog<<"b: +"<<id<<" "<<net->value<<endl;
	}
}

// decide transient step current values
void Circuit::current_tr(Net *net, double &time){
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
void Circuit:: save_tr_nodes(Tran &tran, double *x){
   size_t id=0;
   for(size_t j=0;j<tran.nodes.size();j++){
      if(tran.nodes[j].node != NULL ){
         id = tran.nodes[j].node->rep->rid;
         tran.nodes[j].value.push_back(x[id]);
      }
   }
}

// assign value back to transient nodes
void Circuit:: save_ckt_nodes(double *x){
   size_t id=0;
   for(size_t j=0;j<ckt_nodes.size();j++){
	//if(ckt_nodes[j].node->name == "n0_4_75")
         id = ckt_nodes[j].node->rep->rid;
	
	//if(ckt_nodes[j].node->name == "n0_4_75")
         ckt_nodes[j].value.push_back(x[id]);
      }
}

// link transient nodes with nodelist
void Circuit:: link_ckt_nodes(Tran &tran){
   Node_TR_PRINT nodes_temp;
   
   for(size_t i=0;i<tran.nodes.size();i++){
	  tran.nodes[i].node = get_node(tran.nodes[i].name);
	  
	   if(tran.nodes[i].node == NULL){
		 continue;
	   }
	   nodes_temp.flag = i;
	   nodes_temp.node = tran.nodes[i].node;
	   ckt_nodes.push_back(nodes_temp);
   } 
}

// link transient nodes with nodelist
void Circuit:: link_tr_nodes(Tran &tran){
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

void Circuit:: release_tr_nodes(Tran &tran){
   for(size_t j=0;j<tran.nodes.size();j++){
      if(tran.nodes[j].flag != -1){
         tran.nodes[j].node = NULL;
      }
   }
}

void Circuit:: release_ckt_nodes(){
   for(size_t j=0;j<ckt_nodes.size();j++){
         ckt_nodes[j].node = NULL;
   }
}

void Circuit::make_A_symmetric(double *b){
	int type = RESISTOR;
	NetList & ns = net_set[type];
	NetList::iterator it;
	Node *p=NULL, *q=NULL;

	for(it=ns.begin();it!=ns.end();it++){
           if( (*it) == NULL ) continue;
	   if(fzero((*it)->value)) continue;
           //assert( fzero((*it)->value) == false );
           if(!((*it)->ab[0]->rep->isS()==X || (*it)->ab[1]->rep->isS()==X)) continue;
	   if((*it)->ab[0]->rep->isS()==X && (*it)->ab[1]->rep->isS()==X)  continue;
	   //cout<<"symmetric net: "<<*(*it)<<endl;
           // node p points to X node
           if((*it)->ab[0]->rep->isS()==X){
              p = (*it)->ab[0]->rep; q = (*it)->ab[1]->rep;
           } 
           else if((*it)->ab[1]->rep->isS()==X){
              p = (*it)->ab[1]->rep; q = (*it)->ab[0]->rep;
           }           
           size_t id = q->rid;
           double G = 1.0 / (*it)->value;
           
           b[id] += p->value * G;
        }
}

void Circuit::make_A_symmetric_tr(double *b, double *x, Tran &tran){
	int type = RESISTOR;
	NetList & ns = net_set[type];
	NetList::iterator it;
	Node *p = NULL;
	Node *q = NULL;

	for(it=ns.begin();it!=ns.end();it++){
           if( (*it) == NULL ) continue;
		clog<<"net: "<<*(*it)<<endl;
           if( fzero((*it)->value))
		continue;
           if(!((*it)->ab[0]->rep->isS()==X || (*it)->ab[1]->rep->isS()==X)) continue;
           clog<<"net: "<<*(*it)<<endl;
           // node p points to X node
           if((*it)->ab[0]->rep->isS()==X){
              p = (*it)->ab[0]->rep; q = (*it)->ab[1]->rep;
           } 
           else if((*it)->ab[1]->rep->isS()==X ){
              p = (*it)->ab[1]->rep; q = (*it)->ab[0]->rep;
           }           
           clog<<"p and q: "<<*p<<" "<<*q<<endl;

           size_t id = q->rid;
           double G = tran.step_t / ((*it)->value*2);
           
           //b[id] += p->value * G;
           b[id] += x[p->rid] *G;
           clog<<"stamp p->value, G, b: "<<p->value<<" "<<G<<" "<<b[id]<<endl;
        }
}

bool compare_Node_G(const Node_G *nd_1, const Node_G *nd_2){
   return (nd_1->value < nd_2->value);
 }

void Circuit::update_node_set_bx(){
    int id = 0;
    //clog<<"len_node_set_b. "<<pg.node_set_b.size();
    for(size_t i=0;i<pg.node_set_b.size();i++){
       id = pg.node_set_b[i];
       pg.node_set_b[i] = id_map[id];
       //clog<<"b_old, b_new: "<<id<<" "<<id_map[id]<<endl;
    }
    //clog<<endl;
    //clog<<"len_node_set_x. "<<pg.node_set_x.size()<<endl;
    for(size_t i=0;i<pg.node_set_x.size();i++){
       id = pg.node_set_x[i];
       pg.node_set_x[i] = id_map[id];
       //clog<<"x_old, x_new: "<<id<<" "<<id_map[id]<<endl;
    }
 }


void Circuit::parse_path_table(){
   // build up nodelist info
      Node_G *node;
      for(size_t i=0;i<replist.size();i++){
         node = new Node_G();
         node->value = i;
         pg.nodelist.push_back(node);
      }
}

void Circuit::build_path_graph(){
   build_FFS_path();

   build_FBS_path();

   // only keep the 2 paths, switch from List into array
   len_path_b = pg.path_FFS.get_size();
   len_path_x = pg.path_FBS.get_size();

   path_b = new int[len_path_b];
   path_x = new int [len_path_x];
   
   Node_G *nd;
   nd = pg.path_FFS.first;
   for(int i=0;i<len_path_b;i++){
      path_b[i] = nd->value;
      if(nd->next != NULL)
         nd = nd->next;
   }
   pg.path_FFS.destroy_list();

   nd = pg.path_FBS.first;
   for(int i=0;i<len_path_x;i++){
      path_x[i] = nd->value;
      if(nd->next != NULL)
         nd = nd->next;
   }
   pg.path_FBS.destroy_list();

   pg.nodelist.clear();
   pg.node_set_b.clear();
   pg.node_set_x.clear();
}

void Circuit::build_FFS_path(){
   parse_path_table(); 

   set_up_path_table();

   find_path(pg.node_set_b, pg.path_FFS);

   pg.path_FFS.assign_size();

   for(size_t i=0;i<replist.size();i++)
      pg.nodelist[i]->flag = 0;
}

void Circuit::build_FBS_path(){
  pg.nodelist.clear();
  parse_path_table();
  set_up_path_table();
  find_path(pg.node_set_x, pg.path_FBS);
   pg.path_FBS.assign_size();
}

void Circuit::set_up_path_table(){
   size_t n = L->n;
   //int *Lp, *Li, *Lnz;
   int p, lnz, s, e;
   //Lp = static_cast<int *> (L->p);
   //Lnz = (int *)L->nz;
   //Li = static_cast<int *> (L->i);
   for(size_t i=0;i<n;i++){
      p = Lp[i];
      lnz = Lnz[i];

      s = Li[p];
      e = s;
      if(lnz >1) 
         e = Li[p+1];

      if(s<e)
         pg.nodelist[s]->next = pg.nodelist[e];
   }
}

void Circuit::find_path(vector<size_t> &node_set, List_G &path){
   Node_G* ne = pg.nodelist[pg.nodelist.size()-1];
   vector <Node_G *> insert_list;
   sort(node_set.begin(), node_set.end());
   if(node_set.size() == 0) return;
   
   // build up the first path start with id = min 
   int id = node_set[0];
   do{
      path.add_node(pg.nodelist[id]);
      pg.nodelist[id]->flag =1;
      if(pg.nodelist[id]->next == NULL) break;
      pg.nodelist[id] = pg.nodelist[id]->next;
   }while(pg.nodelist[id]->value != ne->value);
   path.add_node(ne);

   for(size_t i=0; i<node_set.size();i++){
      int id = node_set[i];
      if(pg.nodelist[id]->flag == 1) continue;
      // stops at first place where flag = 0
      // clog<<"stops at i: "<<i<<endl;
      do{
         if(pg.nodelist[id]->flag ==0){
            insert_list.push_back(pg.nodelist[id]);
            pg.nodelist[id]->flag =1;
         }
         if(pg.nodelist[id]->next == NULL || 
           pg.nodelist[id]->next->flag ==1)
            break;
         // else clog<<"next node: "<<*pg.nodelist[id]->next;
         pg.nodelist[id] = pg.nodelist[id]->next;
      }while(pg.nodelist[id]->value != ne->value); 
   }

   //clog<<"insert_list.size: "<<insert_list.size()<<endl;
   sort(insert_list.begin(), insert_list.end(), compare_Node_G);
   //for(int i=0;i<insert_list.size();i++)
      //clog<<"i, insert: "<<i<<" "<<*insert_list[i]<<endl;

   //clog<<"path: "<<&path<<endl;
   // p is the old pointer to the list
   // will be updated into new one
   Node_G *q=NULL;
   Node_G *p=NULL;
//#if 0
   for(size_t k=0;k<insert_list.size();k++){
      if(k ==0) p = path.first;
      else p = q;
      q = path.insert_node(insert_list[k], p);
   }
//#endif

   //clog<<"path: "<<&path<<endl;
   //clog<<endl;
   insert_list.clear();
}

#if 0
void Circuit::solve_eq_set(){
   int p, q, r, lnz, pend;
   int j, n = L->n ;

   // FFS solve
   for (j = 0 ; j < n ; ){
      /* get the start, end, and length of column j */
      p = Lp [j] ;
      lnz = Lnz [j] ;
      pend = p + lnz ;

      if (lnz < 4 || lnz != Lnz [j+1] + 1 || Li [p+1] != j+1)
      {
	s_col_FFS.push_back(1); 
         j++ ;	/* advance to next column of L */

      }
      //#if 0
      else if (lnz != Lnz [j+2] + 2 || Li [p+2] != j+2)
      {
	s_col_FFS.push_back(2); 
         j += 2 ;	    /* advance to next column of L */

      }
      else
      {
	s_col_FFS.push_back(3); 
         j += 3 ;	    /* advance to next column of L */
      }
      //#endif
   }
   //t2 = clock();
   //clog<<"FFS cost: "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;

   //t1 = clock();
   // FBS solve
   for(j = n-1; j >= 0; ){

      /* get the start, end, and length of column j */
      p = Lp [j] ;
      lnz = Lnz [j] ;
      pend = p + lnz ;

      /* find a chain of supernodes (up to j, j-1, and j-2) */

      if (j < 4 || lnz != Lnz [j-1] - 1 || Li [Lp [j-1]+1] != j)
      {
	s_col_FBS.push_back(1); 
         j--;
      }
      else if (lnz != Lnz [j-2]-2 || Li [Lp [j-2]+2] != j)
      {
	s_col_FBS.push_back(2);
        j -= 2 ;	    /* advance to the next column of L */

      }
      else
      {
	s_col_FBS.push_back(3);
        j -= 3 ;	    /* advance to the next column of L */
      }
      //#endif
   }
}
#endif

void Circuit::solve_eq(double *X){
   int p, q, r, lnz, pend;
   int j, n = L->n ;
   // assign xp[i] = bnewp[i]
   //if(n<THRESHOLD){
	   for(int i=0;i<n;i++)
		   X[i] = bnewp[i];
   //}
#if 0
   else{
	int i=0;
#pragma omp parallel for private(i)
	for(i=0;i<n;i++)
		X[i] = bnewp[i];
   }
#endif
   size_t count = 0;
   // FFS solve
   for (j = 0 ; j < n ; ){
      /* get the start, end, and length of column j */
      p = Lp [j] ;
      lnz = Lnz [j] ;
      pend = p + lnz ;

      if (lnz < 4 || lnz != Lnz [j+1] + 1 || Li [p+1] != j+1)
      {

         /* -------------------------------------------------------------- */
         /* solve with a single column of L */
         /* -------------------------------------------------------------- */
	 //s_col_FFS[j] = 1;
         double y = X [j] ;
         if(L->is_ll == true){
            X[j] /= Lx [p] ;
         }

         for (p++ ; p < pend ; p++)
         {
            X [Li [p]] -= Lx [p] * y ;
         }
         j++ ;	/* advance to next column of L */

      }
      //#if 0
      else if (lnz != Lnz [j+2] + 2 || Li [p+2] != j+2)
      {

         /* -------------------------------------------------------------- */
         /* solve with a supernode of two columns of L */
         /* -------------------------------------------------------------- */
	 //s_col_FFS[j] = 2;
         double y [2] ;
         q = Lp [j+1] ;
         if(L->is_ll == true){
            y [0] = X [j] / Lx [p] ;
            y [1] = (X [j+1] - Lx [p+1] * y [0]) / Lx [q] ;
            X [j  ] = y [0] ;
            X [j+1] = y [1] ;
         }

         else{
            y [0] = X [j] ;
            y [1] = X [j+1] - Lx [p+1] * y [0] ;
            X [j+1] = y [1] ;
         }
         for (p += 2, q++ ; p < pend ; p++, q++)
         {
            X [Li [p]] -= Lx [p] * y [0] + Lx [q] * y [1] ;
         }
         j += 2 ;	    /* advance to next column of L */

      }
      else
      {

         /* -------------------------------------------------------------- */
         /* solve with a supernode of three columns of L */
         /* -------------------------------------------------------------- */
	 //s_col_FFS[j] = 3;
         double y [3] ;
         q = Lp [j+1] ;
         r = Lp [j+2] ;
         if(L->is_ll == true){
            y [0] = X [j] / Lx [p] ;
            y [1] = (X [j+1] - Lx [p+1] * y [0]) / Lx [q] ;
            y [2] = (X [j+2] - Lx [p+2] * y [0] - Lx [q+1] * y [1]) / Lx [r] ;
            X [j  ] = y [0] ;
            X [j+1] = y [1] ;
            X [j+2] = y [2] ;
         }

         else{
            y [0] = X [j] ;
            y [1] = X [j+1] - Lx [p+1] * y [0] ;
            y [2] = X [j+2] - Lx [p+2] * y [0] - Lx [q+1] * y [1] ;
            X [j+1] = y [1] ;
            X [j+2] = y [2] ;
         }
         for (p += 3, q += 2, r++ ; p < pend ; p++, q++, r++)
         {
            X [Li [p]] -= Lx [p] * y [0] + Lx [q] * y [1] + Lx [r] * y [2] ;
         }
         j += 3 ;	    /* advance to next column of L */
      }
      //#endif
   }
   //t2 = clock();
   //clog<<"FFS cost: "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;

    count = 0;
   //t1 = clock();
   // FBS solve
   for(j = n-1; j >= 0; ){

      /* get the start, end, and length of column j */
      p = Lp [j] ;
      lnz = Lnz [j] ;
      pend = p + lnz ;

      /* find a chain of supernodes (up to j, j-1, and j-2) */

      if (j < 4 || lnz != Lnz [j-1] - 1 || Li [Lp [j-1]+1] != j)
      {

         /* -------------------------------------------------------------- */
         /* solve with a single column of L */
         /* -------------------------------------------------------------- */
	 //s_col_FBS[j] = 1;
         double d = Lx [p] ;
         if(L->is_ll == false)
            X[j] /= d ;
         for (p++ ; p < pend ; p++)
         {
            X[j] -= Lx [p] * X [Li [p]] ;
         }
         if(L->is_ll == true)
            X [j] /=  d ;
         j--;
      }
      else if (lnz != Lnz [j-2]-2 || Li [Lp [j-2]+2] != j)
      {

         /* -------------------------------------------------------------- */
         /* solve with a supernode of two columns of L */
         /* -------------------------------------------------------------- */
	 //s_col_FBS[j] =2;
         double y [2], t ;
         q = Lp [j-1] ;
         double d [2] ;
         d [0] = Lx [p] ;
         d [1] = Lx [q] ;
         t = Lx [q+1] ;
         if(L->is_ll == false){
            y [0] = X [j  ] / d [0] ;
            y [1] = X [j-1] / d [1] ;
         }
         else{
            y [0] = X [j  ] ;
            y [1] = X [j-1] ;
         }
         for (p++, q += 2 ; p < pend ; p++, q++)
         {
            int i = Li [p] ;
            y [0] -= Lx [p] * X [i] ;
            y [1] -= Lx [q] * X [i] ;
         }
         if(L->is_ll == true){
            y [0] /= d [0] ;
            y [1] = (y [1] - t * y [0]) / d [1] ;
         }
         else
            y [1] -= t * y [0] ;
         X [j  ] = y [0] ;
         X [j-1] = y [1] ;
         j -= 2 ;	    /* advance to the next column of L */

      }
      else
      {

         /* -------------------------------------------------------------- */
         /* solve with a supernode of three columns of L */
         /* -------------------------------------------------------------- */
         double y [3], t [3] ;
         q = Lp [j-1] ;
         r = Lp [j-2] ;
         double d [3] ;
         d [0] = Lx [p] ;
         d [1] = Lx [q] ;
         d [2] = Lx [r] ;
         t [0] = Lx [q+1] ;
         t [1] = Lx [r+1] ;
         t [2] = Lx [r+2] ;
         if(L->is_ll == false){
            y [0] = X [j]   / d [0] ;
            y [1] = X [j-1] / d [1] ;
            y [2] = X [j-2] / d [2] ;
         }
         else{
            y [0] = X [j] ;
            y [1] = X [j-1] ;
            y [2] = X [j-2] ;
         }
         for (p++, q += 2, r += 3 ; p < pend ; p++, q++, r++)
         {
            int i = Li [p] ;
            y [0] -= Lx [p] * X [i] ;
            y [1] -= Lx [q] * X [i] ;
            y [2] -= Lx [r] * X [i] ;
         }
         if(L->is_ll == true){
            y [0] /= d [0] ;
            y [1] = (y [1] - t [0] * y [0]) / d [1] ;
            y [2] = (y [2] - t [2] * y [0] - t [1] * y [1]) / d [2] ;
         }
         else{
            y [1] -= t [0] * y [0] ;
            y [2] -= t [2] * y [0] + t [1] * y [1] ;
         }
         X [j-2] = y [2] ;
         X [j-1] = y [1] ;
         X [j  ] = y [0] ;
         j -= 3 ;	    /* advance to the next column of L */
      }
      //#endif
   }
}

 // find super node columns for path_b and path_x
void Circuit::find_super(){
    int p, lnz;
    int j, k;
    // FFS loop
    for(k=0;k<len_path_b;){
       j = path_b[k];
       p = Lp[j];
       lnz = Lnz[j];
       if (lnz < 4 || path_b[k+1]!=j+1 || lnz != Lnz [j+1] + 1
         || Li [p+1] != j+1){
          s_col_FFS[k]= 1;
          k++;
       }
       else if (path_b[k+2]!= j+2 || lnz != Lnz [j+2] + 2
         || Li [p+2] != j+2){
         s_col_FFS[k]=2;
         k+=2;
       }
       else{
          s_col_FFS[k]=3;
          k+=3;
       }
    }
    //FBS loop
    for(k=len_path_x-1;k>=0;){
       j = path_x[k];
       p = Lp[j];
       lnz = Lnz[j];
       if (j < 4 || path_x[k-1]!=j-1||lnz != Lnz [j-1] - 1
         || Li [Lp [j-1]+1] != j){
         s_col_FBS[k]=1;
         k--;
       }
       else if (path_x[k-2] != j-2 ||lnz != Lnz [j-2]-2 ||
         Li [Lp [j-2]+2] != j){
         s_col_FBS[k]=2;
         k-=2;
       }
       else{
          s_col_FBS[k]=3;
          k-=3;
       }
    }
 }
 
 void Circuit::solve_eq_sp(double *X){
    int p, q, r, lnz, pend;
    int j, k, n = L->n ;
    for(int i=0;i<n;i++){
       X[i] = bnewp[i];
    }
    // FFS solve
    for(k=0; k < len_path_b;){
       j = path_b[k];
       //for (j = 0 ; j < n ; ){
       /* get the start, end, and length of column j */
       p = Lp [j] ;
       lnz = Lnz [j] ;
       pend = p + lnz ;
 
       if (s_col_FFS[k]==1)
	//if(lnz < 4 || lnz != Lnz [j+1] + 1 || Li [p+1] != j+1)
       {
 
          /* -------------------------------------------------------------- */
          /* solve with a single column of L */
          /* -------------------------------------------------------------- */
          double y = X [j] ;
          if(L->is_ll == true){
             X[j] /= Lx [p] ;
          }
          for (p++ ; p < pend ; p++)
          {
             X [Li [p]] -= Lx [p] * y ;
          }
          k++ ;  /* advance to next column of L */
 
       }
       else if (s_col_FFS[k]==2)
	//else if(lnz != Lnz [j+2] + 2 || Li [p+2] != j+2)
       {
          {
             double y [2] ;
             q = Lp [j+1] ;
             if(L->is_ll == true){
                y [0] = X [j] / Lx [p] ;
                y [1] = (X [j+1] - Lx [p+1] * y [0]) / Lx [q] ;
                X [j  ] = y [0] ;
                X [j+1] = y [1] ;
             }
 
             else{
                y [0] = X [j] ;
                y [1] = X [j+1] - Lx [p+1] * y [0] ;
                X [j+1] = y [1] ;
             }
             for (p += 2, q++ ; p < pend ; p++, q++)
             {
                X [Li [p]] -= Lx [p] * y [0] + Lx [q] * y [1] ;
             }
 
          }
          k += 2 ;           /* advance to next column of L */
 
       }
       else
       {
          {
             double y [3] ;
             q = Lp [j+1] ;
             r = Lp [j+2] ;
             //#ifdef LL
             if(L->is_ll == true){
                y [0] = X [j] / Lx [p] ;
                y [1] = (X [j+1] - Lx [p+1] * y [0]) / Lx [q] ;
                y [2] = (X [j+2] - Lx [p+2] * y [0] - Lx [q+1] * y [1]) / Lx [r] ;
                X [j  ] = y [0] ;
                X [j+1] = y [1] ;
                X [j+2] = y [2] ;
             }
 
             else{
                y [0] = X [j] ;
                y [1] = X [j+1] - Lx [p+1] * y [0] ;
                y [2] = X [j+2] - Lx [p+2] * y [0] - Lx [q+1] * y [1] ;
                X [j+1] = y [1] ;
                X [j+2] = y [2] ;
             }
             for (p += 3, q += 2, r++ ; p < pend ; p++, q++, r++)
             {
                X [Li [p]] -= Lx [p] * y [0] + Lx [q] * y [1] + Lx [r] * y [2] ;
             }
          }
             // marched to next 2 columns of L
          k += 3;
       }
    }
    // FBS solve
    for(k = len_path_x - 1; k >=0;){
       j = path_x[k];
       //for(j = n-1; j >= 0; ){
 
       /* get the start, end, and length of column j */
       p = Lp [j] ;
       lnz = Lnz [j] ;
       pend = p + lnz ;
 
       /* find a chain of supernodes (up to j, j-1, and j-2) */
 
        if (s_col_FBS[k]==1)
	//if(j < 4 || lnz != Lnz [j-1] - 1 || Li [Lp [j-1]+1] != j)
       {
 
          /* -------------------------------------------------------------- */
          /* solve with a single column of L */
          /* -------------------------------------------------------------- */
 
          double d = Lx [p] ;
          if(L->is_ll == false)
             X[j] /= d ;
          for (p++ ; p < pend ; p++)
          {
             X[j] -= Lx [p] * X [Li [p]] ;
          }
          if(L->is_ll == true)
             X [j] /=  d ;
          k--;
       }
       else if (s_col_FBS[k]==2)
	//else if(lnz != Lnz [j-2]-2 || Li [Lp [j-2]+2] != j)
       {
          {
             double y [2], t ;
             q = Lp [j-1] ;
             double d [2] ;
             d [0] = Lx [p] ;
             d [1] = Lx [q] ;
             t = Lx [q+1] ;
             if(L->is_ll == false){
                y [0] = X [j  ] / d [0] ;
                y [1] = X [j-1] / d [1] ;
             }
             else{
                y [0] = X [j  ] ;
                y [1] = X [j-1] ;
             }
             for (p++, q += 2 ; p < pend ; p++, q++)
             {
                int i = Li [p] ;
                y [0] -= Lx [p] * X [i] ;
                y [1] -= Lx [q] * X [i] ;
             }
             if(L->is_ll == true){
                y [0] /= d [0] ;
                y [1] = (y [1] - t * y [0]) / d [1] ;
             }
             else
                y [1] -= t * y [0] ;
             X [j  ] = y [0] ;
             X [j-1] = y [1] ;
          }
          k -= 2;
       }
       else
       {
          {
             double y [3], t [3] ;
             q = Lp [j-1] ;
             r = Lp [j-2] ;
             double d [3] ;
             d [0] = Lx [p] ;
             d [1] = Lx [q] ;
             d [2] = Lx [r] ;
             t [0] = Lx [q+1] ;
             t [1] = Lx [r+1] ;
             t [2] = Lx [r+2] ;
             if(L->is_ll == false){
                y [0] = X [j]   / d [0] ;
                y [1] = X [j-1] / d [1] ;
                y [2] = X [j-2] / d [2] ;
             }
             else{
                y [0] = X [j] ;
                y [1] = X [j-1] ;
                y [2] = X [j-2] ;
             }
             for (p++, q += 2, r += 3 ; p < pend ; p++, q++, r++)
             {
                int i = Li [p] ;
                y [0] -= Lx [p] * X [i] ;
                y [1] -= Lx [q] * X [i] ;
                y [2] -= Lx [r] * X [i] ;
             }
             if(L->is_ll == true){
                y [0] /= d [0] ;
                y [1] = (y [1] - t [0] * y [0]) / d [1] ;
                y [2] = (y [2] - t [2] * y [0] - t [1] * y [1]) / d [2] ;
             }
             else{
                y [1] -= t [0] * y [0] ;
                y [2] -= t [2] * y [0] + t [1] * y [1] ;
             }
             X [j-2] = y [2] ;
             X [j-1] = y [1] ;
             X [j  ] = y [0] ;
          }
          k -= 3;
       }
    }
 }

void Circuit::print_ckt_nodes(Tran &tran){
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

// calculate maximum IR drop from the observation nodelist
double Circuit::locate_maxIRdrop_tr(Tran &tran){
	double time = 0;
	size_t j=0;
	double max_IRdrop = 0;
	//clog<<"ckt, nodes: "<<get_name()<<" "<<ckt_nodes.size()<<endl;
	for(size_t i=0;i<ckt_nodes.size();i++){
		time = 0;
		j=0;
		while(time < tran.tot_t){// && iter <1){
			if(VDD - ckt_nodes[i].value[j] > max_IRdrop)
				max_IRdrop = VDD - ckt_nodes[i].value[j];
			j++;
			time += tran.step_t;
		}
	}
	return max_IRdrop;
}

void Circuit::save_ckt_nodes_to_tr(Tran &tran){
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

double Circuit::locate_maxIRdrop(){
	max_IRdrop = 0;
	for(size_t i=0;i<replist.size();i++){
		double IR_drop = VDD - replist[i]->value;		
		if(IR_drop > max_IRdrop)
			max_IRdrop = IR_drop;
	}
	return max_IRdrop;
}

double Circuit::locate_maxIRdrop(double *x, size_t n){
	max_IRdrop = 0;
	for(size_t i=0;i<n;i++){
		double IR_drop = VDD - x[i];		
		if(IR_drop > max_IRdrop)
			max_IRdrop = IR_drop;
	}
	return max_IRdrop;
}

double Circuit::locate_special_maxIRdrop(){
	double max_IRdrop = 0;
	Node *nd;
	for(size_t j=0;j<special_nodes.size(); 
			j++){
		nd = special_nodes[j];
		double IR_drop = VDD - nd->value;	
		if(IR_drop > max_IRdrop)
			max_IRdrop = IR_drop;
	}
	return max_IRdrop;
}

void Circuit::relocate_pads_graph(){
	vector<Node*> pad_set_old;
	double dist = 0;
	//double new_dist = 0;
	pad_set_old.resize(pad_set.size());
	assign_pad_set(pad_set_old);
	// store a original copy of the pad set
	
	// build up the map for nodes in PAD layer
	build_map_node_pt();
	
	vector<double> ref_drop_vec;
	//print_pad_set();
	for(size_t i=0;i<12;i++){
		clog<<"iter for pad move. "<<i<<endl;
		int pad_number = 1;
		origin_pad_set.resize(pad_set.size());
		assign_pad_set(origin_pad_set);
		// build pad connection graph
		build_graph();
		// find control nodes for each pad
		extract_pads(pad_number);
		// find the tune spot for control nodes	
		update_pad_control_nodes(ref_drop_vec, i);
		//print_all_control_nodes();	
		if(i>=6)
			dynamic_update_violate_ref(ref_drop_vec);
		// find new point for all pads	
		dist = update_pad_pos_all(ref_drop_vec);
		// move the low 10% pads into high 10% 
		// pad area 
		if(i==0)
			extract_min_max_pads(ref_drop_vec);
		// update the old pad set value
		assign_pad_set(pad_set_old);
		
		move_violate_pads(ref_drop_vec);
		
		// actual move pads into the new spots
		//project_pads();

		// move pads according to graph contraints
		graph_move_pads(ref_drop_vec);	
		
		clear_flags();
		// actual move pads into the new spots
		// project_pads();
		
		double max_IR = resolve_direct();
		if(max_IR ==0)
			break;
		//resolve_queue(origin_pad_set);
		//solve_GS();
		//clog<<"max_IRS is: "<<max_IRS<<endl<<endl;
	}
	ref_drop_vec.clear();
	map_node_pt.clear();
	origin_pad_set.clear();
	pad_set_old.clear();
	//print_pad_set();
	cout<<nodelist<<endl;
}

void Circuit::assign_pad_set(vector<Node*>&pad_set_old){
	//clog<<"assign pad set."<<endl;
	for(size_t i=0;i<pad_set_old.size();i++){
		pad_set_old[i] = pad_set[i]->node;
		//clog<<"pad: "<<i<<" "<<*pad_set_old[i]<<endl;	
	}
}
// special nodes includes nodes in local grids (including X node for LDO)
void Circuit::mark_special_nodes(){
	special_nodes.clear();
	/*int type = CURRENT;
	Net *net;
	Node *nd;
	// push back all current nodes
	for(size_t i=0;i<net_set[type].size();i++){
		net = net_set[type][i];
		nd = net->ab[0];
		if(nd->is_ground())
			nd = net->ab[1];
		special_nodes.push_back(nd);
	}*/
	//int z;
	//vector<int>::iterator it;
	for(size_t i=0;i<replist.size();i++){
		/*z = replist[i]->pt.z;
		it = find(local_layers.begin(), local_layers.end(), z);
		if(it == local_layers.end())
			continue;*/
		// only consider local grid nodes
		if(replist[i]->isS()!=X)
			special_nodes.push_back(replist[i]);
	}
}

double Circuit::update_pad_pos_all(vector<double> ref_drop_vec){
	double total_dist = 0;
	double dist = 0;
	for(size_t i=0;i<pad_set.size();i++){
		if(pad_set[i]->control_nodes.size()==0)
			continue;

		double ref_drop_value = ref_drop_vec[i];

		dist = update_pad_pos(ref_drop_value, i);
		total_dist += dist;
	}
	return total_dist;
}

void Circuit::modify_newxy(){
	Pad *pad;
	Node *nd;
	double delta_x;
	double delta_y;
	for(size_t i=0;i<pad_set.size();i++){
		pad = pad_set[i];
		nd = pad->node;
		if(pad->newx == nd->pt.x &&
		   pad->newy == nd->pt.y) 
			continue;
		delta_x = pad->newx - nd->pt.x;
		delta_y = pad->newy - nd->pt.y;
		delta_x /=2;
		delta_y /=2;
		pad->newx += delta_x;
		pad->newy += delta_y;
		
		round_data(pad->newx);
		round_data(pad->newy);
	}
}
// decide pad's new pos with the weights
// need to be tuned
double Circuit::update_pad_pos(double ref_drop_value, size_t i){
	//double total_dist=0;
	Node *pad;
	Pad *pad_ptr;
	Node *nd;
	double weight = 0;
	//double distance = 0;
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
		 	  
			if((pad_set[i]->node->name == "n0_30_74" ||
			    pad_set[i]->node->name == "n0_135_104"||
			    pad_set[i]->node->name == "n0_255_59")){
				//clog<<"data: "<<pad_set[i]->data<<endl;
			  //cout<<"control node: "<<*it->first<<" "<<it->second<<endl;
			 printf("%ld %ld  %.5e\n", it->first->pt.y+1, it->first->pt.x+1, it->first->value);
			
			}
			nd = it->first;
			weight = 1.0/it->second;
			weighted_x += weight * nd->pt.x;
			weighted_y += weight * nd->pt.y;
			sum_weight += weight; 	
		}
		
		if(sum_weight !=0){
			pad_newx = weighted_x / sum_weight;
			pad_newy = weighted_y / sum_weight;

			// pad_newx = (pad_newx - pad->pt.x)/2+pad->pt.x;
			//pad_newy = (pad_newy - pad->pt.y)/2+pad->pt.y;
			round_data(pad_newx);
			round_data(pad_newy);

			//if(pad->name == "n0_67_159")
				//clog<<"band pad, new: "<<*pad_ptr->node<<" "<<pad_newx<<" "<<pad_newy<<endl;

			if((pad_ptr->node->pt.x > 300 || pad_ptr->node->pt.y > 150) && (pad_newx <= 300 && pad_newy <= 150)){
				//clog<<"band pad, new: "<<*pad_ptr->node<<" "<<pad_newx<<" "<<pad_newy<<endl;
				pad_ptr->control_nodes.clear();
				pad_ptr->visit_flag = true;
				
			}else{
				

				pad_ptr->newx = pad_newx;
				pad_ptr->newy = pad_newy;
			}
		}else{
			pad_ptr->newx = pad->pt.x;
			pad_ptr->newy = pad->pt.y;
		}
 
		double dist = sqrt(weighted_x*weighted_x 			 + weighted_y*weighted_y);
		//total_dist += temp;
	//}

	//clog<<"dist: "<<total_dist<<endl<<endl;
	return dist;
}

void Circuit::project_pads(){
	Node *pad;
	Node *new_pad;
	Pad *pad_ptr;

	for(size_t i=0;i<pad_set.size();i++){
		if(pad_set[i]->control_nodes.size()==0)
			continue;
		pad_ptr = pad_set[i];
		pad = pad_ptr->node;
		// search for nearest node to (pad_newx,
		// pad_newy)
		
		new_pad = pad_projection(pad_ptr, pad);
		//cout<<"old pad / new pad: "<<*pad<<" "<<*new_pad<<endl;
		// update pad information
		pad_ptr->node = new_pad;
		pad_ptr->control_nodes.clear();	
	}
}

// round double, depends on whether the frac >=0.5
void Circuit::round_data(double &data){
	double fractpart, intpart;
	fractpart = modf(data, &intpart);
	if(fractpart >= 0.5)
		data = ceil(data);
	else
		data = floor(data);
}

// expand from (x,y) to nearest node in grid
// fits for non-uniform grid
Node * Circuit::pad_projection(Pad *pad, Node *nd){
	//Node *nd;
	queue<Point> q;
	Point pt;
	Point pt_cur;
	stringstream sstream;
	string pt_name;
	Node *nd_new=NULL;
	//int gap - 10;
	double dx[4] = {30, 0, -30, 0};
	double dy[4] = {0, 30, 0, -30};

	//nd = pad->node;
	pt.z = nd->get_layer();
	pt.x = pad->newx;
	pt.y = pad->newy;
	if(pad->newx == nd->pt.x && 
		pad->newy == nd->pt.y){
		nd_new = nd;	
		return nd_new;
	}

	sstream<<pt.z<<"_"<<pt.x<<"_"<<pt.y; 
	pt_name = sstream.str();
	// clog<<"pt_name: "<<pt_name<<endl;
	// first see if this node is on grid
	// and if it is occupied by pad or not
	//clog<<"orig pad: "<<*nd<<endl;
	if(has_node_pt(pt_name)){
		nd_new = get_node_pt(pt_name);
		 //cout<<"has new node: "<<*nd_new<<" "<<nd_new->isS()<<endl;
		// if this node is not occupied by pad
		if(nd_new->isS()!=X){
			nd->disableX();
			nd->value = 0;
			//nd_new->enableX();
			//nd_new->value = VDD;
			return nd_new;
		}
	}
	bool return_flag = false;
	// else start to search for node
	q.push(pt);
	// if not, expand it to neighboring area
	while(!q.empty()&& return_flag == false){
		pt_cur = q.front();
		Point pt_nbr = pt_cur;
		//expand_pad_pos(q, pt_cur);	
		for(size_t i=0;i<4;i++){
			pt_nbr.x = pt_cur.x + dx[i];
			pt_nbr.y = pt_cur.y + dy[i];
			stringstream sstream;
			string pt_name;
			sstream << pt_nbr.z<<"_"<<
				pt_nbr.x<<"_"<<
				pt_nbr.y;
			pt_name = sstream.str();
			if(has_node_pt(pt_name)){
				nd_new = get_node_pt(pt_name);

				//cout<<"new name: "<<*nd_new<<" "<<nd_new->isX()<<endl;
				if(nd_new->isS()!=X){
					nd->disableX();
					nd->value = 0;
					//nd_new->enableX();
					//nd_new->value = VDD;
					return_flag = true;
					break;
				}
			}
			q.push(pt_nbr);
		}
		q.pop();
	}
	while(!q.empty()){
		q.pop();
	}
	if(return_flag == true)
		return nd_new;
	clog<<"no point for new pad. return. "<<endl;
	return NULL;
}

void Circuit::build_map_node_pt(){
	if(pad_set.size()==0)
		clog<<"no pad on grid. ERROR"<<endl;
	// ref layer
	int ref_layer = pad_set[0]->node->get_layer();

	Node *nd;
	pair<string, Node*> pt_pair;
	for(size_t i=0;i<nodelist.size();i++){
		nd = nodelist[i];
		if(nd->get_layer()!=ref_layer)
			continue;
		if(nd->isS()==Z || nd->isS() ==X)
			continue;
		stringstream sstream;
		sstream<<ref_layer<<"_"<<nd->pt.x<<
			"_"<<nd->pt.y;
		pt_pair.first = sstream.str();
		//cout<<"string: "<<pt_pair.first<<endl;
		pt_pair.second = nd;
		map_node_pt.insert(pt_pair);
	}
}

// need to be modified: - rebuild voltage nets
void Circuit::restore_pad_set(vector<Node*>&pad_set_old){
	Node *nd_old=NULL;
	Node *nd_new=NULL;
	for(size_t i=0;i<pad_set_old.size();i++){
		// have to use find
		if(pad_set[i]->node !=
				pad_set_old[i]){
			nd_old = pad_set_old[i];
			nd_new = pad_set[i]->node;
			//clog<<"nd_old / nd_new: "<<
			//*nd_old<<" "<<*nd_new<<endl;
			nd_new->disableX();
			nd_new->value = 0;
			//nd_old->enableX();
			//nd_old->value = VDD;		
		}
		pad_set[i]->node = nd_old;	
	}
}

// modify voltage net set with rm_node and add_node
void Circuit::rebuild_voltage_nets(){
	int type = VOLTAGE;
	size_t index_rm_net = 0;
	size_t index_rm_resistor_net = 0;
	Net *net=NULL;
	Net *add_net=NULL;
	Net *resistor_net = NULL;
	Net *add_resistor_net = NULL;
	int type_1 = RESISTOR;
	Node *add_node_new;
	stringstream sstream;
	//Node *nd_ori=NULL;
	//Node *nd_new=NULL;
	Node *rm_node=NULL;
	Node *add_node=NULL;
	vector<Net*> rm_net;

	/*for(size_t i=0;i<origin_pad_set.size();i++){
		cout<<"add: "<<*pad_set[i]->node<<endl;
	}*/
	
	// delete all origin pad set
	// and build nets of new pad set
	for(size_t i=0;i<origin_pad_set.size();i++){
		rm_node = origin_pad_set[i];
		//rm_node = pad_set_old[i];
		add_node = pad_set[i]->node;
		if(rm_node == add_node || add_node->isS()==X)
			continue;

		//cout<<"rm_nod, add_node: "<<*rm_node<<" "<<*add_node<<endl;
		for(size_t i=0;i<net_set[type].size();i++){
			net = net_set[type][i];
			//cout<<"id, net: "<<net->id<<" "<<*net<<endl;

			if(net->ab[0]->name == rm_node->name ||
				net->ab[1]->name == rm_node->name){
				//clog<<endl<<"candi net: "<<*net<<endl;
				index_rm_net = net->id;
				rm_net.push_back(net);
				
				// also remove the resistor net
				if(net->ab[0]->name == rm_node->name){
					resistor_net = net->ab[0]->nbr[BOTTOM];
					rm_net.push_back(resistor_net);
					index_rm_resistor_net = resistor_net->id;
					net->ab[0]->nbr[BOTTOM] = NULL;
					resistor_net->ab[0]->nbr[TOP] = NULL;
				}
				else if(net->ab[1]->name == rm_node->name){
					resistor_net = net->ab[1]->nbr[BOTTOM];
					rm_net.push_back(resistor_net);
					index_rm_resistor_net = resistor_net->id;
					net->ab[1]->nbr[TOP] = NULL;
					resistor_net->ab[0]->nbr[TOP] = NULL;
				}
				break;
			}
		}
		// stores the X node
		// new node id and rep  =  old one
		add_node_new = new Node(*add_node);
		add_node_new->enableX();
		add_node_new->value = VDD;

		sstream.str("");
		// modify add_node_new with _X_xxxx
		sstream<<"_X_"<<add_node_new->name;
		add_node_new->name = sstream.str();

		add_net = new Net(VOLTAGE, VDD, add_node_new, 
				nodelist[nodelist.size()-1]);
		add_net->id = index_rm_net;

		net_set[type][index_rm_net] = add_net;

		double value = resistor_net->value;
		add_resistor_net = new Net(RESISTOR, value, add_node, add_node_new);
		add_resistor_net->id = index_rm_resistor_net;

		// modify the neighboring nets
		add_net->ab[0]->nbr[BOTTOM] = add_resistor_net;
		add_resistor_net->ab[0]->nbr[TOP] = add_resistor_net;
	
		type_1 = RESISTOR;
		net_set[type_1][index_rm_resistor_net] = add_resistor_net;

		// substitue the new pos with the old one
		// may need to free rm node
		//cout<<"rm_id: node, add_node: "<<rm_node->rid<<" "<<*rm_node<<" "<<*add_node_new<<endl;
		nodelist[rm_node->id] = add_node_new;
	}

	for(size_t i=0;i<rm_net.size();i++){
		delete rm_net[i];
	}
	origin_pad_set.clear();
	//free(rm_node);
}

void Circuit::print_pad_set(){
	for(size_t i=0;i<pad_set.size();i++){
		clog<<"pad: "<<*pad_set[i]->node<<endl;
		//printf("%d %d\n", pad_set[i]->node->pt.x+1, pad_set[i]->node->pt.y+1);
	}		
}

void Circuit::extract_min_max_pads_new(vector<double> ref_drop_vec){
	Node *nd;
	Pad *pad;
	size_t max_index = 0;
	double max = 0;

	vector<Node*> min_pads;
	vector<Node*> max_pads;
	vector<bool >temp_flag;

	map<Node *, double>::iterator it;
	temp_flag.resize(pad_set.size());
	for(size_t i=0;i<temp_flag.size();i++)
		temp_flag[i] = false;
	size_t id_minpad = 0;
	double drop = 0;
	double avg_ref = calc_avg_ref(ref_drop_vec);
	double avg_drop = VDD - avg_ref;
	
	for(size_t i=0;i<pad_set.size();i++){
		pad = pad_set[i];
		nd = pad->node;
		drop = VDD - ref_drop_vec[i];
		if(drop>max){
			max = drop;
			max_index = i;
		}
		if(drop < 0.7*avg_drop){
			min_pads.push_back(nd);	
		}	
	}
	double max_id;
	do{
		double max_temp = 0;
		max_id = -1;
		for(size_t j=0;j<ref_drop_vec.size();j++){
			if(VDD - ref_drop_vec[j] < max*0.9)
				continue;
			if(temp_flag[j] ==  false){
				if(max_id == -1)
					max_id = j;
				if(VDD - ref_drop_vec[j] > max_temp){
					max_temp = VDD - ref_drop_vec[j];
					max_id = j;
				}
			}
		}
		if(max_id == -1) break;
		temp_flag[max_id] = true;
		if(max_temp >= max*0.9)	
			max_pads.push_back(pad_set[max_id]->node);
	}while(max_id != -1);

	temp_flag.clear();
	Node *new_pad;
	Pad * pad_ptr;

	// set nd into the weighted center
	// start to map min pads into max pads locations
	for(size_t j=0;j<min_pads.size();j=j+2){
		for(size_t k=0;k<pad_set.size();k++){
			if(pad_set[k]->node->name == 
				min_pads[j]->name){
				id_minpad = k;
			}
		}
	
		size_t i = j % max_pads.size();
		//size_t i = locate_max_pad(max_pads, iter);
		Node * nd = max_pads[i];
		for(size_t k=0;k<pad_set.size();k++){
			if(pad_set[k]->node->name == 
				nd->name){	
				pad_ptr = pad_set[k];
				//double ref_drop_value = ref_drop_vec[k];

				new_pad = pad_projection(pad_ptr, min_pads[j]);
				//cout<<"old pad / new pad: "<<*min_pads[j]<<" "<<*new_pad<<endl;

				size_t m = id_minpad;
				pad_set[m]->node->disableX();
				pad_set[m]->node->value = 0;
				pad_set[m]->node = new_pad;
				pad_set[m]->visit_flag = true;
				// already taken care of
				pad_set[m]->control_nodes.clear();
				break;
			}
		}

	}
	// next step is to insert the min pads into max pads area
	min_pads.clear();
	max_pads.clear();
	
}

void Circuit::extract_min_max_pads(vector<double> ref_drop_vec){
	double min = 0;
	double max = 0;	
	size_t min_index=0;
	size_t max_index=0;
	vector<Node*> min_pads;
	vector<Node*> max_pads;
	size_t id_minpad = 0;

	for(size_t i=0;i<ref_drop_vec.size();i++){
		if(ref_drop_vec[i] == 0)
			continue;
		min = VDD - ref_drop_vec[i];
		min_index = i;
		break;
	}
	
	for(size_t j=0;j<ref_drop_vec.size();j++){	
		if(pad_set[j]->control_nodes.size()==0)
			continue;
		if(VDD-ref_drop_vec[j]>max){
			max = VDD - ref_drop_vec[j];
			max_index = j;
		}
		if(ref_drop_vec[j]!=0 && VDD-ref_drop_vec[j]<min){
			min = VDD - ref_drop_vec[j];
			min_index=j;
		}
	}

	//cout<<"min, max: "<<min<<" "<<max<<endl;
	for(size_t j=0;j<ref_drop_vec.size();j++){
		// skip if a pad has no control nodes
		if(pad_set[j]->control_nodes.size()==0)
			continue;
		if(VDD - ref_drop_vec[j]<= max*0.5){
			min_pads.push_back(pad_set[j]->node);
		}
		if(VDD-ref_drop_vec[j]>=max*0.8){
			max_pads.push_back(pad_set[j]->node);
		}
	}
	
	Node *new_pad;
	Pad * pad_ptr;

	// set nd into the weighted center
	// start to map min pads into max pads locations
	for(size_t j=0;j<min_pads.size();j++){
		for(size_t k=0;k<pad_set.size();k++){
			if(pad_set[k]->node->name == 
				min_pads[j]->name){
				id_minpad = k;
			}
		}
		size_t i = j % max_pads.size();
		Node * nd = max_pads[i];
		for(size_t k=0;k<pad_set.size();k++){
			if(pad_set[k]->node->name == 
				nd->name){	
				pad_ptr = pad_set[k];
				//double ref_drop_value = ref_drop_vec[k];
				new_pad = pad_projection(pad_ptr, min_pads[j]);
				
				//cout<<"old pad / new pad: "<<*min_pads[j]<<" "<<*new_pad<<endl;
				size_t m = id_minpad;		
				pad_set[m]->node->disableX();
				pad_set[m]->node->value = 0;
				pad_set[m]->node = new_pad;
				//pad_set[m]->visit_flag = true;
				// already taken care of
				pad_set[m]->control_nodes.clear();
				break;
			}
		}

	}
	/*for(size_t j=0;j<min_pads.size();j++){
		cout<<"min_pads: "<<*min_pads[j]<<endl;
	}
	for(size_t j=0;j<max_pads.size();j++){
		cout<<"max_pads: "<<*max_pads[j]<<endl;
	}*/

	min_pads.clear();
	max_pads.clear();
}

void Circuit::move_violate_pads(vector<double> ref_drop_vec){
	Pad *pad_ptr;
	Node * pad;
	Node * new_pad;
	for(size_t i=0;i<pad_set.size();i++){
		if(pad_set[i]->control_nodes.size()==0)
			continue;

		pad_ptr = pad_set[i];
		pad = pad_ptr->node;
		// if violate, move this pad
		if(pad_ptr->violate_flag == true){
			new_pad = pad_projection(pad_ptr, pad);
			//cout<<"old pad / new pad: "<<*pad<<" "<<*new_pad<<endl;
			pad_ptr->node = new_pad;
			pad_ptr->control_nodes.clear();
			pad_ptr->visit_flag = true;
		}
	}
}

void Circuit::graph_move_pads(vector<double> ref_drop_vec){
	Node *new_pad;
	int id=0;
	//for(size_t i=0;i<5;i++){
	do{
		id = locate_max_drop_pad(ref_drop_vec);
		if(id==-1) break;
		Pad *pad_ptr = pad_set[id];
		Pad *pad_nbr = NULL;
		Node *pad = pad_ptr->node;
		//clog<<endl<<"pad: "<<*pad<<endl;
		new_pad = pad_projection(pad_ptr, pad);
		
		//cout<<"old pad / new pad: "<<*pad<<" "<<*new_pad<<endl;

		//bool flag = print_flag(pad);
		//if(flag == true || new_pad->name == "n0_86_30" || new_pad->name =="n0_477_17" || pad->name == "n0_477_17")
			//clog<<"old pad / new pad: "<<*pad<<" "<<*new_pad<<endl;

		//clog<<"new_pad: "<<*new_pad<<endl;
		pad_ptr->visit_flag = true;
		for(size_t i=0;i<pad_ptr->nbrs.size();i++){
			pad_nbr = pad_ptr->nbrs[i];
			if(pad_nbr->fix_flag == false){
				pad_nbr->fix_flag = true;
				//clog<<"pad_nbr: "<<*pad_nbr->node<<endl;
			}
		}

		pad_ptr->node = new_pad;
		pad_ptr->control_nodes.clear();
	}while(id != -1);
}

// locate id that has minimum value and not visited or fixed 
int Circuit::locate_max_drop_pad(vector<double> vec){
	int min_id = -1;
	double min_ref = 0;
	bool flag = false;
	for(size_t i=0;i<vec.size();i++){
		if(pad_set[i]->control_nodes.size()==0)
			continue;
		if(pad_set[i]->visit_flag == true ||
			pad_set[i]->fix_flag ==  true)
			continue;
		//clog<<"i, vec: "<<i<<" "<<vec[i]<<endl;
		if(flag == false){
			flag = true;
			min_ref = vec[i];
			min_id = i;
		}
		else if(vec[i] < min_ref){
			min_ref = vec[i];
			min_id = i;
		}
	}	
	//clog<<"min_ref, min_pad: "<<min_ref<<" "<<*pad_set[min_id]->node<<endl;
	return min_id;
}

// locate the tune spot for the control nodes.
double Circuit::locate_ref(size_t i){
	Pad *pad_ptr;
	Node *pad;
	map<Node*, double>::iterator it;
	Node *nd;
	double weight = 0;
	//vector<double> drop_vec;
	pad_ptr = pad_set[i];
	pad = pad_ptr->node;
	pad_ptr->drop_vec.clear();
	// cout<<"pad: "<<*pad<<endl;
	for(it = pad_ptr->control_nodes.begin();
			it != pad_ptr->control_nodes.end();
			it++){
		nd = it->first;
		// cout<<"control: "<<*nd<<endl;
		weight = nd->value;
		if(weight <0)
			weight *=10;

		pad_ptr->control_nodes[nd] = weight;
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

void Circuit::dynamic_update_violate_ref(vector<double> & ref_drop_vec){
	//for(size_t j=0;j<2;j++){
	double avg_drop = calc_avg_ref_drop(ref_drop_vec);
	Pad *pad_ptr;
	Node *pad;
	//cout<<"j: "<<j<<endl;
	for(size_t i=0;i<pad_set.size();i++){
		if(pad_set[i]->control_nodes.size()==0)
			continue;
		pad_ptr = pad_set[i];
		pad = pad_ptr->node;

		if(pad_ptr->data >= 2*avg_drop){
			pad_ptr->violate_flag = true;
			double ratio_new = pad_ptr->ratio * 2;
			size_t id = pad_ptr->drop_vec.size() / ratio_new;
			pad_ptr->ratio = ratio_new;
			double middle_value = pad_ptr->drop_vec[id];
			ref_drop_vec[i] = middle_value;
		}
	}
	
	extract_min_max_pads_new(ref_drop_vec);	
}

bool Circuit::print_flag(Node *pad){
	bool flag = false;
	if(pad->name == "n0_135_104" ||
	   pad->name == "n0_255_59"||
	   pad->name == "n0_30_74")
	   /*pad->name == "n0_67_148" ||
	   pad->name == "n0_114_171" ||
	   pad->name == "n0_162_159" ||
	   pad->name == "n0_216_171" ||
	   pad->name == "n0_268_177")*/
		flag = true;
	return flag;
}

double Circuit::calc_avg_ref_drop(vector<double> &ref_drop_vec){
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

double Circuit::calc_avg_ref(vector<double> ref_drop_vec){
	//Node *pad;
	//Pad *pad_ptr;
	double sum_ref = 0;
	size_t count = 0;
	double ref_drop_value = 0;

	for(size_t i=0;i<pad_set.size();i++){
		if(pad_set[i]->control_nodes.size()==0)
			continue;

		count ++;	
		ref_drop_value = ref_drop_vec[i];
		sum_ref += ref_drop_value;
	}
	double avg_ref = sum_ref / count;
	return avg_ref;
}

void Circuit::clear_flags(){
	Pad *pad;
	for(size_t i=0;i<pad_set.size();i++){
		pad = pad_set[i];
		pad->visit_flag = false;
		pad->fix_flag = false;
		pad->violate_flag = false;
		//nd->region_flag = false;
		//nd->fix_flag = false;
		//nd->visit_flag = false;
	}
}

double Circuit::resolve_direct(){
	clock_t t1, t2;
	t1 = clock();
	//cout<<endl;
	//cout<<"============ a new round ======"<<endl;
	rebuild_voltage_nets();	
	Net *net;	

	solve_DC();
	//solve_LU_core();
	double max_IR = locate_maxIRdrop();	
	//double max_IRS = locate_special_maxIRdrop();
	clog<<"max_IR by cholmod is: "<<max_IR<<endl;
	t2 = clock();
	//clog<<"single solve by cholmod is: "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;
	return max_IR;
}

// build graph for pad nodes
void Circuit::build_graph(){
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
		pad_nbr = find_nbr_pad(pad);
		//cout<<"pad, nbr: "<<*pad->node<<" "<<*pad_nbr->node<<endl;
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

void Circuit::extract_pads(int pad_number){
	vector<Node*> pair_first;
	vector<double> pair_second;
	pair<Node*, double> pair_nd;
	//int pad_number = 5;
	double distance = 0;
	map<Node*, double>::iterator it;

	clear_pad_control_nodes();
	for(size_t i=0;i<special_nodes.size();i++){
		int count = 0;
		pair_first.clear();
		pair_second.clear();
		Node *nd = special_nodes[i];
		// search for closest pads
		for(size_t j=0;j<pad_set.size();j++){
			Node *ptr = pad_set[j]->node;
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

// use Euclidiean distance to locate nearest nbr pad
Pad * Circuit::find_nbr_pad(Pad *pad){
	Pad * nbr;
	double distance=-1;
	double min_dist=0;
	bool flag = false;
	size_t min_index=0;
	for(size_t i=0;i<pad_set.size();i++){
		nbr = pad_set[i];
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

double Circuit::get_distance(Node *na, Node *nb){
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

void Circuit::clear_pad_control_nodes(){
	for(size_t i=0;i<pad_set.size();i++){
		pad_set[i]->control_nodes.clear();
	}
}

// tune 50% nodes with the small IR drops
void Circuit::update_pad_control_nodes(vector<double> & ref_drop_value, size_t iter){
	ref_drop_value.resize(pad_set.size());
	for(size_t i=0;i<pad_set.size();i++){
		if(pad_set[i]->control_nodes.size()==0)
			continue;
		
		double middle_value = locate_ref(i);
		ref_drop_value[i] = middle_value;
		//cout<<"middle value: "<<middle_value<<endl;
	}
}

void Circuit::build_pad_set(){
	pad_set.resize(0);
	//static Pad pad_a;
	for(size_t i=0;i<nodelist.size()-1;i++){
		if(nodelist[i]->isS()==X){
			Pad *pad_ptr = new Pad();
			pad_ptr->node = nodelist[i];
			pad_set.push_back(pad_ptr);
		}
	}
	//for(size_t j=0;j<pad_set.size();j++)
		//cout<<"pad: "<<*pad_set[j]->node<<endl;
}

void Circuit::print_pad_map(){
	Node *nd;
	Pad *pad;
	pair<Node*, double> pad_pair;
	map<Node*, double>::iterator it;
	for(size_t i=0;i<pad_set.size();i++){
		nd = pad_set[i]->node;
		pad = pad_set[i];
		if(nd->name == "n0_477_17"){// ||
		   //nd->name == "n0_255_59" ||
		   //nd->name == "n0_30_74"){
		for(it = pad->control_nodes.begin();it != pad->control_nodes.end();it++){
			printf("%ld %ld  %.5e\n", it->first->pt.y+1, it->first->pt.x+1, it->first->value);
		}
		}
	}
}

void Circuit::get_pad_tr_cur(Tran &tran){
	Net *net;
	Node *na, *nb;
	size_t p=0;
	double time = 0;
	while(time < tran.tot_t){
		// insert all the pad node into tran object
		for(size_t i=0;i<pad_set.size();i++){
			net = pad_set[i]->node->nbr[BOTTOM];
			if(net == NULL) clog<<"no bottom net for pad. ERROR. "<<endl;
			p = 0;
			if(net->ab[0]->isS()==X)
				p = 1;
			na = net->ab[p];
			double current = 0;
			// add neighboring nodes of pad nodes into the list
			for(DIRECTION d=WEST;d<=NORTH;d=DIRECTION(d+1)){
				net = na->nbr[d];
				if(net == NULL) continue;
				nb = net->ab[0];
				if(nb->name == na->name)
					nb = net->ab[1];
				current += fabs(na - nb)/net->value;
			}
			// add bottom current if there is one
			if(na->nbr[BOTTOM]!=NULL){
				net = na->nbr[BOTTOM];
				if(net->type == CURRENT){
					current += net->value;
				}
			}
			pad_set[i]->current.push_back(current);
		}
		time += tran.step_t;
	}
}


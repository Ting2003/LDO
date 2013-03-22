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
	/*for(size_t i=0;i<ldolist.size();i++){
		Node *nd = ldolist[i]->A;

		clog<<"i, ldo_A, ldo_in: "<<i<<" "<<*ldolist[i]->A<<" "<<nd->is_LDO()<<" "<<*ldolist[i]->nd_in<<endl;
	}*/
	//clog<<"before build pad set. "<<endl;
	if(pad_set.size()==0)	
		build_pad_set();

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
			else if((p->nbr[NORTH] !=NULL && p->nbr[SOUTH] !=NULL)			      &&(p->nbr[EAST]==NULL && p->nbr[WEST] ==NULL))
				count++;
		}
	}
	//clog<<"number of nodes can be merged is: "<<count<<endl;
}

void SubCircuit::solve(Tran &tran, bool flag){
	solve_LU(tran, flag);
}


// stamp the matrix and solve
void SubCircuit::solve_LU_core(Tran &tran, bool local_flag){
   size_t n = replist.size();	// replist dosn't contain ground node
   if( n == 0 ) return;		// No node  

   // configure_init();
   cm = &c;
   cholmod_start(cm);
   cm->print = 5;
   b = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
   x = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
   bp = static_cast<double *> (b->x);

   // Matrix A;
   A.clear();
   stamp_by_set(A, bp);

   if(local_flag == false)
   	make_A_symmetric(bp);
   else{
	make_A_symmetric_local(bp);
   }

   A.set_row(n);
   /*cout<<A<<endl;
   for(size_t i=0;i<n;i++)
	   cout<<"i, b: "<<i<<" "<<bp[i]<<endl;*/
   Algebra::solve_CK(A, L, x, b, cm);
   //return;
   xp = static_cast<double *> (x->x);
   /*for(size_t i=0;i<n;i++)
	cout<<"dc solution: "<<
	get_name()<<" "<<xp[i]<<endl;*/
   
   get_voltages_from_LU_sol(xp);
   // print out dc solution
   cout<<"dc solution. "<<endl;
   cout<<nodelist<<endl;

   return; 
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
   make_A_symmetric_tr(bp, xp, tran);
   
   stamp_current_tr(bp, time);
  
   Algebra::CK_decomp(A, L, cm);
   Lp = static_cast<int *>(L->p);
   Lx = static_cast<double*> (L->x);
   Li = static_cast<int*>(L->i) ;
   Lnz = static_cast<int *>(L->nz); 
   A.clear();
//#if 0 
   /*********** the following 2 parts can be implemented with pthreads ***/
   // build id_map immediately after transient factorization
   id_map = new int [n];
   cholmod_build_id_map(CHOLMOD_A, L, cm, id_map);

   temp = new double [n];
   // then substitute all the nodes rid
   for(size_t i=0;i<n;i++){
	   int id = id_map[i];
	   replist[id]->rid = i;
	   temp[i] = bp[i];
   }
   for(size_t i=0;i<n;i++)
	   bp[i] = temp[id_map[i]];
   for(size_t i=0;i<n;i++)
	   temp[i] = xp[i];
   for(size_t i=0;i<n;i++)
	   xp[i] = temp[id_map[i]];
 
   delete [] temp;
   delete [] id_map;
   /*****************************************/ 
   for(size_t i=0;i<n;i++)
	bnewp[i] = bp[i];
     
   set_eq_induc(tran);
   set_eq_capac(tran);
   modify_rhs_tr_0(bnewp, xp);

   // push rhs node into node_set b
   for(size_t i=0;i<n;i++){
	   if(bnewp[i] !=0)
		   pg.node_set_b.push_back(i);
   }

   // push back all nodes in output list
   vector<size_t>::iterator it;
   size_t id;
   for(size_t i=0;i<tran.nodes.size();i++){
	   if(tran.nodes[i].node == NULL) continue;
	   if(!tran.nodes[i].node->rep->is_ground()){
		   id = tran.nodes[i].node->rep->rid;
		   it = find(pg.node_set_x.begin(), pg.node_set_x.end(), 
				   id);
		   if(it == pg.node_set_x.end()){
			   pg.node_set_x.push_back(id);
		   }
	   }
   }
   // get path_b, path_x, len_path_b, len_path_x
   build_path_graph();

   s_col_FFS = new int [len_path_b];
   s_col_FBS = new int [len_path_x];
   find_super();
   solve_eq_sp(xp);
 
   //save_tr_nodes(tran, xp);
   save_ckt_nodes(tran, xp);
   time += tran.step_t;
   // clog<<"before tr solve."<<endl;
   // then start other iterations
   while(time < tran.tot_t){// && iter < 0){
	for(size_t i=0;i<n;i++)
		bnewp[i] = bp[i];
       // only stamps if net current changes
      // set bp into last state
      stamp_current_tr_1(bp, bnewp, time);
     // get the new bnewp
      modify_rhs_tr(bnewp, xp); 
	
      solve_eq_sp(xp);

      save_ckt_nodes(tran, xp);
      time += tran.step_t;
      //iter ++;
   }
   save_ckt_nodes_to_tr(tran);
   //print_ckt_nodes(tran);
   //release_tr_nodes(tran);
   release_ckt_nodes(tran);
   cholmod_free_dense(&b, cm);
   cholmod_free_dense(&bnew, cm);
   cholmod_free_factor(&L, cm);
   cholmod_free_dense(&x, cm);
   cholmod_finish(&c);
  
   delete [] s_col_FFS;
   delete [] s_col_FBS;
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
void SubCircuit::stamp_by_set(Matrix & A, double * b){
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
			for(size_t i=0;i<ns.size();i++)
				stamp_current(b, ns[i]);
			break;
		case VOLTAGE:
			for(size_t i=0;i<ns.size();i++){
				if( fzero(ns[i]->value)  && 
				    !ns[i]->ab[0]->is_ground() &&
				    !ns[i]->ab[1]->is_ground() )
					continue; // it's a 0v via
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
void SubCircuit::stamp_by_set_tr(Matrix & A, double *b, Tran &tran){
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
void SubCircuit::modify_rhs_tr_0(double * b, double *x){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		if(type ==CAPACITANCE){	
			for(size_t i=0;i<ns.size();i++)
				modify_rhs_c_tr_0(ns[i], b, x);
		}
		else if(type == INDUCTANCE){
			for(size_t i=0;i<ns.size();i++){
				modify_rhs_l_tr_0(ns[i], b, x);	
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

// only stamp the resistor node connected to inductance
void SubCircuit::stamp_resistor_tr(Matrix & A, Net * net){
   double G;
   Node * nk = net->ab[0]->rep;
   Node * nl = net->ab[1]->rep;
   size_t k = nk->rid;
   size_t l = nl->rid;
   G = 1./net->value;

   if( nk->isS()!=Y && !nk->is_ground()&& 
     (nk->nbr[TOP]!=NULL && 
      nk->nbr[TOP]->type == INDUCTANCE)) {
      //clog<<"net: "<<*net<<endl;
      A.push_back(k,k, G);
      //clog<<"("<<k<<" "<<k<<" "<<G<<")"<<endl;
      if(!nl->is_ground() && (nl->nbr[TOP]==NULL || 
           nl->nbr[TOP]->type != INDUCTANCE)){
         if(l < k){
            A.push_back(k,l,-G);
            //clog<<"("<<k<<" "<<l<<" "<<-G<<")"<<endl;
         }
         else if(l > k){ 
            A.push_back(l, k, -G);
            //clog<<"("<<l<<" "<<k<<" "<<-G<<")"<<endl;
         }
      }
   }

   if( nl->isS() !=Y && !nl->is_ground()&&
     (nl->nbr[TOP]!=NULL &&
      nl->nbr[TOP]->type == INDUCTANCE)) {

      //clog<<"net: "<<*net<<endl;
      A.push_back(l,l, G);
      //clog<<"("<<l<<" "<<l<<" "<<G<<")"<<endl;
      if(!nk->is_ground() && (nk->nbr[TOP]==NULL ||
           nk->nbr[TOP]->type != INDUCTANCE)){
         if(k < l){
            A.push_back(l,k,-G);
            //clog<<"("<<l<<" "<<k<<" "<<-G<<")"<<endl;
         }
         else if(k > l){
            A.push_back(k, l, -G);
            //clog<<"("<<k<<" "<<l<<" "<<-G<<")"<<endl;
         }
      }
   }
}

void SubCircuit::stamp_inductance_dc(Matrix & A, double *b, Net * net){
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
		if(!nl->is_ground())
		// A.push_back(k,l,-1);
		// make it symmetric
			b[k] = b[l];
		//clog<<"("<<k<<" "<<k<<" "<<1<<")"<<endl;
		//clog<<"("<<k<<" "<<l<<" "<<-1<<")"<<endl;
	}

	if( nl->isS() !=Y && !nl->is_ground()) {
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
	//clog<<"net: "<<*net<<endl;
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
		A.push_back(k,k, Geq-1);
		//clog<<"("<<k<<" "<<k<<" "<<Geq-1<<")"<<endl;
		//clog<<nl->isS()<<endl;
		if(!nl->is_ground()&& nl->isS()!=Y && k>l){
			A.push_back(k,l,-Geq);
		        //clog<<"("<<k<<" "<<l<<" "<<-Geq<<")"<<endl;
		}
	}

	if( nl->isS() !=Y && !nl->is_ground()) {
		// -1 is to clear formal inserted 1 at (l,l)
		A.push_back(l,l, Geq-1);
		//clog<<"("<<l<<" "<<l<<" "<<Geq-1<<")"<<endl;
		if(!nk->is_ground() && nk->isS()!=Y && l>k){
			A.push_back(l,k,-Geq);
			//clog<<"("<<l<<" "<<k<<" "<<-Geq<<")"<<endl;
		}
	}
}

// stamp capacitance Geq = 2C/delta_t
void SubCircuit::stamp_capacitance_tr(Matrix &A, Net *net, Tran &tran){
	//clog<<"net: "<<*net<<endl;
	double Geq = 0;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	// Geq = 2*C / delta_t
	Geq = (2*net->value) / tran.step_t;
	//net->value = Geq;
	//clog<<"C delta_t Geq: "<<net->value<<" "<<tran.step_t<<" "<<Geq<<endl;
	// Ieq = i(t) + 2*C / delta_t * v(t)

	if( nk->isS()!=Y  && !nk->is_ground()) {
		A.push_back(k,k, Geq);
		//clog<<"("<<k<<" "<<k<<" "<<Geq<<")"<<endl;
		if(!nl->is_ground()&& k > l){
			A.push_back(k,l,-Geq);
			//clog<<"("<<k<<" "<<l<<" "<<-Geq<<")"<<endl;
		}
	}

	if( nl->isS() !=Y && !nl->is_ground()) {
		A.push_back(l,l, Geq);
		//clog<<"("<<l<<" "<<l<<" "<<Geq<<")"<<endl;
		if(!nk->is_ground()&& l > k){
			A.push_back(l,k,-Geq);
			//clog<<"("<<l<<" "<<k<<" "<<-Geq<<")"<<endl;
		}
	}
}

// add Ieq into rhs
// Ieq = i(t) + 2*C / delta_t *v(t)
void SubCircuit::modify_rhs_c_tr_0(Net *net, double * rhs, double *x){
	double i_t = 0;
	double temp = 0;
	double Ieq = 0;
	//clog<<"c net: "<<*net<<endl;
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

	Net *r = nk->nbr[TOP];
	Node *a = r->ab[0]->rep;
	Node *b = r->ab[1]->rep;
	// a point to Z node
	if(a->isS()!=Z) {
		swap<Node *>(a, b);
		swap<Node*>(r->ab[0], r->ab[1]);
	}
	//clog<<"a, b: "<<*a<<" "<<*b<<endl;

	size_t id_a = a->rid;
	size_t id_b = b->rid;
	//i_t = (b->value - a->value) / r->value;
	i_t = (x[id_b] - x[id_a]) / r->value;
	//if(b->value != x[id_b] || a->value != x[id_a])
	   //cout<<"a, b, x_a, x_b: "<<a->value<<" "<<b->value<<" "<<
	     //x[id_a]<<" "<<x[id_b]<<endl;
	//clog<<"i_t: "<<i_t<<endl;
	//temp = 2*net->value / tran.step_t * 
		//(nk->value - nl->value);
       
        // push 2 nodes into node_set_x
        //clog<<*nk<<" "<<k<<endl;
//#if 0
        pg.node_set_x.push_back(k);
        if(!nl->is_ground()) {
              //clog<<*nl<<" "<<l<<endl;
           pg.node_set_x.push_back(l);
        }
        else if(!b->is_ground()){
              //clog<<*b<<" "<<id_b<<endl;
           pg.node_set_x.push_back(id_b);
        }
//#endif
	if(nk->is_ground())
	 //temp = 2*net->value/tran.step_t*(0-x[l]);
	 temp = net->value *(-x[l]);
        else if(nl->is_ground()){
         //temp = 2*net->value/tran.step_t *(x[k]);
	 temp = net->value *x[k];
        }
        else
         //temp = 2*net->value/tran.step_t *(x[k] - x[l]);
	 temp = net->value *(x[k]-x[l]);
	//if(nk->value != x[k] || nl->value != x[l])
	   //cout<<"k, l, x_k, x_l: "<<nk->value<<" "<<nl->value<<" "<<
	     //x[k]<<" "<<x[l]<<endl;
	//clog<<"nk-nl "<<(nk->value - nl->value)<<" "<<2*net->value/tran.step_t<<" "<<temp<<endl;
	
	Ieq  = (i_t + temp);
	//clog<< "Ieq is: "<<Ieq<<endl;
	//clog<<"Geq is: "<<2*net->value / tran.step_t<<endl;
	if(!nk->is_ground()&& nk->isS()!=Y){
		 rhs[k] += Ieq;	// for VDD SubCircuit
		//clog<<*nk<<" rhs +: "<<rhs[k]<<endl;
	}
	if(!nl->is_ground()&& nl->isS()!=Y){
		 rhs[l] += -Ieq; 
		//clog<<*nl<<" rhs +: "<<rhs[l]<<endl;
	}
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

void SubCircuit::set_eq_induc(Tran &tran){
	NetPtrVector &ns = net_set[INDUCTANCE];
	for(size_t i=0;i<ns.size();i++)
		ns[i]->value = tran.step_t /(2*ns[i]->value);
}

void SubCircuit::set_eq_capac(Tran &tran){
	NetPtrVector &ns = net_set[CAPACITANCE];
	for(size_t i=0;i<ns.size();i++)
		ns[i]->value = 2*ns[i]->value/tran.step_t;
}

// add Ieq into rhs
// Ieq = i(t) + delta_t / (2*L) *v(t)
void SubCircuit::modify_rhs_l_tr_0(Net *net, double *rhs, double *x){
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
        pg.node_set_x.push_back(k);
        pg.node_set_x.push_back(id_b);
//#endif
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
	//clog<<"net: "<<*net<<endl;
	//clog<<"current: "<<current<<endl;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	if( !nk->is_ground()&& nk->isS()!=Y) { 
		size_t k = nk->rid;
		//clog<<"node, rid: "<<*nk<<" "<<k<<endl;
		b[k] += -net->value;//current;
		//clog<<"time, k, b: "<<time<<" "<<k<<" "<<b[k]<<endl;
	}
	if( !nl->is_ground() && nl->isS()!=Y) {
		size_t l = nl->rid;
		//clog<<"node, rid: "<<*nl<<" "<<l<<endl;
		b[l] +=  net->value;// current;
		//clog<<"time, l, b: "<<time<<" "<<l<<" "<<b[l]<<endl;
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
void SubCircuit::stamp_VDD(Matrix & A, double * b, Net * net){
	// find the non-ground node
	// clog<<"net: "<<*net<<endl;
	Node * X = net->ab[0];
	if( X->is_ground() ) X = net->ab[1];
	size_t id = X->rep->rid;
	A.push_back(id, id, 1.0);
	//clog<<"push id, id, 1: "<<id<<" "<<id<<" "<<1<<endl;
	Net * south = X->rep->nbr[SOUTH];
	if( south != NULL &&
	    south->type == CURRENT ){
		// this node connects to a VDD and a current
		// assert( feqn(1.0, b[id]) ); // the current should be stamped
		b[id] = net->value;	    // modify it
		//clog<<"b: ="<<id<<" "<<net->value<<endl;
	}
	else{
		b[id] += net->value;
		//clog<<"b: +"<<id<<" "<<net->value<<endl;
	}
}

// stamp a voltage source
void SubCircuit::stamp_VDD_tr(double * b, Net * net){
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
void SubCircuit:: save_ckt_nodes(Tran &tran, double *x){
   size_t id=0;
   for(size_t j=0;j<ckt_nodes.size();j++){
	 //cout<<"nodes: "<<ckt_nodes[j].node->name<<endl;
         id = ckt_nodes[j].node->rep->rid;
	 //cout<<"value: "<<x[id]<<endl;
         ckt_nodes[j].value.push_back(x[id]);
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

void SubCircuit:: release_ckt_nodes(Tran &tran){
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

void SubCircuit::make_A_symmetric_tr(double *b, double *x, Tran &tran){
	int type = INDUCTANCE;
	NetList & ns = net_set[type];
	NetList::iterator it;
	Node *p, *q;

	for(it=ns.begin();it!=ns.end();it++){
           if( (*it) == NULL ) continue;
           assert( fzero((*it)->value) == false );
           if(!((*it)->ab[0]->rep->isS()==Y || (*it)->ab[1]->rep->isS()==Y)) continue;
           //clog<<"net: "<<*(*it)<<endl;
           // node p points to Y node
           if((*it)->ab[0]->rep->isS()==Y){
              p = (*it)->ab[0]->rep; q = (*it)->ab[1]->rep;
           } 
           else if((*it)->ab[1]->rep->isS()==Y ){
              p = (*it)->ab[1]->rep; q = (*it)->ab[0]->rep;
           }           
           //clog<<"p and q: "<<*p<<" "<<*q<<endl;

           size_t id = q->rid;
           double G = tran.step_t / ((*it)->value*2);
           
           //b[id] += p->value * G;
           b[id] += x[p->rid] *G;
           //clog<<"stamp p->value, G, b: "<<p->value<<" "<<G<<" "<<b[id]<<endl;
        }
}

bool compare_Node_G(const Node_G *nd_1, const Node_G *nd_2){
   return (nd_1->value < nd_2->value);
 }

void SubCircuit::update_node_set_bx(){
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


void SubCircuit::parse_path_table(){
   // build up nodelist info
      Node_G *node;
      for(size_t i=0;i<replist.size();i++){
         node = new Node_G();
         node->value = i;
         pg.nodelist.push_back(node);
      }
}

void SubCircuit::build_path_graph(){
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

void SubCircuit::build_FFS_path(){
   parse_path_table(); 

   set_up_path_table();

   find_path(pg.node_set_b, pg.path_FFS);

   pg.path_FFS.assign_size();

   for(size_t i=0;i<replist.size();i++)
      pg.nodelist[i]->flag = 0;
}

void SubCircuit::build_FBS_path(){
  pg.nodelist.clear();
  parse_path_table();
  set_up_path_table();
  find_path(pg.node_set_x, pg.path_FBS);
   pg.path_FBS.assign_size();
}

void SubCircuit::set_up_path_table(){
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

void SubCircuit::find_path(vector<size_t> &node_set, List_G &path){
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
void SubCircuit::solve_eq_set(){
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

void SubCircuit::solve_eq(double *X){
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
void SubCircuit::find_super(){
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
 
 void SubCircuit::solve_eq_sp(double *X){
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
 
       if (s_col_FFS[k]==1)//lnz < 4 || lnz != Lnz [j+1] + 1 || Li [p+1] != j+1)
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
       else if (s_col_FFS[k]==2)//lnz != Lnz [j+2] + 2 || Li [p+2] != j+2)
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
 
        if (s_col_FBS[k]==1)//j < 4 || lnz != Lnz [j-1] - 1 || Li [Lp [j-1]+1] != j)
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
       else if (s_col_FBS[k]==2)//lnz != Lnz [j-2]-2 || Li [Lp [j-2]+2] != j)
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
	// then create new ones
	Node *nd;
	for(size_t i=0;i<ldolist.size();i++){
		nd = ldolist[i]->A;
		nd->value = VDD_G;
		Net *net = new Net(VOLTAGE, VDD_G, nd, nodelist[nodelist.size()-1]);
		add_net(net);
		// update top nbr net
		nd->rep->nbr[TOP] = net;
		//clog<<"add net: "<<*net<<endl;
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
	double current;
	for(size_t i=0;i<ldolist.size();i++){
		nd = ldolist[i]->nd_in;
		current = ldolist[i]->current;
		//current = 0.3;
		Net *net = new Net(CURRENT, current, nd, nodelist[nodelist.size()-1]);
		add_net(net);
		// update top nbr net
		nd->rep->nbr[BOTTOM] = net;
		//clog<<"add net: "<<*net<<endl;
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

// stmap matrix and rhs, decomp matrix for DC
void SubCircuit::stamp_decomp_matrix_DC(bool local_flag){
   A.clear();
   size_t n = replist.size();
   b = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
   bp = static_cast<double *> (b->x);

   stamp_by_set(A, bp);

   if(local_flag == false)
   	make_A_symmetric(bp);
   else
	make_A_symmetric_local(bp);
   A.set_row(replist.size());
   Algebra::CK_decomp(A, L, cm);
   A.clear();
}

// stamp matrix and rhs, decomp matrix for TR
void SubCircuit::stamp_decomp_matrix_TR(Tran &tran, double time){
   A.clear();
   stamp_by_set_tr(A, bp, tran);
   make_A_symmetric_tr(bp, xp, tran);
   stamp_current_tr(bp, time);

   A.set_row(replist.size());
   Algebra::CK_decomp(A, L, cm);
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

// update current values for all LDOs
void SubCircuit::update_ldo_current(){
	Node *nd;
	Net *net;
	Node *nb;
	double current=0;
	for(size_t i=0;i<ldolist.size();i++){
		nd = ldolist[i]->A;
		for(DIRECTION d = WEST; d!= BOTTOM; d = (DIRECTION)(d+1)){
			net = nd->nbr[d];
			if(net == NULL) continue;
			// no current net in ldo 
			// locations in input file
			if(net->type == CURRENT) 
				continue;
			// work on resistance net
			if(net->ab[0] == nd)
				nb = net->ab[1];
			else
				nb = net->ab[0];
			current += (nd->value - nb->value ) / net->value;
			// clog<<"net, nd, nb, current: "<<*net<<" "<<*nd<<" "<<*nb<<current;
		}
		// copy old current
		ldolist[i]->current_old = 
			ldolist[i]->current;
		// assign new current
		ldolist[i]->current = current;
		// clog<<"ldo current: "<<current;
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
		rid = nd->rep->rid;
		// clog<<"old bp: "<<bp[rid]<<endl;
		// restore old current
		bp[rid] += ldolist[i]->current_old;
		// change into new one
		bp[rid] -= ldolist[i]->current;
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
	double min_IR = max_IRdrop;	
	// clog<<"min_IR initial is: "<<min_IR<<endl;
	// for local pad movement
	int pad_number = 1;
	// find control nodes for each pad
	extract_pads(pad_number);
	// find the tune spot for control nodes	
	update_pad_control_nodes(ref_drop_vec, 0);	

	// find new point for all pads	
	update_pad_pos_all(ref_drop_vec);

	// move pads according to graph contraints
	graph_move_pads(ref_drop_vec, true);

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

void SubCircuit::build_pad_graph(){
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
void SubCircuit::extract_pads(int pad_number){
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
void SubCircuit::update_pad_control_nodes(vector<double> & ref_drop_value, size_t iter){
	ref_drop_value.resize(pad_set.size());
	for(size_t i=0;i<pad_set.size();i++){
		if(pad_set[i]->control_nodes.size()==0)
			continue;
		
		double middle_value = locate_ref(i);
		ref_drop_value[i] = middle_value;
		//cout<<"middle value: "<<middle_value<<endl;
	}
}

bool compare_values(double a, double b){
	return (a<b);
}

// locate the tune spot for the control nodes.
double SubCircuit::locate_ref(size_t i){
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
		// need to be generated with worst_cur
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

void SubCircuit::dynamic_update_violate_ref(double VDD, vector<double> & ref_drop_vec, bool local_flag){
	//for(size_t j=0;j<2;j++){
	double avg_drop = 0;
	avg_drop = calc_avg_ref_drop(ref_drop_vec);

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
	
	extract_min_max_pads_new(VDD, ref_drop_vec, local_flag);	
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

void SubCircuit::extract_min_max_pads_new(double VDD, vector<double> ref_drop_vec, bool local_flag){
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
	double avg_ref = calc_avg_ref(pad_set, ref_drop_vec);
	double avg_drop = VDD - avg_ref;
	
	for(size_t i=0;i<pad_set.size();i++){
		pad = pad_set[i];
		nd = pad->node;
					
		drop = VDD_G - ref_drop_vec[i];
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

				new_pad = pad_projection(pad_ptr, min_pads[j], local_flag);
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

// expand from (x,y) to nearest node in grid
// fits for non-uniform grid
Node * SubCircuit::pad_projection( 
	Pad *pad, Node *nd,
	bool local_flag){
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
	// clog<<"pt_name: "<<name<<endl;
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

	// then start to expand to find candi LDO nodes
	nd_new_ldo = expand_ldo_location(ref_dist, 
			ref_x, ref_y, ldo_ptr);

	ldo->A = nd_new_ldo;
	// also need to change nd_in of LDO in ckt_g
	return nd_new_ldo; 
}

double SubCircuit::update_pad_pos_all(vector<double> ref_drop_vec){
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

// use queue to search for candi pads around nd
Node * SubCircuit::expand_candi_pads(Node *nd){
	Node *na;
	Node *nd_new;
	queue<Point> q;
	Point pt;
	Point pt_cur;
	stringstream sstream;
	string pt_name;
	bool return_flag = false;
	q.push(nd->pt);
	
	double dx[4] = {3, 0, -3, 0};
	double dy[4] = {0, 3, 0, -3};
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
			sstream <<"n"<<pt_nbr.z<<"_"<<
				pt_nbr.x<<"_"<<
				pt_nbr.y;
			pt_name = sstream.str();
			if(has_node(pt_name)){
				nd_new = get_node(pt_name);
				Net *net = nd_new->nbr[TOP];
				if(net == NULL){
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
	if(return_flag == true){
		return nd_new;
	}
	clog<<"no point for new pad. return. "<<endl;
	return NULL;
}

void SubCircuit::move_violate_pads(vector<double> ref_drop_vec, bool local_flag){
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
			new_pad = pad_projection(pad_ptr, pad, local_flag);
			// clog<<"old pad / new pad: "<<*pad<<" "<<*new_pad<<endl;
			pad_ptr->node = new_pad;
			pad_ptr->control_nodes.clear();
			pad_ptr->visit_flag = true;
		}
	}
}

void SubCircuit::graph_move_pads(vector<double> ref_drop_vec, bool local_flag){
	Node *new_pad;
	int id=0;
	// for(size_t i=0;i<5;i++){
	do{
		id = locate_max_drop_pad(ref_drop_vec);
		if(id==-1) break;
		Pad *pad_ptr = pad_set[id];
		Pad *pad_nbr = NULL;
		Node *pad = pad_ptr->node;
		// clog<<"old_pad: "<<*pad<<endl;
		new_pad = pad_projection(pad_ptr, pad, local_flag);
		// clog<<" old pad / new pad: "<<*pad<<" "<<*new_pad<<endl;

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
	}while(id != -1);
}

// locate id that has minimum value and not visited or fixed 
int SubCircuit::locate_max_drop_pad(vector<double> vec){
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

void SubCircuit::clear_flags(){
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
	Node *nd = ldo_ptr.A;
	Node *nd_new;
	int width = ldo_ptr.width;
	int height = ldo_ptr.height;
	
	int dx[4]; int dy[4];
	dx[0] = width; dx[2] = -width;
	dx[1] = dx[3] = 0;
	dy[1] = height; dy[3] = -height;
	dy[0] = dy[2] = 0;

	queue<Point> q;
	Point pt;
	pt.x = ref_x;
	pt.y = ref_y;

	Point pt_cur;
	stringstream sstream;
	string pt_name;
	int return_flag = 0;;
	// else start to search for node
	q.push(pt);
	int diff_x, diff_y;
	double dist;
	// if not, expand it to neighboring area
	while(!q.empty()&& return_flag == false){
		pt_cur = q.front();
		// clog<<"pt_cur: "<<pt_cur<<endl;
		Point pt_nbr = pt_cur;
		//expand_pad_pos(q, pt_cur);	
		for(size_t i=0;i<4;i++){
			pt_nbr.x = pt_cur.x + dx[i];
			pt_nbr.y = pt_cur.y + dy[i];
			stringstream sstream;
			string name;
			sstream <<"n"<<nd->pt.z<<"_"<<
				pt_nbr.x<<"_"<<
				pt_nbr.y;
			name = sstream.str();
			nd_new = get_node(name);
			if(nd_new == NULL) continue;
			nd_new = nd_new->rep;
			if(nd_new->isS()==Y) continue;
			// clog<<"nd_new: "<<*nd_new<<" "<<nd_new->get_geo_flag()<<endl;
			// get candidate location
			if(nd_new->get_geo_flag()==SBLANK){
				diff_x = nd_new->pt.x - ref_x;
				diff_y = nd_new->pt.y - ref_y;
				dist = sqrt(diff_x*diff_x + diff_y*diff_y);
				if(dist < ref_dist){
					//clog<<"dist, ref_dist: "<<dist<<" "<<ref_dist<<endl;
					// clog<<"new name: "<<*nd_new<<" "<<endl;
					return_flag = 1;
					break;
				}else{// already out of range, break
					return_flag = 2;
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
	if(return_flag == 1){
		return nd_new;
	}
	// no candidate, return original node
	return nd;
}

// perform LDO node change
// modify nets and nodes with the change of LDO
void SubCircuit::rebuild_local_nets(Node *rm_node, Node *add_node){
	//Node *rm_node=NULL;
	//Node *add_node=NULL;

	/*for(size_t i=0;i<origin_pad_set.size();i++){
		clog<<"old/new: "<<*origin_pad_set[i]<<" "<<*pad_set[i]->node<<endl;
	}*/

	//for(size_t i=0;i<origin_pad_set.size();i++){
		// one-one correspondence
		//rm_node = origin_pad_set[i]->rep;	
		//add_node = pad_set[i]->node->rep;

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
		//clog<<"net: "<<*net<<endl;
	//}
}
	
// for global circuit
// rebuild all the RLV nets from ldolist
void SubCircuit::rebuild_global_nets(){
	pad_set.clear();
	Node *nd;
	static Node nd_temp;
	Net *net;
	stringstream sstream;

	double vol_value = 0;
	double induc_value = 0;
	double resis_value = 0;
	double current_value = 0;

	char name[MAX_BUF];

	for(int type =0 ;type <NUM_NET_TYPE;type++){
		NetList &ns = net_set[type];
		if(ns.size()<=0) continue;
		net = ns[0];
		// if(net != NULL)
			// clog<<"net: "<<*net<<endl;
		if(net->type == VOLTAGE)
			vol_value = net->value;
		else if(net->type == INDUCTANCE)
			induc_value = net->value;
		else if(net->type == RESISTOR)
			resis_value = net->value;
		else if(net->type == CURRENT)
			current_value = net->value;
	}

	// then delete all the nets in global circuit
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		for(size_t j=0;j<ns.size();j++) delete ns[j];
		ns.clear();
	}
	Node *gnd = nodelist[nodelist.size()-1];
	// keep the old nd_in node
	int old_z = ldolist[0]->nd_in->pt.z;
	// and delete all the nodes except ground 
	for(size_t i=0;i<nodelist.size()-1;i++) 
		delete nodelist[i];
	replist.clear();
	nodelist.clear();
	map_node.clear();
	nodelist.push_back(gnd);

	NET_TYPE net_type;
	// start to create new nets
	for(size_t i=0;i<ldolist.size();i++){
		nd = ldolist[i]->A;
		sprintf(name, "_Y_n%ld_%ld_%ld", ldolist[i]->nd_in->pt.z, nd->pt.x, nd->pt.y);
		extract_node(name, nd_temp);
		Node *nd_ptr_vol = new Node(nd_temp);
		nd_ptr_vol->rep = nd_ptr_vol;
		add_node(nd_ptr_vol);
		// first create voltage net
		net_type = VOLTAGE;
		Net *net_vol = new Net(net_type, vol_value, nd_ptr_vol, gnd);
		add_net(net_vol);
		update_node(net_vol);
		// then create inductance net
		sprintf(name, "_X_n%ld_%ld_%ld", ldolist[i]->nd_in->pt.z, nd->pt.x, nd->pt.y);
		extract_node(name, nd_temp);
		Node *nd_ptr_induc = new Node(nd_temp);
		nd_ptr_induc->rep = nd_ptr_induc;
		add_node(nd_ptr_induc);
		net_type = INDUCTANCE;
		Net *net_induc = new Net(net_type, induc_value, nd_ptr_induc, nd_ptr_vol);
		add_net(net_induc);
		update_node(net_induc);
		// then create resistance net 
		sprintf(name, "n%ld_%ld_%ld", ldolist[i]->nd_in->pt.z, nd->pt.x, nd->pt.y);
		extract_node(name, nd_temp);
		Node *nd_ptr_resis = new Node(nd_temp);
		nd_ptr_resis->rep = nd_ptr_resis;
		add_node(nd_ptr_resis);
		net_type = RESISTOR;
		Net *net_resis = new Net(net_type, 
			resis_value, nd_ptr_resis, nd_ptr_induc);
		add_net(net_resis);
		update_node(net_resis);
		// then create current net 
		net_type = CURRENT;
		Net *net_current = new Net(net_type, 
			current_value, nd_ptr_resis, nodelist[0]);
		add_net(net_current);
		update_node(net_current);
	}
	// then update ldo->nd_in;
	for(size_t i=0;i<ldolist.size();i++){
		Node *nd = ldolist[i]->A;
		stringstream sstream;
		sstream<<"n"<<old_z<<"_"<<nd->pt.x<<"_"<<nd->pt.y;
		Node *nd_new = get_node(sstream.str());
		ldolist[i]->nd_in = nd_new;
	}
	/*clog<<"to here. "<<endl;
	clog<<nodelist<<endl;
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		for(size_t j=0;j<ns.size();j++)
			clog<<*ns[j]<<endl;
	}*/
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

// build candidate padlist from wspace
void SubCircuit::build_candi_graph(){
	
}

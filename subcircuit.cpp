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
	Node *nd;
	for(size_t i=0;i<ldolist.size();i++){
		/*LDO *ldo_ptr = ldolist[i];
		//clog<<"ldo node: "<<*ldo_ptr->A<<endl;
		Net *net = ldo_ptr->A->nbr[TOP];
		if(net == NULL) {
			//clog<<"null net. "<<endl;
			continue;
		}
		nd = net->ab[0];
		if(nd->name == ldo_ptr->A->name)
			nd = net->ab[1];*/
		nd = ldolist[i]->A;
		nd->enableLDO();
		//clog<<"nd is LDO: "<<*nd<<endl;
	}
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
   cm = &c;
   cholmod_start(cm);
   cm->print = 5;
   b = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
   x = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
   bp = static_cast<double *> (b->x);

   //Matrix A;
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
   //if(replist.size()<THRESHOLD){
   	for(i=0;i<replist.size();i++)
		bp[i] = 0;
   //}
#if 0
   else{
#pragma omp parallel for private(i)
	for(i=0;i<replist.size();i++)
		bp[i] = 0;
   }
#endif
   //link_tr_nodes(tran);
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
   //if(n<THRESHOLD){
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
   //}
#if 0
   else{
	   size_t i=0;
#pragma omp parallel for private(i)
	   for(i=0;i<n;i++){
		   int id = id_map[i];
		   replist[id]->rid = i;
		   temp[i] = bp[i];
	   }
#pragma omp parallel for private(i)
	   for(i=0;i<n;i++)
		   bp[i] = temp[id_map[i]];
#pragma omp parallel for private(i)
	   for(i=0;i<n;i++)
		   temp[i] = xp[i];
#pragma omp parallel for private(i)
	   for(i=0;i<n;i++)
		   xp[i] = temp[id_map[i]];
   }
#endif
   delete [] temp;
   delete [] id_map;
   /*****************************************/ 
//#endif
   // bnewp[i] = bp[i]
   //if(n<THRESHOLD){
	for(size_t i=0;i<n;i++)
		bnewp[i] = bp[i];
   //}
#if 0
   else{
	size_t i=0;
#pragma omp parallel for private(i)
	for(i=0;i<n;i++)
		bnewp[i] = bp[i];
   }
#endif
   //stamp_current_tr(bnewp, tran, time);
   
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
	// bnewp[i] = bp[i];
	//if(n<THRESHOLD){
		for(size_t i=0;i<n;i++)
			bnewp[i] = bp[i];
	//}
#if 0
	else{
		size_t i=0;
#pragma omp parallel for private(i)
		for(i=0;i<n;i++)
			bnewp[i] = bp[i];
	}
#endif
      // only stamps if net current changes
      // set bp into last state
      //stamp_current_tr(bnewp, tran, time);
      stamp_current_tr_1(bp, bnewp, time);
     // get the new bnewp
      modify_rhs_tr(bnewp, xp); 
	
      solve_eq_sp(xp);

      //save_tr_nodes(tran, xp);
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
void SubCircuit::get_voltages_from_LU_sol(double * x){
   size_t i;
   for(i=0;i<nodelist.size()-1;i++){
      Node * node = nodelist[i];
      size_t id = node->rep->rid;  // get rep's id in Vec
      double v = x[id];		// get its rep's value
      node->value = v;

   }
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

           //cout<<"("<<k<<" "<<k<<" "<<G<<")"<<endl;
           if(!nl->is_ground() && nl->isS()!=Y
			   &&(nl->nbr[TOP]==NULL || 
                 nl->nbr[TOP]->type != INDUCTANCE)&&(k > l)){
                    A.push_back(k,l,-G);

                  //cout<<"("<<k<<" "<<l<<" "<<-G<<")"<<endl;
           }
        }

	if( !nl->is_ground() && nl->isS() !=Y && 
			(nl->nbr[TOP] ==NULL 
		||nl->nbr[TOP]->type != INDUCTANCE)) {
		A.push_back(l,l, G);
                //cout<<"("<<l<<" "<<l<<" "<<G<<")"<<endl;
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
}

// build global current nets from LDO
// only one time step
void SubCircuit::build_global_nets(){
	Node *nd;
	double current;
	for(size_t i=0;i<ldolist.size();i++){
		nd = ldolist[i]->nd_in;
		current = ldolist[i]->current;
		// current = 0.3;
		Net *net = new Net(CURRENT, current, nd, nodelist[nodelist.size()-1]);
		add_net(net);
		// update top nbr net
		nd->rep->nbr[BOTTOM] = net;
		//clog<<"add net: "<<*net<<endl;
	}
}
// ----------------------------------------------------------------//
// Filename : circuit.cpp
// Author : Zigang Xiao <zxiao2@illinois.edu>
//          Ting Yu <tingyu1@illinois.edu>
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
#include <sstream>
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
	for(size_t i=0;i<nodelist.size();i++) delete nodelist[i];
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		for(size_t j=0;j<ns.size();j++) delete ns[j];
	}
	//delete [] id_map;
	Lx = NULL;
	Li = NULL;
	Lp = NULL;
	Lnz = NULL;
	table_ldo.clear();
	ldo_vin_vec.clear();
	ldo_iout_vec.clear();
}

void Circuit::check_sys() const{
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
# if 0
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
#endif
string Circuit::get_name() const{return this->name;}

#if 0
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
#endif
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
void Circuit::solve_init(){	
	ckt_g.name = "GLOBAL";
	ckt_l.name = "LOCAL";
	ckt_l.lx = lx; ckt_l.gx = gx; 
	ckt_l.ly = ly; ckt_l.gy = gy;

	ckt_g.lx = lx; ckt_g.gx = gx; 
	ckt_g.ly = ly; ckt_g.gy = gy;
	ckt_l.layers = local_layers;
	ckt_g.layers = global_layers;
	ckt_l.assign_min_max_layers();
	ckt_g.assign_min_max_layers();
	// assign nodes and nets into ckt_g and ckt_l
	build_subcircuit();
	clog<<"after build subcircuit. "<<endl;
	// update and sort nodes
	ckt_l.solve_init(true);
	ckt_g.solve_init(false);
	// configure subckt, start cholmod process
	ckt_l.configure_init();
	ckt_g.configure_init();
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
			else if((p->nbr[NORTH] !=NULL && p->nbr[SOUTH] !=NULL)			      &&(p->nbr[EAST]==NULL && p->nbr[WEST] ==NULL))
				count++;
		}
	}
	//clog<<"number of nodes can be merged is: "<<count<<endl;
}
//if 0
void Circuit::solve(Tran &tran){
	//1. assign nodes and nets to ckt_g and l
	//2. config subckt, start cholmod for ckt_g and l
	solve_init();
	// build pad set: local build LDO list - need update
	ckt_l.build_local_pad_set();
	// global pad set include LDO nets
	ckt_g.build_global_pad_set();

	// build new nets for the single LDo
	ckt_l.build_local_nets();
	ckt_g.build_global_nets();
// #if 0	
	// solving LDO location with DC
	bool flag = ckt_l.optimize_single_ldo();
	// not in the limit, need to op
	if(flag == true)
		ckt_l.add_ldo_DC(tran);
	clog<<"after local op DC LDO. "<<endl;
	// clear flag_visited for the pads
	// ckt_l.clear_flags();
	// clog<<ckt_l.nodelist<<endl;
	int iter = 0;
	clog<<endl;
	flag = false;

	return;
#if 0
	ckt_l.add_ldo_tr(Tran &tran);	
	// only need to stamp matrix once per t step
	ckt_l.stamp_decomp_matrix_TR(tran);
	ckt_g.stamp_decomp_matrix_TR(tran);
	clog<<"after ckt_g stamp tr. "<<endl;

	// double time = 0;
	// add_LDO_all_global(tran, time);
	// go along all time steps
	for(double time =0; time < tran.tot_t;// && iter <1; 
			time += tran.step_t){
		clog<<"===== "<<time<<" ===="<<endl;
		// solve one time step with LDO

		flag = solve_TR(tran, time);
		clog<<"after solve_TR. "<<endl;
		// flag = solve_TR_LDO(tran, time);
		// clog<<"maxIR_l / g: "<<ckt_l.locate_maxIRdrop()<<" "<<ckt_g.locate_g_maxIRdrop()<<endl;
		// ckt_l.clear_flags();
		iter++;
		if(flag == true){
			clog<<"reach maximum LDO size: "<<endl;
			break;
		}
	}
	/*clog<<"print matlab. "<<endl;
	ckt_l.print_matlab_node();
	ckt_l.print_matlab_LDO();
	*/
	// vector<LDO*>::iterator it;
	// it = unique(ckt_l.ldolist.begin(), ckt_l.ldolist.end());
	// ckt_l.ldolist.resize(std::distance(ckt_l.ldolist.begin(), it));
	clog<<"final ldo size: "<<ckt_l.ldolist.size()<<" "<<ckt_l.MAX_NUM_LDO<<endl;
	// for(size_t i=0;i<ckt_l.ldolist.size();i++)
		// cout<<"i, ldo: "<<i<<" "<<*ckt_l.ldolist[i]->A<<endl;
	// return;
	// clog<<endl;
	clog<<"==== entering verify stage ==== "<<endl;
// #endif
	verify_solve(tran);

#endif
	// release resouces
	ckt_l.release_resources();
	ckt_g.release_resources();
}
//endif

// build nodelist and netlist in subcircuit
void Circuit::build_subcircuit(){
	ckt_l.ldolist = ldolist;
	ckt_g.ldolist = ldolist;
	// build nodelist
	for(size_t i=0;i<nodelist.size();i++){
		if(nodelist[i]->is_ground())
			continue;
		int layer = nodelist[i]->get_layer();
		bool flag = false;
		for(size_t j=0;j<local_layers.size();j++){
			if(layer == local_layers[j]){
				flag = true;
				break;
			}
		}
		if(flag == true){
			ckt_l.nodelist.push_back(nodelist[i]);
		}
		else{
			ckt_g.nodelist.push_back(nodelist[i]);
		}
	}
	// then build net_set
	Node *nd;
	Net *net;
	// handles the LDO_NET independently
	for(int type =0;type<NUM_NET_TYPE-1;type++){
		NetPtrVector &ns = net_set[type];
		for(size_t i=0;i<ns.size();i++){
			net = ns[i];
			if(net == NULL) continue;
			nd = net->ab[0];
			if(nd->is_ground())
				nd = net->ab[1];
			int layer = nd->get_layer();
			bool flag = false;
			for(size_t j=0;j<local_layers.size();j++){
				if(layer == local_layers[j]){
					flag = true;
					break;
				}
			}
			if(flag == true){
				// clog<<"local net: "<<*net<<endl;
				ckt_l.net_set[type].push_back(net);
			}
			else{
				// clog<<"global net: "<<*net<<endl;
				ckt_g.net_set[type].push_back(net);
			}
		}
	}

	// write the LDO_net to both l and g
	int type = NUM_NET_TYPE-1;
	NetPtrVector &ns = net_set[type];
	for(size_t i=0;i<ns.size();i++){
		net = ns[i];
		if(net == NULL) continue;
		ckt_l.net_set[type].push_back(net);
		ckt_g.net_set[type].push_back(net);
		
	}
	// build the map between local and global node
	ckt_l.build_map_landg();
	ckt_g.build_map_landg();

	// mark the geometry with occupy info
	ckt_l.mark_geo_occupation();
	clog<<"MAX_NUM_LDO: "<<ckt_l.MAX_NUM_LDO<<endl;
	// ckt_g.build_candi_pad_set();
}

// verify final solution with optimized LDOs
void Circuit::verify_solve(Tran &tran){
   // ckt_l.build_pad_set();
   // ckt_g.build_pad_set();

   ckt_l.reset_b();
   ckt_g.reset_b();
   clog<<ckt_l.ldolist.size()<<" "<<ckt_g.pad_set.size()<<" "<<ckt_l.MAX_NUM_LDO<<endl; 
   // first solve DC solution
   solve_DC();
   clog<<"after verify DC. "<<endl;
   // clog<<"max_IR for l and g: "<<ckt_l.locate_maxIRdrop()<<" "<<ckt_g.locate_g_maxIRdrop()<<endl;
   // clog<<ckt_l.nodelist<<endl;
   // clog<<ckt_g.nodelist<<endl;
   // return;
   // link ckt nodes for local circuit
   ckt_l.link_ckt_nodes(tran);
   
   double time = 0;
   // stamp A and bp, and decompose
   ckt_l.stamp_decomp_matrix_TR(tran);
   ckt_g.stamp_decomp_matrix_TR(tran);
   // clog<<"after stamp matrix. "<<endl;

   double min_IRdrop = 0;
   double min_step = 0;
   bool flag = false;
   int iter =0;
   for(time = 0; time < tran.tot_t;
		   time += tran.step_t){
	   // cout<<endl<<"===== new step: "<<time <<" ====="<<endl;
   	verify_one_LDO_step(tran, time);
	locate_maxIRdrop();
	if(flag == false){
		min_IRdrop = max_IRdrop;
		min_step = time;
		flag = true;
	}else{
		if(max_IRdrop > min_IRdrop){
			min_IRdrop = max_IRdrop;
			min_step = time;
		}
	}
	iter ++;
	// clog<<"max IR after solve_TR is: "<<max_IRdrop<<endl;
   }
   for(time = 0; time < tran.tot_t;
		   time += tran.step_t){
   	if(time - min_step < 1e-11){
		clog<<"print matlab. "<<endl;
		ckt_l.print_matlab_node();
		ckt_l.print_matlab_LDO();
		break;
	}
   }
   clog<<"min_IRdrop, step: "<<min_IRdrop<<" "<<min_step<<endl;
   ckt_l.save_ckt_nodes_to_tr(tran);
   clog<<"after verify solve." <<endl;
}

#if 0
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
   //for(size_t i=0;i<n;i++)
	//cout<<"dc solution from ckt: "<<
	//get_name()<<" "<<xp[i]<<endl;
   
   // print out dc solution
   //cout<<"dc solution. "<<endl;
   //cout<<nodelist<<endl;

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
   //clog<<"before tr solve."<<endl;
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
#endif
// change the wspace into a local variable
// here is simply copy
void Circuit::build_wspacelist(vector<MODULE*> wspace_vec){
	wspacelist.clear();
	wspacelist = wspace_vec;
	ckt_l.wspacelist = wspace_vec;
	ckt_g.wspacelist = wspace_vec;
}

void Circuit::build_ldolist(vector<LDO*> ldo_vec){
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
void Circuit::solve_LU(Tran &tran){
        solve_init();
	solve_LU_core(tran);
}
#endif

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
		case LDO_NET:
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
		case LDO_NET:
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
        if( !nk->is_ground()&& nk->isS()!=Y && 
          (nk->nbr[TOP]== NULL 
           || nk->nbr[TOP]->type != INDUCTANCE)) {
           A.push_back(k,k, G);

           //clog<<"("<<k<<" "<<k<<" "<<G<<")"<<endl;
           if(!nl->is_ground() &&(nl->nbr[TOP]==NULL || 
                 nl->nbr[TOP]->type != INDUCTANCE)&&(k > l)){
                    A.push_back(k,l,-G);

                  //clog<<"("<<k<<" "<<l<<" "<<-G<<")"<<endl;
           }
        }

	if( !nl->is_ground() && nl->isS() !=Y && 
			(nl->nbr[TOP] ==NULL 
		||nl->nbr[TOP]->type != INDUCTANCE)) {
		A.push_back(l,l, G);
                //clog<<"("<<l<<" "<<l<<" "<<G<<")"<<endl;
		if(!nk->is_ground()&& (nk->nbr[TOP]==NULL ||
		  nk->nbr[TOP]->type != INDUCTANCE) && l > k){
			A.push_back(l,k,-G);
		//clog<<"("<<l<<" "<<k<<" "<<-G<<")"<<endl;
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

   if( nk->isS()!=Y && !nk->is_ground()&& 
     (nk->nbr[TOP]!=NULL && 
      nk->nbr[TOP]->type == INDUCTANCE)) {
      //clog<<"net: "<<*net<<endl;
      A.push_back(k,k, G);
      //clog<<"("<<k<<" "<<k<<" "<<G<<")"<<endl;
      if(!nl->is_ground() &&(nl->nbr[TOP]==NULL || 
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
      if(!nk->is_ground()&& (nk->nbr[TOP]==NULL ||
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

void Circuit::stamp_inductance_dc(Matrix & A, double *b, Net * net){
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
void Circuit::stamp_capacitance_dc(Matrix & A, Net * net){
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
void Circuit::stamp_capacitance_tr(Matrix &A, Net *net, Tran &tran){
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
void Circuit::modify_rhs_c_tr_0(Net *net, double * rhs, double *x){
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
		 rhs[k] += Ieq;	// for VDD circuit
		//clog<<*nk<<" rhs +: "<<rhs[k]<<endl;
	}
	if(!nl->is_ground()&& nl->isS()!=Y){
		 rhs[l] += -Ieq; 
		//clog<<*nl<<" rhs +: "<<rhs[l]<<endl;
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
	if(!nk->is_ground()&& nk->isS()!=Y){
		 rhs[k] += Ieq;	// for VDD circuit
	}
	if(!nl->is_ground()&& nl->isS()!=Y){
		 rhs[l] -= Ieq; 
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
        pg.node_set_x.push_back(k);
        pg.node_set_x.push_back(id_b);
//#endif
	Ieq  = i_t + temp;
	//clog<<"Ieq: "<<Ieq<<endl;
	if(nk->isS() !=Y && !nk->is_ground()){
		 rhs[k] += Ieq; // VDD circuit
		//clog<<*nk<<" "<<rhs[k]<<endl;
	}
	if(nl->isS()!=Y && !nl->is_ground()){
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
	if(nk->isS() !=Y && !nk->is_ground()){
		 rhs[k] += Ieq; // VDD circuit
		//clog<<*nk<<" "<<rhs[k]<<endl;
	}
	if(nl->isS()!=Y && !nl->is_ground()){
		 rhs[l] += -Ieq; // VDD circuit
		//clog<<*nl<<" "<<rhs[l]<<endl;
	}
}
// stamp a current source
void Circuit::stamp_current(double * b, Net * net){
	//clog<<"net: "<<*net<<endl;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;

	if( !nk->is_ground() && nk->isS()!=Y){// && 
		size_t k = nk->rid;
		b[k] += -net->value;
		//clog<<"b: "<<k<<" "<<-net->value<<endl;
	}
	if( !nl->is_ground() && nl->isS() !=Y){// &&
		size_t l = nl->rid;
		b[l] +=  net->value;
		//clog<<"b: "<<l<<" "<<-net->value<<endl;
	}
}

void Circuit::stamp_current_tr_net(double * b, Net * net, double &time){
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
	 //cout<<"nodes: "<<ckt_nodes[j].node->name<<endl;
         id = ckt_nodes[j].node->rep->rid;
	 //cout<<"value: "<<x[id]<<endl;
         ckt_nodes[j].value.push_back(x[id]);
      }
}

// link transient nodes with nodelist
void Circuit:: link_ckt_nodes(Tran &tran){
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

void Circuit::make_A_symmetric_tr(double *b, double *x, Tran &tran){
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

#if 0
bool compare_Node_G(const Node_G *nd_1, const Node_G *nd_2){
   return (nd_1->value < nd_2->value);
 }
#endif
# if 0
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
#endif
#if 0
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
#endif
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

void Circuit::build_pad_set_l(vector<LDO*> ldolist){
	pad_set_l.clear();
	//Node *nd;
	for(size_t i=0;i<ldolist.size();i++){
		Pad *pad_ptr  = new Pad();
		pad_ptr->node = ldolist[i]->A;
		pad_set_l.push_back(pad_ptr);
	}
}

void Circuit::build_pad_set(){
	pad_set_g.resize(0);
	pad_set_l.resize(0);
	//static Pad pad_a;
	for(size_t i=0;i<nodelist.size()-1;i++){
		if(nodelist[i]->isS()==X){
			Pad *pad_ptr = new Pad();
			pad_ptr->node = nodelist[i];
			if(nodelist[i]->is_LDO() == false){
			//if(nodelist[i]->get_layer()== global_layers[0])
				//cout<<"global pad: "<<*nodelist[i]<<endl;
				pad_set_g.push_back(pad_ptr);
			}
			else{
				//cout<<"local pad: "<<*nodelist[i]<<endl;
				pad_set_l.push_back(pad_ptr);
			}
		}
	}
}

// solve_DC with LDO locations
void Circuit::solve_DC(){
	int iter=0;
	double diff_l = 1;
	double diff_g = 1;
	
	// first atamp matrix: already symmetric
	cout<<endl<<"<==== local ckt stamping: ====>"<<endl;
	ckt_l.stamp_decomp_matrix_DC();
	cout<<endl<<"<==== global ckt stamping: ====>"<<endl;
	ckt_g.stamp_decomp_matrix_DC();

	
		
	// solve DC with fixed number of LDOs
	while((diff_l > 1e-4 || diff_g > 1e-4) && iter <1){
		ckt_l.reset_b();
		ckt_l.stamp_rhs_DC(true);	
		diff_l = ckt_l.solve_CK_with_decomp();
		
		//clog<<ckt_l.nodelist<<endl;
		//clog<<"local max_IR: "<<ckt_l.locate_maxIRdrop()<<endl;
		ckt_g.reset_b();
		ckt_g.stamp_rhs_DC(false);	
		diff_g = ckt_g.solve_CK_with_decomp();
		// clog<<ckt_g.nodelist<<endl;
		//clog<<"global max_IR: "<<ckt_g.locate_maxIRdrop()<<endl;
		// then generate the ldo info for SPICE
		extract_ldo_info();
		
		// clog<<"iter, diff_l, diff_g: "<<iter<<" "<<diff_l<<" "<<diff_g<<endl<<endl;
		iter++;
	}
	// locate_maxIRdrop();		
}

// solve one transient step with different LDO numbers
bool Circuit::solve_TR_LDO(Tran &tran, double time){
	// clog<<ckt_l.nodelist<<endl;
	// 1. first modify L and C nets
	// 2. no need to build new nets, but modify
	bool max_flag = solve_TR(tran, time);
	return max_flag;
}

// solve circuit with fixed LDO for one time step
bool Circuit::solve_TR(Tran &tran, double time){
	int iter=0;
	double diff_l = 1;
	double diff_g = 1;

	// only need to stamp matrix once per t step
	// ckt_l.stamp_decomp_matrix_TR(tran);
	// ckt_g.stamp_decomp_matrix_TR(tran);
	// stores the Ieq for both C and L nets 
	ckt_l.modify_rhs_tr_0(tran);
	ckt_g.modify_rhs_tr_0(tran);
	int local_bad_flag = 0;
	double THRES_l = VDD_G*0.1;
	// double THRES_g = 2.2*0.1;
	bool max_flag = false;
	double diff_old_l = 1;
	double diff_old_g = 1;
	double delta_diff_l = 1;
	double delta_diff_g = 0;
	bool flag = true;

	// int count_log = 0;
	ckt_l.calculate_tr_current(time);
	while((delta_diff_l > 1e-3 || delta_diff_g > 1e-3)|| flag == true){//&& iter < 20){
		// if(count_log++ >=5) break;
		if(delta_diff_l <=1e-3 && delta_diff_g <= 1e-3 && ckt_l.max_IRdrop <= THRES_l)
			flag = false;
		// clog<<endl<<"time - iter: "<<time<<" "<<iter<<endl;
		// then update netlist
		ckt_l.modify_local_nets();
		ckt_l.stamp_rhs_tr(true, time, tran);
		// solve eq with decomped matrix
		diff_old_l = diff_l;
		diff_l = ckt_l.solve_CK_with_decomp_tr();
		delta_diff_l = fabs(diff_l - diff_old_l);
		ckt_l.locate_maxIRdrop();
		
		if(iter ==0)
			clog<<"ckt_l max_IR: "<<ckt_l.max_IRdrop<<" "<<THRES_l<<endl;
		if(ckt_l.max_IRdrop >= THRES_l){
			/*if(fabs(time - 6.2e-9) < 1e-15 && iter ==0){
				clog<<"print matlab. "<<endl;
				ckt_l.print_matlab_node();
				ckt_l.print_matlab_LDO();
			}*/
			// clog<<"time: "<<time - 9.3e-9<<endl;
			
			// clog<<"start to optimize ckt_l max_IR: "<<ckt_l.max_IRdrop<<endl;
			local_bad_flag = 1;
			max_flag = optimize_local_LDO(local_bad_flag, THRES_l, tran, time);
			iter = 0;
			// cannot improve anymore
			/*if(fabs(op_local_IR - temp_local_IR) < 1e-4){
				flag = false;
				clog<<"done with optimization. "<<endl;
			}*/
			clog<<"optimized new local max_IR: "<<ckt_l.max_IRdrop<<" "<<ckt_l.ldolist.size()<<" "<<flag<<endl;
		}
		// calculate ldo current from ckt_l
		ckt_l.update_ldo_current();
		// clog<<"after update current. "<<endl;
		// break;	
//#if 0
		// restamp global rhs with ldo current
		// ckt_g.modify_ldo_rhs();	
		ckt_g.stamp_rhs_tr(false, time, tran);
		diff_old_g = diff_g;
		diff_g = ckt_g.solve_CK_with_decomp_tr();
		delta_diff_g = fabs(diff_g - diff_old_g);
		// clog<<"iter, delta_diff_l / g; "<<iter<<" "<<delta_diff_l<<" "<<delta_diff_g<<endl;
		ckt_g.locate_g_maxIRdrop();
		// clog<<"global max_IR: "<<ckt_g.max_IRdrop<<endl;
		//clog<<ckt_g.nodelist<<endl;
		//update_ldo_vout();
				
// #endif
		iter++;
	}
	/*if(fabs(time - 6.2e-9) < 1e-15){
		clog<<"print matlab. "<<endl;
		ckt_l.print_matlab_node();
		ckt_l.print_matlab_LDO();
	}*/
	return max_flag;
	// clog<<"final ldo size: "<<ckt_l.ldolist.size()<<endl;
	// cout<<"final local: "<< ckt_l.nodelist<<endl;
	// cout<<"final global: "<<ckt_g.nodelist<<endl;
}

// readin LDO lookup table
void Circuit::Readin_LDO(){
	FILE *f;
	key_LDO key;
	pair<key_LDO, double> LDO_pair;
	ldo_vin_vec.clear();
	ldo_iout_vec.clear();
	f = fopen("../data/LDO/LDO_lookupTable.txt", "r");
	if( f == NULL ) 
		report_exit("LDO table file not exist!\n");
	char line[MAX_BUF];
	double vin, iout, vout;
	vector<double>::iterator it;
	while(fgets(line, MAX_BUF, f) != NULL){
		if(line[0]== '*') continue;
		sscanf(line, "%lf %lf %lf", &vin, &iout, &vout);
		key.first = vin;
		key.second = iout;
		LDO_pair.first = key;
		LDO_pair.second = vout;
		/* clog<<LDO_pair.first.first<<" "<<
			LDO_pair.first.second<<" "<<
			LDO_pair.second<<endl;
		*/
		table_ldo.insert(LDO_pair);
		// build ldo_vin list
		it = find(ldo_vin_vec.begin(), 
			ldo_vin_vec.end(), vin);
		if(it == ldo_vin_vec.end()){
			ldo_vin_vec.push_back(vin);
		}
		// build ldo_iout list
		it = find(ldo_iout_vec.begin(), 
			ldo_iout_vec.end(), iout);
		if(it == ldo_iout_vec.end()){
			ldo_iout_vec.push_back(iout);
		}
	}
	fclose(f);
	// sort ldo_iout_vec and vin_vec 
	sort(ldo_vin_vec.begin(), ldo_vin_vec.end());
	sort(ldo_iout_vec.begin(), ldo_iout_vec.end());
}

// extract voltage values of ldo nets (for SPICE)
void Circuit::extract_ldo_info(){
	vector <Node *> va;
	vector <Node *> vb;
	// extract va and vb values from solutions
	ckt_g.extract_ldo_vol(va);
	ckt_l.extract_ldo_vol(vb);
	for(size_t i=0;i<va.size();i++){
		clog<<"ldo node: "<<*va[i]<<" "<<*vb[i]<<endl;
	}
	// then call SPICE
	// first need to modify the Vin1 and Vin2 in the spice 
	// file with new values
	FILE *f;
	f = fopen("c_out.inc", "w");
	int id = 1;
	for(size_t i=0;i<va.size();i++){
		fprintf(f, "X%d Vout%d1 Vout%d2 ldo_subckt Vol_Vin1=%f Vol_Vin2=%f\n", id, id, id, va[i]->value, vb[i]->value);
		id++;
	}
	fclose(f);
	system("hspice ldo_top.spice");
	clog<<"after calling spice once. "<<endl;
}

// find correlated elements in the lookup table
void Circuit::find_table_elements(double &vin, 
	double &vin_1, double &vin_2, 
	vector<double> vec){
	size_t n = vec.size();
	// if out of range
	if(vin < vec[0]){
		vin_1 = vec[0];
		vin_2 = vin_1;
		vin = vin_1;
	}
	else if(vin > vec[n-1]){
		vin_1 = vec[n-1];
		vin_2 = vin_1;
		vin = vin_1;
	}
	else{
		for(size_t j=0;j<n-1;j++){
			if(vec[j] <= vin && 
				vec[j+1] >= vin){
				vin_1 = vec[j];
				vin_2 = vec[j+1];
				break;
			}
		}
	}
}

// first interpolate a line
double Circuit::inter_1D_table(double x, double x1, double y1, double x2, double y2){
	return ((x-x1)*y2 + (x2-x)*y1) / (x2-x1);
}
// perform interpolation operation
double Circuit::inter_table_ldo(double vin, double iout, double vin_1, double vin_2, double iout_1, double iout_2){
	double z11=0;
	double z12=0;
	double z21=0;
	double z22=0;
	double temp1 = 0;
	double temp2 = 0;

	z11 = table_ldo[make_pair(vin_1, iout_1)];
	z12 = table_ldo[make_pair(vin_2, iout_1)];
	z21 = table_ldo[make_pair(vin_1, iout_2)];
	z22 = table_ldo[make_pair(vin_2, iout_2)];
	// clog<<"z11, z12, z21, z22: "<<z11<<" "<<z12<<" "<<z21<<" "<<z22<<endl;

	// if only 1 element, return
	if(vin_1 == vin_2 && iout_1 == iout_2)
		return z11;
	// if 2 elements
	if(iout_1 != iout_2 && vin_1 != vin_2){
		temp1 = inter_1D_table(iout, 
			iout_1, z11, iout_2, z21);
		temp2 = inter_1D_table(iout, 
			iout_1, z12, iout_2, z22);
		double vout = inter_1D_table(vin, 
			vin_1, temp1, vin_2, temp2);
		return vout;
	}
	if(iout_1 == iout_2){
		double vout = inter_1D_table(vin, 
			vin_1, z11, vin_2, z12);
		return vout;
	}
	if(vin_1 == vin_2){
		double vout = inter_1D_table(iout, 
			iout_1, z11, iout_2, z21);
		return vout;
	}
	return z11;
}

double Circuit::locate_maxIRdrop(){
	max_IRdrop = ckt_l.locate_maxIRdrop();
	return max_IRdrop;
}

// perform DC LDO optimization
void Circuit::solve_DC_LDO(){
	map<Node*, double> ldo_best;
	pair<Node*, double> ldo_pair;
	vector<Node *> nd_out_vec;

	// first get IR drop values
	solve_DC();

	ckt_l.print_matlab_LDO();
	// ckt_l.print_matlab_node();
	double max_IRdrop = locate_maxIRdrop();
	clog<<"local initial max_IRdrop: "<<ckt_l.max_IRdrop<<endl;
	// clog<<"global max_IRdrop: "<<ckt_g.max_IRdrop<<endl;
	// cout<<ckt_l.nodelist<<endl;
	// cout<<ckt_g.nodelist<<endl;
	Node *nd = ckt_l.ldolist[0]->A;
	ldo_pair.first = nd;
	nd_out_vec.push_back(ckt_l.ldolist[0]->nd_out);
	ldo_pair.second = max_IRdrop;
	ldo_best.insert(ldo_pair);

	// clog<<"first max_IR: "<<max_IRdrop<<endl;
	double THRES = VDD_G * 0.1;
	// safe area, return with current LDO location
	// if(max_IRdrop <= THRES) return;

	// need to keep best ldolist in the trial
	// stores the best ldolist	
	for(int i=0;i<5;i++){
		// optimize the locations of LDO and rebuild nets
		relocate_LDOs();
		solve_DC();
		max_IRdrop = locate_maxIRdrop();
		Node *nd = ckt_l.ldolist[0]->A;

		clog<<"optimized max_IR: "<<max_IRdrop<<" "<<*nd<<endl;
		ldo_pair.first = nd;
		nd_out_vec.push_back(ckt_l.ldolist[0]->nd_out);
		ldo_pair.second = max_IRdrop;
		ldo_best.insert(ldo_pair);
	}
	map<Node*, double>::iterator it;
	Node *nd_min = NULL;
	double min_IR=0;
	int i=0;
	int j=0;
	for(it = ldo_best.begin();it!=ldo_best.end();it++){
		// clog<<"ldo best list: "<<*it->first<<" "<< it->second<<" "<<*nd_out_vec[i]<<endl;
		if(it ==  ldo_best.begin()){
			min_IR = it->second;
			nd_min = it->first;
			j=i;
		}else if(it->second < min_IR){
			min_IR = it->second;
			nd_min = it->first;
			j=i;
		}
		i++;
	}
	clog<<"best ldo: "<<*nd_min<<" "<<*nd_out_vec[j]<<endl;
	// switch to the best ldo
	ckt_l.ldolist[0]->nd_out = nd_out_vec[j];
	nd_out_vec.clear();
	recover_best_ldo(nd_min);
	// clog<<"after recover. "<<endl;
	ldo_best.clear();
	solve_DC();

	// cout<<ckt_l.nodelist<<endl;
	// cout<<endl;
	// cout<<ckt_g.nodelist<<endl;
	// max_IRdrop = locate_maxIRdrop();
	// clog<<"final max_IR with single LDO relocation: "<<max_IRdrop<<endl;
	// clog<<"final g max: "<<ckt_g.locate_g_maxIRdrop()<<endl;
	// clog<<"final ldo: "<<*ckt_l.ldolist[0]->nd_out<<" "<<*ckt_l.ldolist[0]->A<<endl;

	// ckt_l.print_matlab_LDO();
	// ckt_l.print_matlab_node();
	// clog<<"MAX_NUM_LDO: "<<ckt_l.MAX_NUM_LDO<<endl;
	// need to add more LDOs
	/*while(max_IRdrop > THRES){
		add_LDO_DC();
		clog<<"after add_LDO_DC. "<<endl;
	}*/
}

void Circuit::relocate_LDOs(){
	ckt_l.relocate_pads();
	// update ckt global nets with the change of LDO
	ckt_g.rebuild_global_nets();
}

// recover the SubCircuit with best ldo location
void Circuit::recover_best_ldo(Node *nd_min){
	// clog<<"nd_min: "<<*nd_min<<endl;
	// clog<<"ldo node: "<<*ckt_l.ldolist[0]->A<<endl;
	ckt_l.rebuild_local_nets(ckt_l.ldolist[0]->A, nd_min);
	clog<<"final best ldo: "<<*ckt_l.ldolist[0]->A<<" "<<*ckt_l.ldolist[0]->nd_out<<endl;
	ckt_g.rebuild_global_nets();
}

// add more LDO into circuit: new_ldo_flag = true
void Circuit::add_LDO_DC(){
	int iter=0;
	double diff_l = 1;
	double diff_g = 1;
	int local_bad_flag = 0;
	double THRES_l = VDD_G*0.1;
	// double THRES_g = 2.2*0.1;
	bool max_flag = false;
	double diff_old_l = 1;
	double diff_old_g = 1;
	double delta_diff_l = 1;
	double delta_diff_g = 0;

	// first atamp matrix: already symmetric
	ckt_l.stamp_decomp_matrix_DC();
	ckt_g.stamp_decomp_matrix_DC();
	// solve DC with fixed number of LDOs
	while((delta_diff_l > 1e-3 || delta_diff_g > 1e-3) && iter <10){
		// then update netlist
		ckt_l.modify_local_nets();
		ckt_l.stamp_rhs_DC(true);
		diff_old_l = diff_l;
		// solve eq with decomped matrix
		diff_l = ckt_l.solve_CK_with_decomp();
		delta_diff_l = fabs(diff_l - diff_old_l);
		clog<<"local max_IR: "<<ckt_l.locate_maxIRdrop()<<endl;
		if(ckt_l.max_IRdrop > THRES_l){
			local_bad_flag = 1;
			max_flag = optimize_local_LDO_DC(local_bad_flag, THRES_l);
			
			clog<<"optimized new local max_IR: "<<ckt_l.max_IRdrop<<" "<<ckt_l.ldolist.size()<<endl;
		}
		
		// calculate ldo current from ckt_l
		// update global current net
		ckt_l.update_ldo_current();
	
		// restamp global rhs with ldo current
		ckt_g.stamp_rhs_DC(false);
		diff_old_g = diff_g;
		// ckt_g.modify_ldo_rhs();
		
		diff_g = ckt_g.solve_CK_with_decomp();
		delta_diff_g = fabs(diff_g - diff_old_g);
		// then throw into ldo lookup table
		//update_ldo_vout();
		
		// clog<<"iter, diff_l, diff_g: "<<iter<<" "<<diff_l<<" "<<diff_g<<endl<<endl;
		iter++;
	}
	// locate_maxIRdrop();
}

// create new LDOs in circuit
void Circuit::create_new_LDOs(vector<Pad*> LDO_pad_vec){
	LDO *ldo_ptr;
	stringstream sstream;
	int pt_z = ldolist[0]->nd_in->pt.z;
	for(size_t i=0;i<LDO_pad_vec.size();i++){
		ldo_ptr = new LDO();
		ldo_ptr->A = LDO_pad_vec[i]->node;
		ldo_ptr->nd_out = LDO_pad_vec[i]->node;
		ldo_ptr->width = ldolist[0]->width;
		ldo_ptr->height = ldolist[0]->height;
		// build nd_in node
		sstream.str("");
		sstream<<"n"<<pt_z<<"_"<<ldo_ptr->A->pt.x<<"_"<<ldo_ptr->A->pt.y;
		ldo_ptr->nd_in = ckt_g.get_node(sstream.str());
		// ldolist.push_back(ldo_ptr);
		// synchronize ckt_l and ckt_g with ckt
		ckt_l.ldolist.push_back(ldo_ptr);
		// ckt_g.ldolist.push_back(ldo_ptr);
	}
	// clog<<"ckt ldolist size: "<<ldolist.size()<<endl;
	// clog<<"ckt_l ldolist size: "<<ckt_l.ldolist.size()<<endl;
}

#if 0
// working on add LDOs to TR step
void Circuit::add_LDO_all_global(Tran &tran, double time){
	// clog<<"start add LDO tr global. "<<endl;
	// temporary storing newly added LDOs
	vector<Pad*> LDO_pad_vec;
	Node *nd;
	for(size_t i=0;i<ckt_g.candi_pad_set.size();i++){
		nd = ckt_g.candi_pad_set[i]->node;
		if(nd->isS()!=Y)
			LDO_pad_vec.push_back(ckt_g.candi_pad_set[i]);
	}
		
	// find the node where new LDOs shoudl go to
	// ckt_g.extract_add_pad_dc_info(LDO_pad_vec, local_bad_flag);
	// if no room to add new LDO pad, return
	// if(LDO_pad_vec.size()==0)
		// return;
	
	// rebuild local and global net
	ckt_g.create_global_pad_new_nets(LDO_pad_vec);
	ckt_g.build_pad_set();
	/*for(size_t i=0;i<ckt_l.pad_set.size();i++)
		clog<<"i, ckt_l pad set: "<<i<<" "<<*ckt_l.pad_set[i]->node<<" "<<ckt_l.pad_set[i]->node->isS()<<endl;
		*/
	solve_global(tran, time);
	// clog<<"after solve_ITR. "<<endl;
	max_IRdrop = locate_maxIRdrop();
	// clog<<"best global max_IR drop for TR: "<<max_IRdrop<<endl;
	LDO_pad_vec.clear();
}
#endif

#if 0
// working on add LDOs to TR step
void Circuit::add_LDO_TR_global(Tran &tran, double time, bool local_bad_flag){
	// clog<<"start add LDO tr global. "<<endl;
	// temporary storing newly added LDOs
	vector<Pad*> LDO_pad_vec;
	// find the node where new LDOs shoudl go to
	ckt_g.extract_add_pad_dc_info(LDO_pad_vec, local_bad_flag);
	// if no room to add new LDO pad, return
	if(LDO_pad_vec.size()==0)
		return;
	
	// rebuild local and global net
	ckt_g.create_global_pad_new_nets(LDO_pad_vec);
	ckt_g.build_pad_set();
	/*for(size_t i=0;i<ckt_l.pad_set.size();i++)
		clog<<"i, ckt_l pad set: "<<i<<" "<<*ckt_l.pad_set[i]->node<<" "<<ckt_l.pad_set[i]->node->isS()<<endl;
		*/
	solve_global(tran, time);
	// clog<<"after solve_ITR. "<<endl;
	max_IRdrop = locate_maxIRdrop();
	// clog<<"final max_IR drop for TR: "<<max_IRdrop<<endl;
	LDO_pad_vec.clear();
}
#endif

// solve local
void Circuit::solve_global(Tran &tran, double time){
	ckt_g.reset_b();
	ckt_g.stamp_decomp_matrix_TR(tran);
	// Ieq already there
	// ckt_g.modify_global_nets();
	ckt_g.stamp_rhs_tr(false, time, tran);
	// solve eq with decomped matrix
	ckt_g.solve_CK_with_decomp();
	ckt_g.locate_g_maxIRdrop();
	// clog<<"new global max IR is: "<<ckt_g.locate_g_maxIRdrop()<<endl;
}


// working on add LDOs to TR step
void Circuit::add_LDO_TR_local(Tran &tran, double time){
	// clog<<"initial max_IR: "<<ckt_l.locate_maxIRdrop()<<endl;
	// temporary storing newly added LDOs
	vector<Pad*> LDO_pad_vec;
	// find the node where new LDOs shoudl go to
	ckt_l.extract_add_LDO_dc_info(LDO_pad_vec, tran);
	// if no room to add new LDO pad, return
	if(LDO_pad_vec.size()==0){
		return;
	}
	// rebuild local and global net
	ckt_l.create_local_LDO_new_nets(LDO_pad_vec);
	ckt_g.create_global_LDO_new_nets(LDO_pad_vec);

	// create new LDOs in circuit
	create_new_LDOs(LDO_pad_vec);
	// recover all the ldo voltages to 1.8V
	for(size_t i=0;i<ckt_l.ldolist.size();i++){
		ckt_l.ldolist[i]->A->value = VDD_G;
	}
	ckt_l.add_pad_set(LDO_pad_vec);
	ckt_l.solve_local(tran, time);
	//ckt_l.solve_GS(tran);
	//clog<<"after solve GS. "<<endl;
	max_IRdrop = locate_maxIRdrop();

	Node *nd = ckt_l.extract_maxIR_node();
	//clog<<"optimized max_IR drop for TR: "<<max_IRdrop<<" "<<*nd<<endl;
	LDO_pad_vec.clear();
}

// solve local
void Circuit::solve_local(Tran &tran, double time){
	ckt_l.stamp_decomp_matrix_TR(tran);
	// Ieq already there
	ckt_l.modify_local_nets();
	ckt_l.stamp_rhs_tr(true, time, tran);
	// solve eq with decomped matrix
	// ckt_l.solve_CK_with_decomp_tr();
	ckt_l.solve_CK_op_tr();
	ckt_l.locate_maxIRdrop();
	// clog<<"new local max IR in add LDO is: "<<ckt_l.locate_maxIRdrop()<<endl;
}

// working on add LDOs to TR step
void Circuit::add_LDO_TR(Tran &tran, double time){
	// temporary storing newly added LDOs
	vector<Pad*> LDO_pad_vec;
	// find the node where new LDOs shoudl go to
	//ckt_l.extract_add_LDO_dc_info(LDO_pad_vec);
	// clog<<"LDO_pad size: "<<LDO_pad_vec.size()<<endl;
	// if no room to add new LDO pad, return
	if(LDO_pad_vec.size()==0)
		return;
	
	// rebuild local and global net
	ckt_l.create_local_LDO_new_nets(LDO_pad_vec);
	ckt_g.create_global_LDO_new_nets(LDO_pad_vec);

	// for(size_t i=0;i<LDO_pad_vec.size();i++)
		// clog<<"i, LDO_pad: "<<i<<" "<<*LDO_pad_vec[i]->node<<endl;
	// create new LDOs in circuit
	create_new_LDOs(LDO_pad_vec);
	// clog<<"orig ckt_l.pad_set.size(): "<<ckt_l.pad_set.size()<<endl;
	ckt_l.add_pad_set(LDO_pad_vec);
	// clog<<"new ckt_l.pad_set.size(): "<<ckt_l.pad_set.size()<<endl;
	/*for(size_t i=0;i<ckt_l.pad_set.size();i++)
		clog<<"i, ckt_l pad set: "<<i<<" "<<*ckt_l.pad_set[i]->node<<" "<<ckt_l.pad_set[i]->node->isS()<<endl;
		*/
	// clog<<"create new LDOs. "<<endl;
	solve_TR(tran, time);
	// clog<<"after solve_ITR. "<<endl;
	max_IRdrop = locate_maxIRdrop();
	clog<<"final max_IR drop for TR: "<<max_IRdrop<<endl;
	LDO_pad_vec.clear();
}

void Circuit::verify_one_LDO_step(Tran &tran, double time){
	int iter=0;
	double diff_l = 1;
	double diff_g = 1;
	// stores the Ieq for both C and L nets 
	ckt_l.modify_rhs_tr_0(tran);
	ckt_g.modify_rhs_tr_0(tran);
	double diff_old_l = 1;
	double diff_old_g = 1;
	double delta_diff_l = 1;
	double delta_diff_g = 0;

	while((delta_diff_l > 1e-4 || delta_diff_g > 1e-4)){
		// clog<<endl<<"time: "<<time<<" iter: "<<iter<<endl;
		// update ldo vol net values  
		ckt_l.modify_local_nets();
		// restamp bp
		ckt_l.stamp_rhs_tr(true, time, tran);
		diff_old_l  = diff_l;
		// solve eq with decomped matrix
		diff_l = ckt_l.solve_CK_with_decomp_tr();
		delta_diff_l = fabs(diff_l - diff_old_l);
		ckt_l.locate_maxIRdrop();
		// calculate ldo current from ckt_l
		ckt_l.update_ldo_current();
		// clog<<"ckt_l IR: "<<ckt_l.locate_maxIRdrop()<<endl;

		// modify global bp
		ckt_g.stamp_rhs_tr(false, time, tran);
		diff_old_g = diff_g;
		diff_g = ckt_g.solve_CK_with_decomp_tr();
		delta_diff_g = fabs(diff_g - diff_old_g);
		ckt_g.locate_g_maxIRdrop();

		// clog<<"ckt_g IR: "<<ckt_g.locate_maxIRdrop()<<endl;
		// clog<<"iter, diff_l, g: "<<iter<<" "<<diff_l<<" "<<diff_g<<endl;
		// cout<<ckt_l.nodelist<<endl;
		// cout<<ckt_g.nodelist<<endl;
		// then throw into ldo lookup table
		// update_ldo_vout();
		iter++;
	}
	// save ckt_nodes value for t step
	ckt_l.save_ckt_nodes();
	// locate_maxIRdrop();		
	// clog<<"max IR after solve_TR is: "<<max_IRdrop<<endl;
}

bool Circuit::optimize_local_LDO(int flag, double THRES_l, Tran &tran, double time){
	if(flag != 1)
		return false;
	bool max_flag = false;
	// add LDO to lcoal grid
	int iter_i = 0;
	while(flag ==1 && iter_i <5){
		if((int)ckt_l.ldolist.size() >= ckt_l.MAX_NUM_LDO){
			max_flag = true;
			return max_flag;
			break;
		}
		// clog<<endl<<"iter_i: "<<iter_i<<" "<<ckt_l.ldolist.size()<<" "<<ckt_l.locate_maxIRdrop()<<endl;
		double temp_IR = ckt_l.locate_avgIRdrop();
		add_LDO_TR_local(tran, time);
		/*double op_IR = ckt_l.locate_avgIRdrop();
		if(fabs(op_IR - temp_IR) < 1e-4)
			break;
			*/

		// clog<<"local opti IR: "<<ckt_l.max_IRdrop<<endl;
		if(ckt_l.locate_maxIRdrop() <THRES_l){
			flag = 0;
			break;
		}
		iter_i++;
	}
	// ckt_l.update_ldo_current();
	return false;
}

#if 0
// add power pads into global grid
void Circuit::optimize_global_LDO(int flag, double THRES_g, Tran &tran, double time, double THRES_l){
	// if(flag !=2) return;
	// clog<<"enter optimize global"<<endl;
	// add LDO to lcoal grid
	int iter_i = 0;
	// while(flag ==2 || iter_i <3){
	while(ckt_l.max_IRdrop>THRES_l || ckt_g.max_IRdrop >THRES_g || iter_i <3){
		if((int)ckt_g.pad_set.size() >= ckt_g.MAX_NUM_PAD)
			break;
		// clog<<endl<<"iter_i: "<<iter_i<<" "<<ckt_g.locate_g_maxIRdrop()<<endl;;
		add_LDO_TR_global(tran, time, flag);
		// ckt_g.max_IRdrop = ckt_g.locate_g_maxIRdrop();
		if(flag ==2 && ckt_g.locate_g_maxIRdrop() <= THRES_g)
			flag = 0;
		if(flag ==1 || ckt_g.max_IRdrop < THRES_g){
			update_ldo_vout();
			solve_local(tran, time);
			// ckt_l.max_IRdrop = ckt_l.locate_maxIRdrop();
			// clog<<"local max IR is: "<<ckt_l.locate_maxIRdrop()<<endl;
			if(ckt_l.locate_maxIRdrop() < THRES_l)
				break;
		}
		iter_i++;
	}
}
#endif
bool Circuit::optimize_local_LDO_DC(int flag, double THRES_l){
	if(flag != 1)
		return false;
	bool max_flag = false;
	// add LDO to lcoal grid
	int iter_i = 0;
	while(flag ==1 && iter_i <1){
		if((int)ckt_l.ldolist.size() >= ckt_l.MAX_NUM_LDO){
			max_flag = true;
			return max_flag;
			break;
		}
		// clog<<endl<<"iter_i: "<<iter_i<<" "<<ckt_l.ldolist.size()<<" "<<ckt_l.locate_maxIRdrop()<<endl;
		//add_LDO_DC_local();
		// clog<<"opti IR: "<<ckt_l.max_IRdrop<<endl;
		if(ckt_l.locate_maxIRdrop() <= THRES_l){
			flag = 0;
			break;
		}
		iter_i++;
	}
	ckt_l.update_ldo_current();
	return false;
}

// working on add LDOs to TR step
void Circuit::add_LDO_DC_local(Tran tran){
	// temporary storing newly added LDOs
	vector<Pad*> LDO_pad_vec;
	// find the node where new LDOs shoudl go to
	ckt_l.extract_add_LDO_dc_info(LDO_pad_vec, tran);
	// if no room to add new LDO pad, return
	if(LDO_pad_vec.size()==0){
		// clog<<"no add new ldo. "<<endl;
		return;
	}
	// clog<<"new LDO: "<<*LDO_pad_vec[0]->node<<endl;	
	// rebuild local and global net
	ckt_l.create_local_LDO_new_nets(LDO_pad_vec);
	ckt_g.create_global_LDO_new_nets(LDO_pad_vec);

	// create new LDOs in circuit
	create_new_LDOs(LDO_pad_vec);
	ckt_l.add_pad_set(LDO_pad_vec);
	solve_local_DC();
	max_IRdrop = locate_maxIRdrop();
	ckt_l.update_ldo_current();
	// clog<<"final max_IR drop for TR: "<<max_IRdrop<<endl;
	LDO_pad_vec.clear();
}

// solve local
void Circuit::solve_local_DC(){
	ckt_l.reset_b();
	ckt_l.stamp_decomp_matrix_DC();
	// Ieq already there
	ckt_l.modify_local_nets();
	ckt_l.stamp_rhs_DC(true);
	// solve eq with decomped matrix
	ckt_l.solve_CK_with_decomp();
	ckt_l.locate_maxIRdrop();
	// cout<<"new diff_l, local max IR is: "<<diff_l<<" "<<ckt_l.locate_maxIRdrop()<<endl;
}




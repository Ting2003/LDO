// ----------------------------------------------------------------//
// Filename : circuit.h
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// declaration of Circuit class
// used to construct the circuit network
// ----------------------------------------------------------------//
// - Ting Yu - Tue Feb 8 5:45 pm 2011
//   * added the ostream<< func in .h file
// - Zigang Xiao - Tue Jan 18 21:46:13 CST 2011
//   * added solve() and related function
// - Zigang Xiao - Sun Jan 16 16:04:03 CST 2011
//   * added this log

#ifndef __SUBCIRCUIT_H__
#define __SUBCIRCUIT_H__

#include <string>
#include <vector>
#include <tr1/unordered_map>
#include <map>
#include <cmath>
#include "global.h"
#include "node.h"
#include "net.h"
#include "vec.h"
#include "transient.h"
#include "cholmod.h"
#include <algorithm>
#include "sp_graph_table.h"
#include "sp_node.h"
#include "pad.h"
#include "ldo.h"
#include <queue>

using namespace std;
using namespace std::tr1;

typedef vector<double> DoubleVector;
typedef vector<Net *> NetPtrVector;
typedef vector<Node *> NodePtrVector;
typedef NetPtrVector NetList;

// functor of translating Node * to void *
namespace std{ 
	namespace tr1{
		template<> struct hash< Node * >{
			size_t operator()( const Node * x ) const {
				return hash< const char* >()( (char*) x );
			}
		};
	}
}

class SubCircuit{
public:
	SubCircuit(string name="");
	~SubCircuit();
	void check_sys() const;
	// can be written as inline to speed up
	Node * get_node(string name);
	Net * get_net(string name);
	string get_name() const;

	static size_t get_total_num_layer();

	// add a node into nodelist
	bool add_node(Node * nd);

	// add a net into netset
	bool add_net(Net * net);

	bool has_node(string name) const;
	bool has_net(string name) const;

	// sort nodes according to predefined order
	void sort_nodes();
        
        void make_A_symmetric(double *bp);
	void make_A_symmetric_local(double *b);
	void make_A_symmetric_tr(double *b, Tran &tran);
	//void make_A_symmetric_block();
      //void make_A_symmetric_block_tr(Tran &tran);

      // solve for node voltage
	void solve(Tran &tran, bool flag);
	
	//void set_blocklist(Node * nd);
	friend ostream & operator << (ostream & os, const SubCircuit & ckt);
	friend class Parser;
	friend class Circuit;

	///// new functions for pad /////
	double locate_maxIRdrop();
	double locate_g_maxIRdrop();
	double locate_maxIRdrop(double *x, size_t n);
	double locate_maxIRdrop_tr(Tran &tran);
	double locate_special_maxIRdrop();
	// void mark_special_nodes();
	bool set_ldo(double ref_dist, double ref_x, double ref_y, double &x0, double &y0, LDO &ldo);
	void recover_local_pad(Tran &tran, vector<LDO*> &ldolist_best);
	void recover_global_pad(Tran &tran, vector<Node*> &pad_set_best);
	bool adjust_ldo_pos(double ref_dist, double ref_x, double ref_y, LDO &ldo, MODULE *wspace);
	void build_pad_set();
	void get_pad_tr_cur(vector<Pad*> &pad_set, Tran &tran);
	////// new member for pad //////
	
	double max_IRdrop;
	vector<Pad*> pad_set;
	vector<Node*> origin_pad_set;
	vector<Pad*> candi_pad_set;
	////// new functions for pad /////
	void assign_distance(Node *nds, Node *nd, double dist);
	void print_pad_map(vector<Pad*> &pad_set);
	void release_resource();
	void clear_flags();
	double update_pad_pos(double ref_drop_value, size_t i);
	double update_pad_pos_all(vector<Pad*> pad_set);
	void round_data(double &data);
	Node * pad_projection(Pad *pad, Node *nd);
	void project_pads(unordered_map<string, Node*> map_node_pt, vector<Pad*> &pad_set);
	void build_ldolist(vector<LDO*> ldo_vec);
	void build_wspacelist(vector<MODULE*> wspace_vec);
	void relocate_pads();
	void relocate_pads_graph_global(Tran &tran, vector<LDO*> &ldo_vec, vector<MODULE*> &wspace_vec);
	void relocate_pads_graph(Tran &tran, vector<LDO*> &ldo_vec, vector<MODULE*> &wspace_vec);
	void relocate_pads(Tran &tran, vector<LDO*> &ldolist, vector<MODULE*> &wspace_vec);

	void restore_pad_set(vector<Pad*> &pad_set, vector<Node*>&pad_set_old);
	void assign_pad_set(vector<Pad*> pad_set, vector<Node*>&pad_set_old);
	void rebuild_local_nets(Node *rm_node, Node *add_node);
	void rebuild_global_nets();
	void modify_local_nets();
	void modify_global_nets();

	void build_candi_graph();
	double compute_stand_dev();
	void rebuild_voltage_nets_g(vector<Pad*> pad_set, vector<Node*> origin_pad_set);
	void rebuild_voltage_nets_l(vector<Pad*> pad_set, vector<Node*> origin_pad_set);
	void compute_ldo_current();
	void update_ldo_voltage(char *filename);
	void verify_ldo(Tran &tran, char *filename);
	void print_pad_set(vector<Pad*> &pad_set);
	void extract_pads(vector<Pad*> pad_set);
	Pad* locate_candi_pad_maxIR(vector<Pad*> pad_set);
	void print_matlab();
	void clear_pad_control_nodes(vector<Pad*> &pad_set);
	void update_pad_control_nodes(vector<Pad*> pad_set);
	void extract_add_LDO_dc_info(vector<Pad*> &LDO_pad_vec);
	void create_current_LDO_graph();
	void create_local_LDO_new_nets(vector<Pad*> LDO_pad_vec);
	void create_global_LDO_new_nets(vector<Pad*> LDO_pad_vec);
	void extract_min_max_pads_new(double VDD, vector<double> ref_drop_vec, bool local_flag);

	void build_pad_graph(vector<Pad*> &pad_set);
	void modify_graph(bool flag);
	void print_ldo_list();
	Pad *find_nbr_pad(vector<Pad*> &pad_set, Pad *pad);
	double get_distance(Node *na, Node *nb);
	void graph_move_pads();
	int locate_max_drop_pad(vector<Pad*> pad_set);
	double calc_avg_ref_drop(vector<double> &ref_drop_vec);
	double calc_avg_ref(vector<Pad*> &pad_set, vector<double> ref_drop_vec);
	double locate_ref(vector<Pad*> pad_set, size_t i);
	void update_single_pad_flag(Pad* pad);
	bool print_flag(Node *nd);
	void modify_newxy(vector<Pad*> &pad_set);
	void update_node(Net *net);
	void print_all_control_nodes();
	//////// end functions for pad ////////

	// C style output
	void print();

private:
	// member functions
	// void solve_LU(Tran &tran, bool flag);
	// void solve_LU_core(Tran &tran, bool flag);
	void build_local_nets();
	void build_global_nets();
	void stamp_decomp_matrix_DC(bool local_flag);
	void stamp_decomp_matrix_TR(Tran &tran, double time, bool local_flag);

	double solve_CK_with_decomp();
	double solve_CK_with_decomp_tr();
	void update_ldo_current();
	void modify_ldo_rhs();
	void modify_ldo_rhs_TR();
		
	// add solve with ADI method
	//void solve_partial_ADI();
	//double solve_ADI_DC();
	// initialize things before solve_iteration
	void solve_init(bool flag);
	void configure_init();
	void reconfigure_TR();
	void mark_geo_occupation();
	void extract_node(char * str, Node & nd);

	void count_merge_nodes();

	// methods of stamping the matrix
	void stamp_by_set(Matrix & A);
	void stamp_resistor(Matrix & A, Net * net);
	void stamp_current(double * b, Net * net);
	void stamp_VDD(Matrix & A, Net * net);
	void stamp_rhs_DC(bool local_flag);
	void stamp_rhs_VDD(double *bp, Net *net);
	void stamp_induc_rhs_dc(double *b, Net * net);
	void release_resources();
	void stamp_inductance_dc(Matrix & A, Net * net);
	void stamp_capacitance_dc(Matrix & A, Net * net);

	void stamp_by_set_tr(Matrix & A, Tran &tran);
	void stamp_rhs_tr(bool local_flag, double time, Tran &tran);
	void stamp_resistor_tr(Matrix & A, Net * net);
	void current_tr(Net *net, double &time);
	
	void stamp_current_tr_1(double *bp, double *b, double &time);
	void stamp_current_tr_net_1(double *bp, double *b, Net *net, double &time);

	void stamp_current_tr(double *b, double &time);
	void stamp_current_tr_net(double *b, Net * net, double &time);
	void stamp_capacitance_tr(Matrix & A, Net * net, Tran &tran);
	void stamp_inductance_tr(Matrix & A, Net * net, Tran &tran);
	void modify_rhs_tr_0(Tran &tran);

	void modify_rhs_tr(double *b, double *xp);
	void set_eq_induc(Tran &tran);
	void set_eq_capac(Tran &tran);
	void modify_rhs_c_tr_0(Net *net, Tran &tran);
	void modify_rhs_l_tr_0(Net *net, Tran &tran);

	void modify_rhs_c_tr(Net *net, double *xp, double *x);
	void modify_rhs_l_tr(Net *net, double *xp, double *x);
	void modify_rhs_Ieq(double *rhs);
	void modify_rhs_Ieq_c(Net *net, double *rhs);
	void modify_rhs_Ieq_l(Net *net, double *rhs);

	void reset_b();

	bool qualify_pad(Node *nd_new, LDO *ldo, unordered_map<string, Node*> map_node_pt);
	void locate_ldo_region_bound(int a, int b, int &min, int &max);
	Node * project_local_pad(Node *nd, Node *nd_new, LDO *ldo);
	void project_ldo_node(int &ref_x, int &ref_y, LDO &ldo);
	bool node_in_ldo_or_block(double ref_x, double ref_y);
	Node * expand_ldo_location(double ref_dist, int ref_x, int ref_y, LDO &ldo_ptr);
	Node * expand_pad(Node *nd_new, LDO *ldo, unordered_map<string, Node*> map_node_pt);
	void get_candi_wspace(double ref_x, double ref_y, double ref_dist, vector<int> &candi_wspace);
	//bool place_ldo(double ref_dist, double ref_x, double ref_y, LDO &ldo_ptr, vector<int> &candi_wspace);

	void release_tr_nodes(Tran &tran);
	void release_ckt_nodes();
	void print_ckt_nodes(Tran &tran);
	void save_ckt_nodes_to_tr(Tran &tran);
	void link_tr_nodes(Tran &tran);
	void link_ckt_nodes(Tran &tran);
	void save_tr_nodes(Tran &tran, double *x);
	void save_ckt_nodes();

	void print_tr_nodes(Tran &tran);

	void copy_node_voltages(double *x, size_t &size, bool from=true);

	// after solving, copy node voltage from replist to nodes
	double get_voltages_from_LU_sol(double *x);
	void select_omega();

	void set_type(CIRCUIT_TYPE type){circuit_type = type;};
        // ************* functions and members for thread **********

        double *temp;	
        int *id_map;
	Matrix A;
        cholmod_factor *L;
	double *Lx;
	int *Li, *Lp, *Lnz;
        cholmod_common c, *cm;

        cholmod_dense *b, *x;
        double *bp, *xp;

	// ********* sparse vectors ******
	int flag_ck;
	vector<Node_TR_PRINT> ckt_nodes;
	// ************** member variables *******************
	NodePtrVector nodelist;		// a set of nodes
	NodePtrVector replist;		// a set of representative nodes
	vector<LDO*> ldolist;
	vector<MODULE*> wspacelist;
	NetPtrVector net_set[NUM_NET_TYPE];// should be the same as size of NET_TYPE
	// defines the net direction in layers
	static vector<LAYER_DIR> layer_dir;
	vector<int> layers;
	vector<bool> local_layers;
	vector<bool> global_layers;
	
	// mapping from name to Node object pointer
	unordered_map<string, Node*> map_node;
	int MAX_NUM_LDO;

	// mapping from Net name to object pointer
	// unordered_map<string, Net*> map_net;

	// mapping from Node pointer to their index in nodelist
	unordered_map<Node *, size_t> node_id;
	//unordered_map<Node *, size_t> rep_id;

	// circuit name
	string name;
	// boundary of the circuit
	double lx, ly, gx, gy;

	// blocks
	//BlockInfo block_info;
	//size_t x_min, y_min, x_max, y_max;

	// control variables

	CIRCUIT_TYPE circuit_type;

	double VDD;
	//size_t num_blocks;
};

inline size_t SubCircuit::get_total_num_layer(){return layer_dir.size();}

// adds a node into nodelist
inline bool SubCircuit::add_node(Node * node){
	nodelist.push_back(node);
	map_node[node->name] = node;
	return true;
}

// adds a net into netset
inline bool SubCircuit::add_net(Net * net){
	//map_net[net->name] = net;
	net_set[net->type].push_back(net);
	return true;
}

// fina a node by name
inline bool SubCircuit::has_node(string name) const{
	if( map_node.find(name) != map_node.end() ) return true;
	return false;
}

// get a node by name
inline Node * SubCircuit::get_node(string name){
	unordered_map<string, Node*>::const_iterator it = map_node.find(name);
	if( it != map_node.end() ) return it->second;
	else return NULL;
}

/*
// find a net by name
inline bool Circuit::has_net(string name) const{
	if( map_net.find(name) != map_net.end() ) return true;
	return false;
}


// get a net by name
inline Net * Circuit::get_net(string name){return map_net[name];}
*/


ostream & operator << (ostream & os, const NodePtrVector & nodelist);
ostream & operator << (ostream & os, const NetPtrVector & nets);
//ostream & operator << (ostream & os, const vector<Block > & block_info);
#endif

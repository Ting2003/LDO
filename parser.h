// ----------------------------------------------------------------//
// Filename : parser.h
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// declaration of Parser class
// ----------------------------------------------------------------//
// - Zigang Xiao - Sun Jan 16 16:04:03 CST 2011
//   * added this log

#ifndef __PARSER_H__
#define __PARSER_H__

#include <vector>
#include "global.h"
#include "circuit.h"
#include "transient.h"
#include "ckt_top.h"
using std::vector;

// given an input file, parse it and store corresponding result into Circuit
class Parser{
public:
	// supply Circuit objects
	Parser(vector<CKT_TOP*> * ckts, map<Circuit*, CKT_TOP*> *map_top);
	~Parser();

	// parser a input file and construct the circuit
	void parse(char * filename, Tran &tran);

	int get_num_layers() const;

private:
	int create_circuits();		// parse the file and create circuits

	void try_change_via(Net *);

	//void insert_net_node(string line);
	void insert_net_node(char * line, int *count);
	void extract_node(char * str, Node & nd);
	void update_node(Net * net);
	void parse_dot(char *line, Tran &tran);
	void add_node(Circuit *ckt, Node &nd, Node *nd_ptr);

	char * filename;		  // input file name
	int n_layer;			  // total number of layers
	map<int, Circuit*> layer_map_ckt;	  // which circuit a layer belong
	vector<CKT_TOP*> * p_ckts;	  // pointer to circkt list
	map<Circuit*, CKT_TOP*> map_ckt_top;
};

// Trick: try to modify the net
inline void Parser::try_change_via(Net * net){
	// is it a via?
	 if( net->ab[0]->get_layer() == net->ab[1]->get_layer() )
		return;

	// make it a zero voltage via
	if( net->type == RESISTOR && net->value < 1e-4 ){
		net->type = VOLTAGE;
		net->value = 0.0;
	}
}

#endif

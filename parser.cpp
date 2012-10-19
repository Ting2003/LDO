// ----------------------------------------------------------------//
// Filename : parser.cpp
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// implementation file of parser.h
// ----------------------------------------------------------------//
// - Zigang Xiao - Sun Jan 16 16:04:03 CST 2011
//   * added this log

#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstring>
#include "util.h"
#include "parser.h"
using namespace std;

// store the pointer to circuits
Parser::Parser(vector<CKT_TOP> * ckts):n_layer(0), p_ckts(ckts){
}

// Trick: do not descruct to speed up
Parser::~Parser(){
	//for(size_t i=0;i<(*p_ckts).size();i++) delete (*p_ckts)[i];
}

// _X_n2_19505_20721 
// X:
// n2: layer 2
// 19505 20721: coordinate
void Parser::extract_node(char * str, Node & nd){
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

// given a line, extract net and node information
void Parser::insert_net_node(char * line, int *count){
	char *chs, *saveptr;
	const char* sep = " (),\n";
	char sname[MAX_BUF];
	char sa[MAX_BUF];
	char sb[MAX_BUF];
	static Node nd[2];
	Node * nd_ptr[2];	// this will be set to the two nodes found
	double value;
	int ckt_id;
	sscanf(line, "%s %s %s %lf", sname, sa, sb, &value);
	clog<<endl<<"origin net: "<<line;	
	extract_node(sa, nd[0]);
	extract_node(sb, nd[1]);	

	NET_TYPE net_type = RESISTOR;
	// find net type
	switch(sname[0]){
	case 'r': // resistor
	case 'R':
		net_type = RESISTOR;
		break;
	case 'v': // VDD
	case 'V':
		net_type = VOLTAGE;
		break;
	case 'i': // current
	case 'I':
		net_type = CURRENT;
		break;
	case 'c': // capacitance
	case 'C':
		net_type = CAPACITANCE;
		break;
	case 'l':
	case 'L':
		net_type = INDUCTANCE;
		break;
	case 'x':
	case 'X':
		net_type = LDO;
		break;
	default:
		report_exit("Invalid net type!\n");
		break;
	}

	int layer;
	int i = 0;
	Circuit *ckt;

	// create a Net
	//Net * net = new Net(net_type, name, value, nd_ptr[0], nd_ptr[1]);
	Net * net = new Net(net_type, value, nd_ptr[0], nd_ptr[1]);
	cout<<net_type<<" "<<value<<" "<<*nd_ptr[0]<<" "<<*nd_ptr[1]<<endl;
	cout<<"new net: "<<*net<<endl;	

	// assign pulse paramter for pulse input
	chs = strtok_r(line, sep, &saveptr);
	for(int i=0;i<3;i++)
		chs = strtok_r(NULL, sep, &saveptr);
	if(chs != NULL)
		chs = strtok_r(NULL, sep, &saveptr);
	if(chs != NULL){
		chs = strtok_r(NULL, sep, &saveptr);
		net->V1 = atof(chs);
		chs = strtok_r(NULL, sep, &saveptr);
		net->V2 = atof(chs);
		chs = strtok_r(NULL, sep, &saveptr);
		net->TD = atof(chs);
		chs = strtok_r(NULL, sep, &saveptr);
		net->Tr = atof(chs);
		chs = strtok_r(NULL, sep, &saveptr);
		net->Tf = atof(chs);
		chs = strtok_r(NULL, sep, &saveptr);
		net->PW = atof(chs);
		chs = strtok_r(NULL, sep, &saveptr);
		net->Period = atof(chs);
	}
	// trick: when the value of a resistor via is below a threshold,
	// treat it as a 0-voltage via
	try_change_via(net);

	// modify the nets and insert into two circuits
	if(nd[0].name == "0" || nd[1].name == "0" || nd[0].get_layer() == nd[1].get_layer()){
		if(nd[0].name == "0"){
			layer = nd[1].get_layer();
			ckt = layer_map_ckt[layer];
			i = 1;
		}
		/*else if(nd[1].name == "0"){
			layer = nd[0].get_layer();
			ckt = layer_map_ckt[layer];
			i = 0;
		}*/
		else{
			layer = nd[0].get_layer();
			ckt = layer_map_ckt[layer];
			i = 0;
		}
		// insert this net into circuit
		cout<<"net: "<<*net<<endl;	
		add_node(ckt, nd[i], nd_ptr[i]);
		ckt->add_net(net);
		cout<<"add net: "<<*net<<endl;
	}
	else{ // boundary nets
		for(int i=0;i<2;i++){
			layer = nd[i].get_layer();
			ckt = layer_map_ckt[layer];
			add_node(ckt, nd[i], nd_ptr[i]);
			// both add net
			ckt->add_net(net);
		}
	}	
		
	// IMPORTANT: set the relationship between node and net
	// update node voltage if it is an X node
	// set node to be X node if it connects to a voltage source
	update_node(net);
}

// Given a net with its two nodes, update the connection information
// for thet two nodes
void Parser::update_node(Net * net){
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
		if(a->is_ground()) swap<Node*>(a,b);
		a->flag = Z;
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
		&& b->isS() != -1){
		if(a->isS()== X){
			a->set_nbr(BOTTOM, net);
			b->set_nbr(TOP, net);
		}
		else if(b->isS() == X){
			b->set_nbr(BOTTOM, net);
			a->set_nbr(TOP, net);
		}	
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
			Circuit::layer_dir[layer] = HR;
		}
		else if( a->pt.x == b->pt.x ){// vertical
			if(a->pt.y > b->pt.y) swap<Node*>(a,b);
			a->set_nbr(NORTH, net);
			b->set_nbr(SOUTH, net);
			Circuit::layer_dir[layer] = VT;
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
		a->flag = X;		// set a to be X node
		a->set_nbr(TOP, net);	// X -- VDD -- Ground
		a->set_value(net->value);
	}
	
	else{// if( net->type == CURRENT ){// current source
		// let a be ground node
		if( !a->is_ground() ) swap<Node*>(a,b);
		b->set_nbr(BOTTOM, net);
	}
}

// parse the file and create circuits
int Parser::create_circuits(){
	FILE * fp;		// used for popen/pclose
	int status;		// return status of popen/pclose
	const char grep[]="grep 'layer' ";
	//const char rest[]="|sort -t ',' -k 2 -r |cut -d ',' -f 2 |cut -d ' ' -f 1,3";
	const char rest[]="|sort -t ',' -k 2 -r |cut -d ',' -f 2 |cut -d ' ' -f 1,3,4";

	char cmd[MAX_BUF], name[MAX_BUF]="", GL[MAX_BUF] = "";
	string prev_ckt_name("");
	int layer, n_circuit=0;

	// extract useful information about layers
	sprintf(cmd, "%s %s %s", grep, filename, rest);
	if( (fp = popen(cmd, "r")) == NULL ) report_exit("popen error!\n");

	CKT_TOP * p_last_circuit=NULL;
	pair <int, Circuit*> layer_ckt_pair;
	// now read filename.info to create circuits (they are SORTED)
	while( fscanf(fp, "%s %d %s", name, &layer, GL) != EOF ){
		string name_string(name);
		string GL_string(GL);
		//cout<<name_string<<" "<<layer<<" "<<GL_string<<endl;
		string fullname = name_string;
		// compare with previous circuit name 
		if( prev_ckt_name == "" ||
		    name_string != prev_ckt_name ){
			fullname.append("G");
			//Circuit * circuit = new Circuit(name_string);
			CKT_TOP * ckt_top = new CKT_TOP(name_string);
			Circuit * ckt1 = new Circuit(fullname);
			fullname = name_string;
			fullname.append("L");
			ckt1->type = "GLOBAL";
			Circuit * ckt2 = new Circuit(fullname);
			ckt_top->ckt1 = ckt1;
			ckt_top->ckt2 = ckt2;
			
			ckt2->type = "LOCAL";
			(*p_ckts).push_back(*ckt_top);
			++n_circuit;
			prev_ckt_name = name_string;
			p_last_circuit = ckt_top;
		}
		
		layer_ckt_pair.first = layer;
		if(GL_string == "G"){
			layer_ckt_pair.second = p_last_circuit->ckt1;
			p_last_circuit->ckt1->layers.push_back(layer);
		}
		else if(GL_string == "L"){
			layer_ckt_pair.second = p_last_circuit->ckt2;
			p_last_circuit->ckt2->layers.push_back(layer);
		}		
		layer_map_ckt.insert(layer_ckt_pair);
		this->n_layer++;
	}
	
	if( (status = pclose(fp)) == -1 )    report_exit("pclose error!\n");

	Circuit::layer_dir.resize(1);
	return n_circuit;
}

// parse the file
// Note: the file will be parsed twice
// the first time is to find the layer information
// and the second time is to create nodes
void Parser::parse(char * filename, Tran & tran){
	//filename = "../data/netlist_2M.txt";
	this->filename = filename;
	int count[6];
	for(size_t i=0;i<6;i++)
		count[i] = 0;
	//clog<<"open "<<filename<<endl;

	FILE * f;
	f = fopen(filename, "r");
	if( f == NULL ) 
		report_exit("Input file not exist!\n");
	// first time parse:
	create_circuits();
	
	// second time parser:
	char line[MAX_BUF];
	string l;
	while( fgets(line, MAX_BUF, f) != NULL ){
		char type = line[0];
		switch(type){
		case 'r': // resistor
		case 'R':
		case 'v': // VDD
		case 'V':
		case 'i': // current
		case 'I':
		case 'c':
		case 'C':
		case 'l':
		case 'L':
			insert_net_node(line, count);
			break;
		case '.': // command
			parse_dot(line, tran);	
		case '*': // comment
		case ' ':
		case '\n':
			break;
		default:
			printf("Unknown input line: ");
			report_exit(line);
			break;
		}
	}
	fclose(f);

	// release map_node resource
	/*for(size_t i=0;i<(*p_ckts).size();i++){
		Circuit * ckt = (*p_ckts)[i];
		ckt->map_node.clear();
	}*/
}// end of parse

int Parser::get_num_layers() const{ return n_layer; }

void Parser::parse_dot(char *line, Tran &tran){
	char *chs;
	char *saveptr;
	char sname[MAX_BUF];
	Node_TR_PRINT item;
	const char *sep = "= v() \n";
	switch(line[1]){
		case 't': // transient steps
			sscanf(line, "%s %lf %lf", sname, 
				&tran.step_t, &tran.tot_t);
			tran.isTran = 1; // do transient ana;
			//clog<<"step: "<<tran.step_t<<" tot: "<<tran.tot_t<<endl;
			break;
		case 'w': // output length
			chs = strtok_r(line, sep, &saveptr);
			chs = strtok_r(NULL, sep, &saveptr);
			chs = strtok_r(NULL, sep, &saveptr);
			tran.length = atoi(chs);
			//clog<<"out len: "<<tran.length<<endl;
			break;
		case 'p': // print
			chs = strtok_r(line, sep, &saveptr);
			chs = strtok_r(NULL, sep, &saveptr);
			while(chs != NULL){
				chs = strtok_r(NULL, sep, &saveptr);
				if(chs == NULL) break;
				item.name = chs;
				//item.node = get_node(chs);
				tran.nodes.push_back(item);
			};
			break;
		default: 
			break;
	}
}

void Parser::add_node(Circuit *ckt, Node &nd, Node *nd_ptr){
	if ( (nd_ptr = ckt->get_node(nd.name) ) == NULL ){
		// create new node and insert
		nd_ptr = new Node(nd); // copy constructor
		nd_ptr->rep = nd_ptr;  // set rep to be itself
		ckt->add_node(nd_ptr);
		if( nd_ptr->isS()== X)    // determine circuit type
			ckt->set_type(WB);
		// find the coordinate max and min	
	}	
}

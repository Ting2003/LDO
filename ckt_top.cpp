#include "ckt_top.h"

CKT_TOP::CKT_TOP(string _name){
	ckt1 = NULL;
	ckt2 = NULL;
	name = _name;
}

CKT_TOP::~CKT_TOP(){
	boundary_net.clear();
	ckt1->release_ckt_nodes();
	ckt2->release_ckt_nodes();
}

// link global and local interface nodes
void CKT_TOP::link_ckt1(){
	Node_TR_PRINT node_temp;
	Node *na, *nb;
	Net *net;
	Node *nd_g;
	// get global ckt_nodes
	for(size_t k=0;k<ckt2->ckt_nodes.size();k++){
		Node *nd_target = ckt2->ckt_nodes[k].node;
		//clog<<"nd_local: "<<*nd_target<<endl;
		for(size_t i=0;i<boundary_net.size();i++){
			net = boundary_net[i];
			if(net->ab[0]->name == nd_target->name || net->ab[1]->name == nd_target->name){
				for(int j=0;j<2;j++){
					na = net->ab[j];
					nd_g = ckt1->get_node(na->name);
					if(nd_g == NULL) continue;
					node_temp.node = nd_g;
					//clog<<"nd_global: "<<*nd_g<<endl;
					ckt1->ckt_nodes.push_back(node_temp);
				}
				break;
			}
		}	
	}
}

// assign the local current into global circuit
void CKT_TOP::assign_gcur(){
	Node *nd;
	int type = CURRENT;
	Node *nd_new;
	Net *net;
	Node *nd_global;
	// first clear all the transient values
	for(size_t i=0;i<ckt1->ckt_nodes.size();i++){
		ckt1->ckt_nodes[i].value.clear();
		ckt1->ckt_nodes[i].value_cur.clear();
	}

	for(size_t i=0;i<ckt2->pad_set.size();i++){
		double current_nd = 0;
		nd = ckt2->pad_set[i]->node;
		//clog<<"pad: "<<*nd<<endl;
		net = nd->nbr[BOTTOM];
		if(net == NULL) continue;
		nd_new = net->ab[0];
		if(nd_new->name == nd->name)
			nd_new = net->ab[1];
		//clog<<"nd_new: "<<*nd_new<<endl;
		// get nd_new->value from tran
		for(size_t j=0;j<ckt2->ckt_nodes.size();j++){
			Node *nd_temp = ckt2->ckt_nodes[j].node;	
			// locate this node
			if(nd_temp->name != nd_new->name)
				continue;
			//nd_global = ckt1->ckt_nodes[j].node;
			//clog<<"nd_new, nd_global: "<<*nd_new<<" "<<*nd_global<<endl;
			vector<double> vol = ckt2->ckt_nodes[j].value;
			for(size_t k=0;k<vol.size();k++){
				current_nd = (nd->value - vol[k])/net->value;
				ckt1->ckt_nodes[j].value_cur.push_back(current_nd);
			}
			//for(size_t k=0;k<ckt1->ckt_nodes[j].value.size();k++)
				//clog<<"nd, k, value: "<<*nd_global<<" "<<ckt1->ckt_nodes[j].value[k]<<endl;
			break;
		}
	}
}

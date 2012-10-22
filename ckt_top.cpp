#include "ckt_top.h"

CKT_TOP::CKT_TOP(string _name){
	ckt1 = NULL;
	ckt2 = NULL;
	name = _name;
}

CKT_TOP::~CKT_TOP(){
	boundary_net.clear();
}

/*
// assign the local current into global circuit
CKT_TOP::assign_g_cur(Tran &tran){
	Circuit *ckt_ptr;
	if(ckt1->type == "GLOBAL")
		ckt_ptr = ckt1;
	else
		ckt_ptr = ckt2;
	int type = CURRENT;
	for(size_t i=0;i<tran.nodes.size();i++){
		
	}
}*/

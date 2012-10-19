// ----------------------------------------------------------------//
// Filename : block.h
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// declaration of block class
// ----------------------------------------------------------------//

#ifndef __BLOCK_H__
#define __BLOCK_H__
#include <fstream>
#include "triplet.h"
#include "global.h"
#include "vec.h"
#include "net.h"
#include "util.h"
#include "cholmod.h"
using namespace std;

class Block{
private:
	vector<Node *> nodelist;
        vector<Net *> netlist;
	// block id;
        size_t id;
public:
        Block();
        ~Block();
        void free_block_cholmod(cholmod_common *cm);
        void CK_decomp(Matrix & A, cholmod_common *cm);
        void LU_decomposition();
        void solve_CK(cholmod_common *cm);	// solve for the current matrix

        // get column compressed form from a Matrix
        void get_col_compressed(Matrix & A);

        // allocate space for the matrix and vectors accoding to count
        void allocate_resource(cholmod_common *cm);

        void update_rhs(int flag_tr);
 
        cholmod_factor * L;


        // vector b
        cholmod_dense * b_ck, *b_tr, *b_new_ck;
        // pointer to b_ck, b_new_ck, and x_ck;
        double *bp, *b_trp, *bnewp, *xp, *x_old;
        // solution
        cholmod_dense *x_ck;
};

class BlockInfo{
public:
	friend class Circuit;	
private:
	Block block_g;
	Block block_l;
};

#endif

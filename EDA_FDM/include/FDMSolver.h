/*
 * FDMSolver.h
 *
 *  Created on: 8 Jul, 2016
 *      Author: ming
 */

#ifndef FDMSOLVER_H_
#define FDMSOLVER_H_

#include <vector>
#include "umfpack.h"
#include "cs.h"
#include <assert.h>
#include "matrix.h"
#include "Configure.h"
#include "Debugger.h"

class FDM_Solver {
	std::string ans_path = "answer/C.txt";
	Configure macro_conf;
	int inum;	// inner node number
	int bnum;	// boundary node number
	double uscale[3];    // 1 / the unit scale of each axis
	double uscale2[3];    // the square of ( 1 / unit scale of each axis)
	int*** points;    // store the grid points, and he value is the index of the point in A. If the value == -1, then it means it's not a valid point, such as the 8 vertex of cube
                      // First dim -> axis1, Second dim -> axis2, Third dim -> axis3.
	Matrix_Sparse A11;
	Matrix_Sparse A12;
	Matrix_Sparse A21;
	Matrix_Sparse A22;
	cs_di Cap;

	Debugger debg;

public:
	FDM_Solver();
	FDM_Solver(Configure mconf, Matrix_Sparse a11, Matrix_Sparse a12, Matrix_Sparse a21, Matrix_Sparse a22);
	void Construct_Matrix_A11();
	void Construct_Matrix_A12();
	void Construct_Matrix_A21();
	void Construct_Matrix_A22();
	void solve();
	void write();
	virtual ~FDM_Solver();
	void Debug();
};

#endif /* FDMSOLVER_H_ */

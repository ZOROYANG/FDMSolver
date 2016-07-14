/*
 * Matrix.h
 *
 *  Created on: 8 Jul, 2016
 *      Author: ming
 */

#ifndef _MATRIX_H
#define _MATRIX_H

#include "Debugger.h"
#include <stddef.h>
#include <vector>

class Matrix_Sparse{
	Debugger debg;
public:
	size_t n_row; //the row of the matrix
	size_t n_col; //the column of the matrix
	int nz;  //the number of nonzero elements in  the matrix
	//int status; //"umfpack_di_triplet_to_col" return value;
	/*........................*/
	//    the triplet form    //
	/*........................*/
	std::vector<int> i; //the row index
	std::vector<int> j; //the column index
	std::vector<double> x; // numerical value
	/*........................*/
	//    the CSC form        //
	/*........................*/
	std::vector<int> Ap; // column pointer, size is n_col + 1
	std::vector<int> Ai; //row index
	std::vector<double> Ax; //numerical value

	Matrix_Sparse();
	Matrix_Sparse(int r,int c);
	void matrixFree();
	void Debug(std::string s);
};

#endif

#include "matrix.h"

Matrix_Sparse::Matrix_Sparse() : n_row(0), n_col(0), nz(0){

}

Matrix_Sparse::Matrix_Sparse(int r, int c) : n_row(r), n_col(c), nz(0){

}

void Matrix_Sparse::matrixFree(){
	// construct a local variable which is a null vector, then swap it and the one that we want to clear. In this way, we can collect the memory without destroy the object.
	std::vector<int>().swap(this->i);
	std::vector<int>().swap(this->j);
	std::vector<double>().swap(this->x);
	std::vector<int>().swap(this->Ap);
	std::vector<int>().swap(this->Ai);
	std::vector<double>().swap(this->Ax);
	this->n_row = 0 ;
	this->n_col = 0;
	this->nz = 0 ;
	return;
}

void Matrix_Sparse::Debug(std::string s){
	this->debg.Matrix_print(s, this->n_row, this->n_col, this->nz, this->i, this->j, this->x, this->Ai, this->Ap, this->Ax);
	return;
}

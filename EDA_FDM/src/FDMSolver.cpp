/*
 * FDMSolver.cpp
 *
 *  Created on: 8 Jul, 2016
 *      Author: ming
 */

#include "FDMSolver.h"
FDM_Solver::FDM_Solver() : inum(0), bnum(0){
	// TODO Auto-generated constructor stub

}

FDM_Solver::FDM_Solver(Configure mconf, Matrix_Sparse a11, Matrix_Sparse a12, Matrix_Sparse a21, Matrix_Sparse a22) : macro_conf(mconf), A11(a11), A12(a12), A21(a21), A22(a22){
	// TODO Auto-generated constructor stub
	int* sa = macro_conf.getSegAxis();
	double* len = macro_conf.getLength();
	inum = sa[0] * sa[1] * sa[2];
	bnum = 2 * (sa[0] * sa[1] + sa[0] * sa[2] + sa[2] * sa[1]);

	for(int i = 0; i < 3; ++ i){
		uscale[i] = double(sa[i]) * 2. / len[i];
		uscale2[i] = uscale[i] * uscale[i];
	}

	for(int i = 0; i < 3; ++ i){
			sa[i] += 2;
	}


	points = new int**[sa[0]];
	for(int i = 0; i < sa[0]; ++ i){
		points[i] = new int*[sa[1]];
		for(int j = 0; j < sa[1]; ++ j){
			points[i][j] = new int[sa[2]];
		}
	}

	// assign the index of every point in matrix A
	int bidx = 0;	// boundary points index in A
	int iidx = 0;	// inner points index in A
	for(int i = 0; i < sa[0]; ++ i){
		for(int j = 0; j < sa[1]; ++ j){
			for(int k = 0; k < sa[2]; ++ k){
				if((i == 0 && j == 0) || (i == 0 && k == 0) || (k == 0 && j == 0)){
					points[i][j][k] = -1;
					continue;
				}

				if((i == (sa[0] - 1) && j == (sa[1] - 1)) || (i == (sa[0] - 1) && k == 0) || (k == 0 && j == (sa[1] - 1))){
					points[i][j][k] = -1;
					continue;
				}
				if((i == (sa[0] - 1) && j == 0) || (i == (sa[0] - 1) && k == (sa[2] - 1)) || (k == (sa[2] - 1) && j == 0)){
					points[i][j][k] = -1;
					continue;
				}
				if((i == 0 && j == (sa[1] - 1)) || (i == 0 && k == (sa[2] - 1)) || (k == (sa[2] - 1) && j == (sa[1] - 1))){
					points[i][j][k] = -1;
					continue;
				}

				if(i == 0 || i == (sa[0] - 1) || j == 0 || j == (sa[1] - 1) || k == 0 || k == (sa[2] - 1)){
					points[i][j][k] = inum + bidx;
					++ bidx;
					continue;
				}

				points[i][j][k] = iidx;
				++ iidx;
			}
		}
	}

	this->Construct_Matrix_A11();
	this->A11.Debug("A11.txt");

	this->Construct_Matrix_A12();
	this->A12.Debug("A12.txt");

	this->Construct_Matrix_A21();
	this->A21.Debug("A21.txt");

	this->Construct_Matrix_A22();
	this->A22.Debug("A22.txt");

	delete[] sa;
	delete[] len;
}

FDM_Solver::~FDM_Solver() {
	// TODO Auto-generated destructor stub

	int* sa = macro_conf.getSegAxis();
	for(int i = 0; i < 3; ++ i){
		sa[i] += 2;
	}

	for(int i = 0; i < sa[0]; ++ i){
		for(int j = 0; j < sa[1]; ++ j){
			delete[] points[i][j];
		}
		delete[] points[i];
	}
	delete[] points;

	delete[] sa;

	this->A11.matrixFree();
	this->A12.matrixFree();
	this->A21.matrixFree();
	this->A22.matrixFree();
/*
	delete[] Cap.i;
	delete[] Cap.p;
	delete[] Cap.x;
*/
}

void FDM_Solver::Construct_Matrix_A11(){
	int* sa = macro_conf.getSegAxis();
	for(int i = 0; i < 3; ++ i){
		sa[i] += 2;
	}
	A11.n_row = inum;
    A11.n_col = inum;

    for(int i = 0; i < sa[0]; ++ i){
        for(int j = 0; j < sa[1]; ++ j){
            for(int k = 0; k < sa[2]; ++ k){
                if(points[i][j][k] < 0){
                    continue;
                }
                
                // boundary points
                if(points[i][j][k] >= inum){
                    
                    continue;
                }
				
                // inner points

                double center = 0.;
                if((i + 1) == (sa[0] - 1)){
                    if(i == 1){
                        center += -2. * this->uscale2[0];
                    }else{
                        center += -1. * this->uscale2[0];
                        
                        A11.i.push_back(points[i][j][k]);
                        A11.j.push_back(points[i - 1][j][k]);
                        A11.x.push_back(this->uscale2[0] / 3.);
                    }
                }else{
                    if(i == 1){
                        center += -1. * this->uscale2[0];
                        
                        A11.i.push_back(points[i][j][k]);
                        A11.j.push_back(points[i + 1][j][k]);
                        A11.x.push_back(this->uscale2[0] / 3.);
                    }else{
                        center += this->uscale2[0] / -2.;
                    
                        A11.i.push_back(points[i][j][k]);
                        A11.j.push_back(points[i + 1][j][k]);
                        A11.x.push_back(this->uscale2[0] / 4.);
                    
                        A11.i.push_back(points[i][j][k]);
                        A11.j.push_back(points[i - 1][j][k]);
                        A11.x.push_back(this->uscale2[0] / 4.);
                    }
                }

                if((j + 1) == (sa[1] - 1)){
                    if(j == 1){
                        center += -2. * this->uscale2[1];
                    }else{
                        center += -1. * this->uscale2[1];
                        
                        A11.i.push_back(points[i][j][k]);
                        A11.j.push_back(points[i][j - 1][k]);
                        A11.x.push_back(this->uscale2[1] / 3.);
                    }
                }else{
                    if(j == 1){
                        center += -1. * this->uscale2[1];
                        
                        A11.i.push_back(points[i][j][k]);
                        A11.j.push_back(points[i][j + 1][k]);
                        A11.x.push_back(this->uscale2[1] / 3.);
                    }else{
                        center += this->uscale2[1] / -2.;
                    
                        A11.i.push_back(points[i][j][k]);
                        A11.j.push_back(points[i][j + 1][k]);
                        A11.x.push_back(this->uscale2[1] / 4.);
                    
                        A11.i.push_back(points[i][j][k]);
                        A11.j.push_back(points[i][j - 1][k]);
                        A11.x.push_back(this->uscale2[1] / 4.);
                    }
                }
                
                if((k + 1) == (sa[2] - 1)){
                    if(k == 1){
                        center += -2. * this->uscale2[2];
                    }else{
                        center += -1. * this->uscale2[2];
                        
                        A11.i.push_back(points[i][j][k]);
                        A11.j.push_back(points[i][j][k - 1]);
                        A11.x.push_back(this->uscale2[2] / 3.);
                    }
                }else{
                    if(k == 1){
                        center += -1. * this->uscale2[2];
                        
                        A11.i.push_back(points[i][j][k]);
                        A11.j.push_back(points[i][j][k + 1]);
                        A11.x.push_back(this->uscale2[2] / 3.);
                    }else{
                        center += this->uscale2[2] / -2.;
                    
                        A11.i.push_back(points[i][j][k]);
                        A11.j.push_back(points[i][j][k + 1]);
                        A11.x.push_back(this->uscale2[2] / 4.);
                    
                        A11.i.push_back(points[i][j][k]);
                        A11.j.push_back(points[i][j][k - 1]);
                        A11.x.push_back(this->uscale2[2] / 4.);
                    }
                }

                A11.i.push_back(points[i][j][k]);
                A11.j.push_back(points[i][j][k]);
                A11.x.push_back(center);

            }
        }
    }


    int status;
	A11.nz = A11.i.size();
	A11.Ap.resize(A11.n_col + 1, 0);
	A11.Ai.resize(A11.nz, 0);
	A11.Ax.resize(A11.nz, 0);
	status = umfpack_di_triplet_to_col(A11.n_row, A11.n_col, A11.nz, &A11.i[0], &A11.j[0], &A11.x[0], &A11.Ap[0], &A11.Ai[0], &A11.Ax[0], NULL);
	if(status != 0){
		printf("A11: triplet_to_col goes wrong!\n");
	}

    delete[] sa;
    return;
}

void FDM_Solver::Construct_Matrix_A12(){
	int* sa = macro_conf.getSegAxis();
	for(int i = 0; i < 3; ++ i){
		sa[i] += 2;
	}
	A12.n_row = inum;
	A12.n_col = bnum;

	for(int i = 0; i < sa[0]; ++ i){
		for(int j = 0; j < sa[1]; ++ j){
			for(int k = 0; k < sa[2]; ++ k){
				if(points[i][j][k] < 0){
					continue;
				}

				// boundary points
				if(points[i][j][k] >= inum){
					continue;
				}

				// inner points
				if((i + 1) == (sa[0] - 1)){
					if(i == 1){
						A12.i.push_back(points[i][j][k]);
						A12.j.push_back(points[i + 1][j][k] - inum);
						A12.x.push_back(this->uscale2[0]);

						A12.i.push_back(points[i][j][k]);
						A12.j.push_back(points[i - 1][j][k] - inum);
						A12.x.push_back(this->uscale2[0]);
					}else{
						A12.i.push_back(points[i][j][k]);
						A12.j.push_back(points[i + 1][j][k] - inum);
						A12.x.push_back(2. * this->uscale2[0] / 3.);
					}
				}else{
					if(i == 1){
						A12.i.push_back(points[i][j][k]);
						A12.j.push_back(points[i - 1][j][k] - inum);
						A12.x.push_back(2. * this->uscale2[0] / 3.);
					}else{
					}
				}

				if((j + 1) == (sa[1] - 1)){
					if(j == 1){
						A12.i.push_back(points[i][j][k]);
						A12.j.push_back(points[i][j + 1][k] - inum);
						A12.x.push_back(this->uscale2[1]);

						A12.i.push_back(points[i][j][k]);
						A12.j.push_back(points[i][j - 1][k] - inum);
						A12.x.push_back(this->uscale2[1]);
					}else{
						A12.i.push_back(points[i][j][k]);
						A12.j.push_back(points[i][j + 1][k] - inum);
						A12.x.push_back(2. * this->uscale2[1] / 3.);
					}
				}else{
					if(j == 1){
						A12.i.push_back(points[i][j][k]);
						A12.j.push_back(points[i][j - 1][k] - inum);
						A12.x.push_back(2. * this->uscale2[1] / 3.);
					}else{
					}
				}

				if((k + 1) == (sa[2] - 1)){
					if(k == 1){
						A12.i.push_back(points[i][j][k]);
						A12.j.push_back(points[i][j][k + 1] - inum);
						A12.x.push_back(this->uscale2[2]);

						A12.i.push_back(points[i][j][k]);
						A12.j.push_back(points[i][j][k - 1] - inum);
						A12.x.push_back(this->uscale2[2]);
					}else{
						A12.i.push_back(points[i][j][k]);
						A12.j.push_back(points[i][j][k + 1] - inum);
						A12.x.push_back(2. * this->uscale2[2] / 3.);
					}
				}else{
					if(k == 1){
						A12.i.push_back(points[i][j][k]);
						A12.j.push_back(points[i][j][k - 1] - inum);
						A12.x.push_back(2. * this->uscale2[2] / 3.);
					}else{
					}
				}
			}
		}
	}

	int status;
	A12.nz = A12.i.size();
	A12.Ap.resize(A12.n_col + 1, 0);
	A12.Ai.resize(A12.nz, 0);
	A12.Ax.resize(A12.nz, 0);
	status = umfpack_di_triplet_to_col(A12.n_row, A12.n_col, A12.nz, &A12.i[0], &A12.j[0], &A12.x[0], &A12.Ap[0], &A12.Ai[0], &A12.Ax[0], NULL);
	if(status != 0){
		printf("A12: triplet_to_col goes wrong!\n");
	}

	delete[] sa;
	return;
}

void FDM_Solver::Construct_Matrix_A21(){
	int* sa = macro_conf.getSegAxis();
	for(int i = 0; i < 3; ++ i){
		sa[i] += 2;
	}
	A21.n_row = bnum;
	A21.n_col = inum;

	for(int i = 0; i < sa[0]; ++ i){
		for(int j = 0; j < sa[1]; ++ j){
			for(int k = 0; k < sa[2]; ++ k){
				if(points[i][j][k] < 0){
					continue;
				}

				// boundary points, take care du/dx = -E
				if(points[i][j][k] >= inum){
					A21.i.push_back(points[i][j][k] - inum);
					if(i == 0){
						A21.j.push_back(points[i + 1][j][k]);
						A21.x.push_back(-1. * this->uscale[0]);
						continue;
					}

					if(i == (sa[0] - 1)){
						A21.j.push_back(points[i - 1][j][k]);
						A21.x.push_back(-1. * this->uscale[0]);
						continue;
					}

					if(j == 0){
						A21.j.push_back(points[i][j + 1][k]);
						A21.x.push_back(-1. * this->uscale[1]);
						continue;
					}

					if(j == (sa[1] - 1)){
						A21.j.push_back(points[i][j - 1][k]);
						A21.x.push_back(-1. * this->uscale[1]);
						continue;
					}

					if(k == 0){
						A21.j.push_back(points[i][j][k + 1]);
						A21.x.push_back(-1. * this->uscale[2]);
						continue;
					}

					if(k == (sa[2] - 1)){
						A21.j.push_back(points[i][j][k - 1]);
						A21.x.push_back(-1. * this->uscale[2]);
						continue;
					}
				}
			}
		}
	}

	int status;
	A21.nz = A21.i.size();
	A21.Ap.resize(A21.n_col + 1, 0);
	A21.Ai.resize(A21.nz, 0);
	A21.Ax.resize(A21.nz, 0);
	status = umfpack_di_triplet_to_col(A21.n_row, A21.n_col, A21.nz, &A21.i[0], &A21.j[0], &A21.x[0], &A21.Ap[0], &A21.Ai[0], &A21.Ax[0], NULL);
	if(status != 0){
		printf("A21: triplet_to_col goes wrong!\n");
	}

	delete[] sa;
	return;
}

void FDM_Solver::Construct_Matrix_A22(){
	int* sa = macro_conf.getSegAxis();
	for(int i = 0; i < 3; ++ i){
		sa[i] += 2;
	}
	A22.n_row = bnum;
	A22.n_col = bnum;

	for(int i = 0; i < sa[0]; ++ i){
		for(int j = 0; j < sa[1]; ++ j){
			for(int k = 0; k < sa[2]; ++ k){
				if(points[i][j][k] < 0){
					continue;
				}

				// boundary points, take care du/dx = -E
				if(points[i][j][k] >= inum){
					A22.i.push_back(points[i][j][k] - inum);
					A22.j.push_back(points[i][j][k] - inum);
					if(i == 0){
						A22.x.push_back(this->uscale[0]);
						continue;
					}

					if(i == (sa[0] - 1)){
						A22.x.push_back(this->uscale[0]);
						continue;
					}

					if(j == 0){
						A22.x.push_back(this->uscale[1]);
						continue;
					}

					if(j == (sa[1] - 1)){
						A22.x.push_back(this->uscale[1]);
						continue;
					}

					if(k == 0){
						A22.x.push_back(this->uscale[2]);
						continue;
					}

					if(k == (sa[2] - 1)){
						A22.x.push_back(this->uscale[2]);
						continue;
					}
				}
			}
		}
	}

	int status;
	A22.nz = A22.i.size();
	A22.Ap.resize(A22.n_col + 1, 0);
	A22.Ai.resize(A22.nz, 0);
	A22.Ax.resize(A22.nz, 0);
	status = umfpack_di_triplet_to_col(A22.n_row, A22.n_col, A22.nz, &A22.i[0], &A22.j[0], &A22.x[0], &A22.Ap[0], &A22.Ai[0], &A22.Ax[0], NULL);
	if(status != 0){
		printf("A22: triplet_to_col goes wrong!\n");
	}

	delete[] sa;
	return;
}
 // dielectric value is diel, S may be a matrix
void FDM_Solver::solve(){
	// calculate A22 - A21 * inv(A11) * A12
    // construct a22
    cs_di a22;
    a22.m = A22.n_row;
    a22.n = A22.n_col;
    a22.p = &(A22.Ap[0]);
    a22.i = &(A22.Ai[0]);
    a22.x = &(A22.Ax[0]);
    a22.nzmax = A22.nz;
    a22.        nz = -1;
    // construct a21
    cs_di a21;
    a21.m = A21.n_row;
    a21.n = A21.n_col;
    a21.p = &(A21.Ap[0]);
    a21.i = &(A21.Ai[0]);
    a21.x = &(A21.Ax[0]);
    a21.nzmax = A21.nz;
    a21.nz = -1;
   // construct a11
    cs_di a11;
    a11.m = A11.n_row;
    a11.n = A11.n_col;
    a11.p = &(A11.Ap[0]);
    a11.i = &(A11.Ai[0]);
    a11.x = &(A11.Ax[0]);
    a11.nzmax = A11.nz;
    a11.nz = -1;
    // construct a12
    cs_di a12;
    a12.m = A12.n_row;
    a12.n = A12.n_col;
    a12.p = &(A12.Ap[0]);
    a12.i = &(A12.Ai[0]);
    a12.x = &(A12.Ax[0]);
    a12.nzmax = A12.nz;
    a12.nz = -1;
        
    cs_di pCap; // part Capacitency
    pCap.nz = -1;
    pCap.m = inum;
    pCap.n = bnum;
    int countnz = 0;

    //pCap.p = new int[pCap.n + 1];
    std::vector<int> p_temp;
    std::vector<int> row_idx;
    std::vector<double> value;

    /* Here is supposed to apply UMFPACK_A, because
        	 * Int sys ;		Input argument, not modified.

    	Defines which system to solve.  (') is the linear algebraic transpose
    	(complex conjugate if A is complex), and (.') is the array transpose.

    	    sys value	    system solved
    	    UMFPACK_A       Ax=b
    	    UMFPACK_At      A'x=b
    	    UMFPACK_Aat     A.'x=b
    	    UMFPACK_Pt_L    P'Lx=b
    	    UMFPACK_L       Lx=b
    	    UMFPACK_Lt_P    L'Px=b
    	    UMFPACK_Lat_P   L.'Px=b
    	    UMFPACK_Lt      L'x=b
    	    UMFPACK_U_Qt    UQ'x=b
    	    UMFPACK_U       Ux=b
    	    UMFPACK_Q_Ut    QU'x=b
    	    UMFPACK_Q_Uat   QU.'x=b
    	    UMFPACK_Ut      U'x=b
    	    UMFPACK_Uat     U.'x=b
    */
    for(int k = 0; k < a12.n; ++ k){
    	// x is as long as the width of a11, a.s. sa1 * sa2 * sa3
    	std::vector<double> x(pCap.m, 0);

    	std::vector<double> b(pCap.m, 0);
    	for(int j = a12.p[k]; j < a12.p[k + 1]; ++ j){
    		b[a12.i[j]] = a12.x[j];
    	}

    	void *Symbolic , *Numeric;
		double *null = (double *) NULL ;
		umfpack_di_symbolic(a11.m, a11.n, a11.p, a11.i, a11.x, &Symbolic, null, null);
		umfpack_di_numeric(a11.p, a11.i, a11.x, Symbolic, &Numeric, null, null);
		umfpack_di_free_symbolic(&Symbolic);
    	// reference to above
    	umfpack_di_solve(UMFPACK_A, a11.p, a11.i, a11.x, &x[0], &b[0], Numeric, null, null);
    	umfpack_di_free_numeric(&Numeric);

    	//pCap.p[k] = countnz;
    	p_temp.push_back(countnz);
    	for(int u = 0; u < pCap.m; ++ u){
    		if(x[u] != 0){
    			++ countnz;
    			row_idx.push_back(u);
    			value.push_back(x[u]);
    		}
    	}
    }
        
    //pCap.p[pCap.n] = countnz;
    p_temp.push_back(countnz);
    // assure it's valid
    assert(row_idx.size() == value.size());
        
    //pCap.i = new int[row_idx.size()];
    //pCap.x = new double[value.size()];
    //int vsize = value.size();
    //for(int k = 0; k < vsize; ++ k){
    	//pCap.i[k] = row_idx[k];
    	//pCap.x[k] = value[k];
    //}
    pCap.i = &(row_idx[0]);
    pCap.p = &(p_temp[0]);
    pCap.x = &(value[0]);

    pCap.nzmax = countnz;
        
    cs_di* res1 = cs_di_multiply(&a21, &pCap);
    cs_di* res2 = cs_di_add(&a22, res1, 1, -1);

    // need to delete the mem of Cap

    Cap.m = res2->m;
    Cap.n = res2->n;
    Cap.nzmax = res2->nzmax;
    Cap.p = new int[Cap.n + 1];
    for(int k = 0; k <= Cap.n; ++ k){
        Cap.p[k] = res2->p[k];
    }

    Cap.i = new int[Cap.nzmax];
    for(int k = 0; k < Cap.nzmax; ++ k){
        Cap.i[k] = res2->i[k];
    }
    Cap.x = new double[Cap.nzmax];
    for(int k = 0; k < Cap.nzmax; ++ k){
    	Cap.x[k] = res2->x[k];
    }
    Cap.nz = -1;
        
    cs_di_spfree(res1);
    cs_di_spfree(res2);

    return;
}

void FDM_Solver::write(){
	FILE* output;
	output = fopen(ans_path.c_str(), "w");
	if(output == NULL){
		printf("Cann't open %s\n", ans_path.c_str());
		exit(1);
	}

	/*fprintf(output, "rows:\n    %d\ncols:\n    %d\nnone zeros:\n    %d\n", Cap.m, Cap.n, Cap.nzmax);

	fprintf(output, "\nAi:\n");
	for(int u = 0; u < Cap.nzmax; ++ u){
		fprintf(output, "%d ", Cap.i[u]);
	}
	fprintf(output, "\n");

	fprintf(output, "\nAp:\n");
	for(int u = 0; u <= Cap.n; ++ u){
		fprintf(output, "%d ", Cap.p[u]);
	}
	fprintf(output, "\n");

	fprintf(output, "\nAx:\n");
	for(int u = 0; u < Cap.nzmax; ++ u){
		fprintf(output, "%2f ", Cap.x[u]);
	}
	fprintf(output, "\n");
*/
	// matrix
	double** A = new double*[Cap.m];
	for(int u = 0; u < Cap.m; ++ u){
		A[u] = new double[Cap.n];
	}

	for(int v = 0; v < Cap.n; ++ v){
		for(int u = 0; u < Cap.m; ++ u){
			A[u][v] = 0;
		}
	}

	for(int k = 0; k < Cap.n; ++ k){
		for(int j = Cap.p[k]; j < Cap.p[k + 1]; ++ j){
			A[Cap.i[j]][k] = Cap.x[j];
		}
	}

	fprintf(output, "\n");
	for(int u = 0; u < Cap.m; ++ u){
		for(int v = 0; v < Cap.n; ++ v){
			fprintf(output, "%2f ", A[u][v]);
		}
		fprintf(output, ";\n");
	}

	fclose(output);
	for(int u = 0; u < Cap.m; ++ u){
		delete[] A[u];
	}
	delete[] A;
	return;
}

void FDM_Solver::Debug(){
	int* sa = macro_conf.getSegAxis();
	for(int i = 0; i < 3; ++ i){
		sa[i] += 2;
	}
	this->debg.Order_print(this->points, sa[0], sa[1], sa[2], this->inum, this->bnum, this->uscale, this->uscale2);
	delete[] sa;
	return;
}

/*
 * Debugger.cpp
 *
 *  Created on: 12 Jul, 2016
 *      Author: ming
 */

#include "Debugger.h"

Debugger::Debugger(){

}

Debugger::~Debugger(){

}

void Debugger::Conf_print(double* ori, double* len, int* sa, double d){
	FILE*  confs;
	confs = fopen(this->fconf.c_str(), "w");
	if(confs == NULL){
		printf("Cann't open %s\n", this->fconf.c_str());
		exit(1);
	}
	fprintf(confs, "origin point:\n    (%f, %f, %f)\n", ori[0], ori[1], ori[2]);
	fprintf(confs, "length in axis1, axis2, axis3:\n    (%f, %f, %f)\n", len[0], len[1], len[2]);
	fprintf(confs, "segmentation in axis1, axis2, axis3\n    (%d, %d, %d)\n", sa[0], sa[1], sa[2]);
	fprintf(confs, "dielectric: %f\n", d);
	fclose(confs);
	return;
}

void Debugger::Order_print(int ***p, int x, int y, int z, int inum, int bnum, double* h, double* h2){
	FILE*  ords;
	ords = fopen(this->ford.c_str(), "w");
	if(ords == NULL){
		printf("Cann't open %s\n", this->ford.c_str());
		exit(1);
	}

	fprintf(ords, "number of inner points: %d\nnumber of boundary points: %d\n", inum, bnum);
	fprintf(ords, "the reciprocal of minimum unit: (%2f, %2f, %2f)", h[0], h[1], h[2]);
	fprintf(ords, "the reciprocal of the square of minimum unit: (%2f, %2f, %2f)", h2[0], h2[1], h2[2]);

	fprintf(ords, "Order of every point in cubic.\n");
	for(int i = 0; i < x; ++ i){
		for(int j = 0; j < y; ++ j){
			for(int k = 0; k < z; ++ k){
				fprintf(ords, "Point at (%d, %d, %d) is encoded as: %d\n", i, j, k, p[i][j][k]);
			}
		}
	}

	fclose(ords);
	return;
}

void Debugger::Matrix_print(std::string s, int r, int c, int nz, std::vector<int> i, std::vector<int> j, std::vector<double> x, std::vector<int> Ai, std::vector<int> Ap, std::vector<double> Ax){
	FILE*  as;
	as = fopen((this->fa + s).c_str(), "w");
	if(as == NULL){
		printf("Cann't open %s\n", (this->fa + s).c_str());
		exit(1);
	}
/*
	fprintf(as, "rows:\n    %d\ncols:\n    %d\nnone zeros:\n    %d\n", r, c, nz);
	fprintf(as, "i.size():\n    %lu\nj.size():\n    %lu\nx.size():\n    %lu\n", i.size(), j.size(), x.size());
	fprintf(as, "Ai.size():\n    %lu\nAp.size():\n    %lu\nAx.size():\n    %lu\n", Ai.size(), Ap.size(), Ax.size());

	fprintf(as, "\ni:\n");
	for(int u = 0; u < nz; ++ u){
		fprintf(as, "%d ", i[u]);
	}
	fprintf(as, "\n");

	fprintf(as, "\nj:\n");
	for(int u = 0; u < nz; ++ u){
		fprintf(as, "%d ", j[u]);
	}
	fprintf(as, "\n");

	fprintf(as, "\nx:\n");
	for(int u = 0; u < nz; ++ u){
		fprintf(as, "%f ", x[u]);
	}
	fprintf(as, "\n");
*/
	// matrix
	printf("r: %d, c: %d\n", r, c);

	double** A = new double*[r];
	for(int u = 0; u < r; ++ u){
		A[u] = new double[c];
	}

	for(int u = 0; u < r; ++ u){
		for(int v = 0; v < c; ++ v){
			A[u][v] = 0;
		}
	}

	for(int u = 0; u < nz; ++ u){
		A[i[u]][j[u]] = x[u];
	}

	fprintf(as, "\n");
	for(int u = 0; u < r; ++ u){
		for(int v = 0; v < c; ++ v){
			fprintf(as, "%f ", A[u][v]);
		}
		fprintf(as, ";\n");
	}

	/*
	fprintf(as, "\nAi:\n");
	for(int u = 0; u < nz; ++ u){
		fprintf(as, "%d ", Ai[u]);
	}
	fprintf(as, "\n");

	fprintf(as, "\nAp:\n");
	for(int u = 0; u <= c; ++ u){
		fprintf(as, "%d ", Ap[u]);
	}
	fprintf(as, "\n");

	fprintf(as, "\nAx:\n");
	for(int u = 0; u < nz; ++ u){
		fprintf(as, "%f ", Ax[u]);
	}
	fprintf(as, "\n");
	*/

	for(int u = 0; u < r; ++ u){
		delete[] A[u];
	}
	delete[] A;

	fclose(as);
	return;
}

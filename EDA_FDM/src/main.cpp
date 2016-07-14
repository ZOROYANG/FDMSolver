#include <iostream>
#include "matrix.h"
#include "Configure.h"
#include "FDMSolver.h"

// origin point
double ori[3] = {0., 0., 0.};

// length of 3 edges
double len[3] = {0.8, 0.8, 0.4};

// segmentation of 3 edges
int sa[3] = {1, 1, 1};
double diele = 1.;

int main(){
    Configure macro_conf(ori, len, sa, diele);
    macro_conf.Debug();
    Matrix_Sparse A11, A12, A21, A22;
    FDM_Solver fdms(macro_conf, A11, A12, A21, A22);
    fdms.Debug();
    fdms.solve();
    fdms.write();
    return 0;
}

#ifndef _FDMSOLVER_H_
#define _FDMSOLVER_H_

#include <iostream>
#include <fstream>

class FDMsolver{
private:
	Geometry ge;
	int nx, ny, nz;
	double epsilon;
	std::string filename;
	ofstream fout;
public:
	FDMSOLVER(Geometry gem, int x, int y, int z, double eps, std::string fn): ge(gem), nx(x), ny(y), nz(z), epsilon(eps), filename(fn){
		fout.open(fn);
	}
	void Calculator();
	~FDMSOLVER(){
		fout.close();
	}
}

#endif // _FDMSOLVER_H_
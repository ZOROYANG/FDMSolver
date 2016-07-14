/*
 * Debugger.h
 *
 *  Created on: 12 Jul, 2016
 *      Author: ming
 */

#ifndef DEBUGGER_H_
#define DEBUGGER_H_

#include <stdio.h>
#include <string>
#include <vector>

class Debugger{
	std::string fconf = "debugfiles/confs.txt";
	std::string ford = "debugfiles/ords.txt";
	std::string fa = "debugfiles/";
public:
	Debugger();
	~Debugger();
	void Conf_print(double* ori, double* len, int* sa, double d);
	void Order_print(int ***p, int x, int y, int z, int inum, int bnum, double* h, double* h2);
	void Matrix_print(std::string s,int r, int c, int nz, std::vector<int> i, std::vector<int> j, std::vector<double> x, std::vector<int> Ai, std::vector<int> Aj, std::vector<double> Ax);
};


#endif /* DEBUGGER_H_ */

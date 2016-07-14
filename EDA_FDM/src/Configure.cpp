/*
 * Configure.cpp
 *
 *  Created on: 8 Jul, 2016
 *      Author: ming
 */

#include "Configure.h"

Configure::Configure() {
	// TODO Auto-generated constructor stub
	for(int i = 0; i < 3; ++ i){
		this->origin[i] = 0.;
		this->length[i] = 0.;
		this->seg_axis[i] = 0;
	}
	this->diel = 1.;
}

Configure::Configure(double* ori, double* len, int* sa, double d){
	for(int i = 0; i < 3; ++ i){
			this->origin[i] = ori[i];
			this->length[i] = len[i];
			this->seg_axis[i] = sa[i];
		}
	this->diel = d;
}

Configure::~Configure() {
	// TODO Auto-generated destructor stub
}

double* Configure::getOp(){
	double* ori = new double[3];
	for(int i = 0; i < 3; ++ i){
		ori[i] = origin[i];
	}
	return ori;
}

double* Configure::getLength(){
	double* len = new double[3];
		for(int i = 0; i < 3; ++ i){
			len[i] = length[i];
		}
	return len;
}

int* Configure::getSegAxis(){
	int* sa = new int[3];
	for(int i = 0; i < 3; ++ i){
		sa[i] = seg_axis[i];
	}
	return sa;
}

void Configure::Debug(){
	this->debg.Conf_print(this->origin, this->length, this->seg_axis, this->diel);
	return;
}

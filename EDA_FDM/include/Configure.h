/*
 * Configure.h
 *
 *  Created on: 8 Jul, 2016
 *      Author: ming
 */

#ifndef CONFIGURE_H_
#define CONFIGURE_H_

#include "Debugger.h"
#include <stddef.h>

class Configure {
	double origin[3];  // Macro origin, origin[0]->axis1, origin[1]->axis2, origin[2]->axis3
	double length[3];	// Macro size, default for cubic, length[0]-> axis1, length[1]-> axis2, length[2]-> axis3
	int seg_axis[3];  // FDM segments in direction axis1, axis2, axis3(corresponding to axis[0], axis[1], axis3)
	double diel;	// dielectric value

	Debugger debg;
public:
	Configure();
	Configure(double* ori, double* len, int* sa, double d);
	// we copy to the call for get*, so don't forget to delete them
	double* getOp();
	double* getLength();
	int* getSegAxis();
	virtual ~Configure();
	void Debug();
};

#endif /* CONFIGURE_H_ */

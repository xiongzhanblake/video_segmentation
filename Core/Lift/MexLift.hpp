/*
 * MexLift.hpp
 *
 *  Created on: May 29, 2017
 *      Author: xiong
 */

#ifndef SRC_CORE_LIFT_MEXLIFT_HPP_
#define SRC_CORE_LIFT_MEXLIFT_HPP_

// Only Map Video pixel-space onto feature-space and down-sample to a grid;
void MexLift(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nlhs !=4 ) {printf("Output is not 2!\n");exit(-1);}
	if (nrhs !=6 ) {printf("Output is not 2!\n");exit(-1);}
}

#endif /* SRC_CORE_LIFT_MEXLIFT_HPP_ */

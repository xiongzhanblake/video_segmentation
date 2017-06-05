/*
 * FillIntersection.hpp
 *
 *  Created on: Feb 18, 2017
 *      Author: xiong
 */

#ifndef SRC_FILLINTERSECTION_HPP_
#define SRC_FILLINTERSECTION_HPP_

#include "../../Common/MapOperation.hpp"

//-------------------------------------------------------------------------------------------------
// Function MexIntersection: Set multiple columns using Map container.
// 1. Create a matrix, M, which has the same size of IndiceIn2;
// 2. Map DataIn1 onto M according to M, according to the intersection 
// between IndiceIn1 and IndiceIn2;
// 3. Set the rest values of M to zero;
//-------------------------------------------------------------------------------------------------
void MexIntersection( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
//-------------------------------------------------------------------------------------------------
// Check validation:
//-------------------------------------------------------------------------------------------------
	if (nrhs != 2)
		mexErrMsgTxt(Error_INPUT_NUMBER);
	if (nlhs != 1)
		mexErrMsgTxt(Error_OUTPUT_NUMBER);
	if (!mxIsCell(prhs[0]))
		mexErrMsgTxt("Mask is not a cell matrix!\n");
	if (mxGetN(prhs[0]) != 1 && mxGetM(prhs[0]) != 1)
		mexErrMsgTxt("Mask has no size as Mx1 or 1xN!\n");
//-------------------------------------------------------------------------------------------------
// Function input and out parameters:
//-------------------------------------------------------------------------------------------------
// Inputs:
	const mxArray *Mask = prhs[0];
	const mxArray *DataIndices = prhs[1];
// Outputs:
//	mxArray *IndiceOut;
	mxArray *DataOut;
//-------------------------------------------------------------------------------------------------
// Local variables:
//-------------------------------------------------------------------------------------------------
	mwSize rows = mxGetM(DataIndices);
	mwSize nClass = mxGetM(Mask)>mxGetN(Mask)?mxGetM(Mask):mxGetN(Mask);
	std::map<double, double>**Accu = new std::map<double, double>*[nClass];
	int inNum = 2;
	int outNum = 1;
	const mxArray **inArrays =new const mxArray*[inNum];
	mxArray **outArrays = new mxArray*[outNum];
//	mxArray *doubleScalar = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
//	double *doublescalar = (double *)mxGetData(doubleScalar);
//	mxArray *uint64Scalar = mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
//	uint64_t *uint64scalar = (uint64_t *)mxGetData(doubleScalar);
//	mxArray *tepInd = NULL;
//	mxArray *tepDat = NULL;
//-------------------------------------------------------------------------------------------------
//	Initialization:
//-------------------------------------------------------------------------------------------------
	for (mwSize n = 0; n < nClass; n++)
		Accu[n] = NULL;
	for (int i = 0; i < inNum; i ++)
		inArrays[i] = NULL;
	for (int i = 0; i < outNum; i ++)
		outArrays[i] = NULL;
	if (plhs[0] == NULL)
		plhs[0] = mxCreateNumericMatrix(rows,nClass,mxDOUBLE_CLASS,mxREAL);
//	if (plhs[1] == NULL)
//		plhs[1] = mxCreateNumericMatrix(rows,nClass,mxDOUBLE_CLASS,mxREAL);
//	IndiceOut = plhs[0];
	DataOut = plhs[0];
//-------------------------------------------------------------------------------------------------
//	Compute intersection and fill values:
//-------------------------------------------------------------------------------------------------
	for (mwSize n = 0; n < nClass; n++)
		inArrays[n] = mxGetCell(Mask,n);
	CreateMapByKeysAndValues<double,double>(nClass,inArrays,Accu);

	inArrays[0] = DataIndices; outArrays[0] = DataOut;
	AssignmentByIndicesMap<double,double,double>(1, outArrays, 1, inArrays, Accu);
	DataOut = outArrays[0]; outArrays[0] = NULL;

//	outArrays[0] = tepInd; outArrays[1] = tepDat;
//	SplitMap<uint64_t,double>(nlhs,plhs,1,Accu);
//	tepInd = outArrays[0]; tepDat = outArrays[1]; outArrays[0] = NULL; outArrays[1] = NULL;

//-------------------------------------------------------------------------------------------------
//	Free memories:
//-------------------------------------------------------------------------------------------------
	for (mwSize n = 0; n < nClass; n++)
		Accu[n]->clear();
	delete Accu;
	for (int i = 0; i < inNum; i ++)
		inArrays[i] = NULL;
	for (int i = 0; i < outNum; i ++)
		outArrays[i] = NULL;
	delete inArrays;
	delete outArrays;
//	mxDestroyArray(doubleScalar);
}
#endif /* SRC_FILLINTERSECTION_HPP_ */


/*
 * ComputeIndicesAndWeights.hpp
 *
 *  Created on: Feb 14, 2017
 *      Author: xiong
 */

#ifndef SRC_COMPUTEINDICESANDWEIGHTS_HPP_
#define SRC_COMPUTEINDICESANDWEIGHTS_HPP_

 #define SRC_COMPUTEINDICESANDWEIGHTS_HPP_DEBUG
// #define SRC_COMPUTEINDICESANDWEIGHTS_HPP_DETAIL

#include "../../Common/MexUtilities.hpp"
void ComputeIndicesAndWeights( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
//-------------------------------------------------------------------------------------------------
// Check validation:
//-------------------------------------------------------------------------------------------------
	if (nrhs != 5)
		mexErrMsgTxt(Error_INPUT_NUMBER);
	if (nlhs != 2)
		mexErrMsgTxt(Error_OUTPUT_NUMBER);
	const mxArray **Weights = new const mxArray *[2];;
//-------------------------------------------------------------------------------------------------
// Function input and out parameters:
//-------------------------------------------------------------------------------------------------
// Inputs:
	const mxArray *Index = prhs[0];
	const mxArray *Floors = prhs[1];
	const mxArray *GridSize = prhs[2];
	Weights[0] = prhs[3];
	Weights[1] = prhs[4];
//	const mxArray *NDims = prhs[5];
	
// Outputs:
	mxArray *Indices;
	mxArray *weight;
//-------------------------------------------------------------------------------------------------
// Local variables:
//-------------------------------------------------------------------------------------------------
	int inNum = 4;
	int outNum = 2;
	const mxArray **inArrays;
	mxArray **outArrays;
	mxArray *NDims;
	uint64_t *nDims;
	mxArray *temp = NULL;
	mxArray *binary = NULL;
	mxArray *uint64Scalar;
	uint64_t *uint64scalar;
	mxArray *doubleScalar;
	double *doublescalar;
	mxArray *NPoints;
	uint64_t *nPoints;
//-------------------------------------------------------------------------------------------------
// Initialization:
//-------------------------------------------------------------------------------------------------
	NDims = mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
	nDims = (uint64_t *)mxGetData(NDims);
	nDims[0] = mxGetN(Floors);
	NPoints = mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
	nPoints = (uint64_t *)mxGetData(NPoints);
	nPoints[0] = mxGetM(Floors);
	uint64Scalar = mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
	uint64scalar = (uint64_t *)mxGetData(uint64Scalar);
	doubleScalar = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
	doublescalar = (double *)mxGetData(doubleScalar);
	mwSize dimenstion = mxGetN(Weights[0]);


	inArrays = new const mxArray*[inNum];
	outArrays = new mxArray*[outNum];
	for (int i = 0; i < inNum; i ++)
		inArrays[i] = NULL;
	for (int i = 0; i < outNum; i ++)
		outArrays[i] = NULL;

	if (plhs[0] == NULL)
	{
		uint64scalar[0] = 0;
		inArrays[0]=NPoints; inArrays[1]=uint64Scalar;
		mexAllocateMemoryWithInitials<uint64_t>(1,outArrays,2,inArrays);
		plhs[0] = outArrays[0]; outArrays[0] = NULL;
	}
	if (plhs[1] == NULL)
	{
		doublescalar[0] = 1;
		inArrays[0]=NPoints; inArrays[1]=doubleScalar;
		mexAllocateMemoryWithInitials<double>(1,outArrays,2,inArrays);
		plhs[1] = outArrays[0]; outArrays[0] = NULL;
	}

	Indices = plhs[0];

//	#ifdef SRC_COMPUTEINDICESANDWEIGHTS_HPP_DEBUG
//	std::cout<<"Indices:\n"; SplitLine; inArrays[0] = Indices; mexDisplayArray<uint64_t>(0,NULL,1,inArrays);
//	#endif

	weight = plhs[1];

//	#ifdef SRC_COMPUTEINDICESANDWEIGHTS_HPP_DEBUG
//	std::cout<<"Weight:\n"; SplitLine; inArrays[0] = weight; mexDisplayArray<float>(0,NULL,1,inArrays);
//	#endif

//-------------------------------------------------------------------------------------------------
// Compute indices and weights according to binary string:
//-------------------------------------------------------------------------------------------------
	inArrays[0] = Index; inArrays[1] = NDims;
	mexDec2Bin<uint64_t>(1,outArrays,2,inArrays);
	binary = outArrays[0]; outArrays[0] = NULL;

	inArrays[0] = binary; outArrays[0] = binary;
	mexReverse<uint64_t>(1,outArrays,1,inArrays);

//	#ifdef SRC_COMPUTEINDICESANDWEIGHTS_HPP_DEBUG
//	std::cout<<"binary:\n"; SplitLine; inArrays[0] = binary; mexDisplayArray<uint64_t>(0,NULL,1,inArrays);
//	#endif

//	#ifdef SRC_COMPUTEINDICESANDWEIGHTS_HPP_DEBUG
//	inArrays[0] = Floors; inArrays[1] = binary;
//	MatrixAddition<float>(1,outArrays,2,inArrays);
//	std::cout<<"Indices Array:\n"; SplitLine; inArrays[0] = outArrays[0]; mexDisplayArray<float>(0,NULL,1,inArrays);
//	mxDestroyArray(outArrays[0]);
//	#endif

	inArrays[0] = Floors; inArrays[1] = binary; inArrays[2] = GridSize; outArrays[0] = Indices;
	mexComputeIndices(1, outArrays, 3, inArrays);
	Indices = outArrays[0]; outArrays[0] = NULL;

//	#ifdef SRC_COMPUTEINDICESANDWEIGHTS_HPP_DEBUG
//	std::cout<<"Indices:\n"; SplitLine; inArrays[0] = Indices; mexDisplayArray<uint64_t>(0,NULL,1,inArrays);
//	#endif

	uint64_t *bin = (uint64_t *)mxGetData(binary);
	for (mwSize j = 0; j < dimenstion; j++)
	{
		uint64scalar[0] = j;
		inArrays[0] = Weights[bin[j]]; inArrays[1] = uint64Scalar; outArrays[0] = temp;
		mexExtractColumn<double,uint64_t>(1,outArrays,2,inArrays);
		temp = outArrays[0]; outArrays[0] = NULL;

//		#ifdef SRC_COMPUTEINDICESANDWEIGHTS_HPP_DEBUG
//		std::cout<<"Weights["<<bin[j]<<"]["<<j<<"]:\n";
//		SplitLine; inArrays[0] = temp; mexDisplayArray<float>(0,NULL,1,inArrays);
//		#endif

		inArrays[0] = weight; inArrays[1] = temp; outArrays[0] = weight;
		mexDotProduction<double>(1,outArrays,2,inArrays);
		weight = outArrays[0]; outArrays[0] = NULL;

//		#ifdef SRC_COMPUTEINDICESANDWEIGHTS_HPP_DEBUG
//		std::cout<<"weight:\n"; inArrays[0] = weight; SplitLine; mexDisplayArray<float>(0,NULL,1,inArrays);
//		#endif

	}

//	#ifdef SRC_COMPUTEINDICESANDWEIGHTS_HPP_DEBUG
//	std::cout<<"weight:\n"; inArrays[0] = weight; SplitLine; mexDisplayArray<float>(0,NULL,1,inArrays);
//	#endif

	#ifdef SRC_COMPUTEINDICESANDWEIGHTS_HPP_DETAILS
	SplitLine;
	#endif
//-------------------------------------------------------------------------------------------------
// Free memories:
//-------------------------------------------------------------------------------------------------
	Weights[0] = NULL; Weights[1] = NULL;
	mxDestroyArray(NPoints);
	mxDestroyArray(NDims);
	mxDestroyArray(temp);
	mxDestroyArray(doubleScalar);
	mxDestroyArray(uint64Scalar);
	mxDestroyArray(binary);
	delete inArrays;
	delete outArrays;
	delete Weights;
}
#endif /* SRC_COMPUTEINDICESANDWEIGHTS_HPP_ */

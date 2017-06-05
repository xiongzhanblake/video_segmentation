/*
 * slice_Adjacent.hpp
 *
 *  Created on: Feb 13, 2017
 *      Author: xiong
 */

#ifndef SRC_SLICE_ADJACENT_CPP_
#define SRC_SLICE_ADJACENT_CPP_

// #define SRC_SLICE_ADJACENT_HPP_DEBUG
 #define SRC_SLICE_ADJACENT_HPP__DETAIL

#include "ComputeIndicesAndWeights.hpp"
#include "../../Common/MapOperation.hpp"

void MexSlice_Adjacent( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
//-------------------------------------------------------------------------------------------------
// Check validation:
//-------------------------------------------------------------------------------------------------
	if (nrhs != 6)
		mexErrMsgTxt(Error_INPUT_NUMBER);
	if (nlhs != 1)
		mexErrMsgTxt(Error_OUTPUT_NUMBER);
	if (!mxIsCell(prhs[0]))
		mexErrMsgTxt("Input is not a cell matrix!\n");
	if (mxGetN(prhs[0]) != 1 && mxGetM(prhs[0]) != 1)
		mexErrMsgTxt("Mask has no size as Mx1 or 1xN!\n");
	#ifdef SRC_SLICE_ADJACENT_HPP__DETAIL
	SplitLine; std::cout<<"Function Splice_Adjacent...";
    EndLine; SplitLine;
	#endif
//-------------------------------------------------------------------------------------------------
// Function input and out parameters:
//-------------------------------------------------------------------------------------------------
// Inputs:
    const mxArray **Weights = new const mxArray *[2];
	const mxArray *Data = prhs[0];
	const mxArray *Floors = prhs[1];
	const mxArray *GridSize = prhs[2];
    Weights[0] = prhs[3];
    Weights[1] = prhs[4];
    uint64_t nDims = (uint64_t)((double *)mxGetData(prhs[5]))[0];
// Outputs:
	mxArray *Slice;
//-------------------------------------------------------------------------------------------------
// Local variables:
//-------------------------------------------------------------------------------------------------
	mwSize nClass = mxGetM(Data)>mxGetN(Data)?mxGetM(Data):mxGetN(Data);
	int inNum = 5;
	int outNum = 2;
	const mxArray **inArrays;
	mxArray **outArrays;
	mxArray *NPoints;
	uint64_t *nPoints;
 	mxArray *doubleScalar;
 	double *doublescalar;
	mxArray *uint64Scalar;
	uint64_t *uint64scalar;
	mxArray *Indices;
	mxArray *weight;
	mxArray *Index;
	uint64_t *index;
	mxArray * Grid;
	uint64_t dimensions;
	std::map<double, double>**Accumulator;
	mxArray *temp;
//-------------------------------------------------------------------------------------------------
// Initialization:
//-------------------------------------------------------------------------------------------------
	NPoints = mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
	nPoints = (uint64_t *)mxGetData(NPoints);
	nPoints[0] = mxGetM(Floors);
 	doubleScalar = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
 	doublescalar = (double *)mxGetData(doubleScalar);
	Index = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
	index = (uint64_t *)mxGetData(Index);
	uint64Scalar = mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
	uint64scalar = (uint64_t *)mxGetData(uint64Scalar);
	Grid = NULL;
	temp = NULL;
	Indices = NULL;
	weight = NULL;
	Accumulator = new std::map<double, double>*[nClass];
	for (mwSize i = 0; i < nClass; i ++)
		Accumulator[i] = NULL;

	inArrays = new const mxArray*[inNum];
	outArrays = new mxArray*[outNum];
	for (int i = 0; i < inNum; i ++)
		inArrays[i] = NULL;
	for (int i = 0; i < outNum; i ++)
		outArrays[i] = NULL;

	if (plhs[0] == NULL)
	{
		doublescalar[0] = 0;
		inArrays[0]=NPoints; inArrays[1]=doubleScalar;
		mexAllocateMemoryWithInitials<double>(1,outArrays,2,inArrays);
		plhs[0] = outArrays[0]; outArrays[0] = NULL;
	}

	Slice = plhs[0];

	#ifdef SRC_SLICE_ADJACENT_HPP_DEBUG
	std::cout<<"Slice:\n"; SplitLine; inArrays[0] = Slice; mexDisplayArray<double>(0,NULL,1,inArrays);
	#endif

//-------------------------------------------------------------------------------------------------
// Preparation procedure:
//-------------------------------------------------------------------------------------------------
//    Compute Floors and Weights required;
//-------------------------------------------------------------------------------------------------
	inArrays[0] = GridSize; outArrays[0] = Grid;
	mexFactorialDotProduct<double>(1,outArrays,1,inArrays);
	Grid = outArrays[0]; outArrays[0] = NULL;

	#ifdef SRC_SLICE_ADJACENT_HPP_DEBUG
	std::cout<<"GridSize:\n"; inArrays[0] = Grid; mexDisplayArray<double>(0,outArrays,1,inArrays);
	#endif
//-------------------------------------------------------------------------------------------------
// Computing procedure:
//-------------------------------------------------------------------------------------------------
//    Accumulate indices and weight of vertices along all dimensions:
//-------------------------------------------------------------------------------------------------

//	#ifdef SRC_SLICE_ADJACENT_HPP__DETAIL
//	std::cout<<"Create sparse matrix for Labels:\n";
//	#endif

	for (mwSize n = 0; n < nClass; n++)
		inArrays[n] = mxGetCell(Data,n);
	CreateMapByKeysAndValues<double,double>(nClass, inArrays,Accumulator);

	dimensions = pow(2,nDims);

//	#ifdef SRC_SLICE_ADJACENT_HPP__DETAIL
//	std::cout<<"Compute Weight of video pixels for "<<dimensions<<" dimension:\n"; SplitLine;
//	#endif

	for (uint64_t i = 0; i < dimensions; i ++)
	{
//		#ifdef SRC_SLICE_ADJACENT_HPP__DETAIL
//		SplitLine; std::cout<<i+1<<"th dimension (out of "<<dimensions<<" dimensions):\n";
//		#endif

		*index = i;
		inArrays[0] = Index; inArrays[1] = Floors; inArrays[2] = Grid; inArrays[3] = Weights[0]; inArrays[4] = Weights[1];
		outArrays[0] = Indices; outArrays[1] = weight;
		ComputeIndicesAndWeights( 2, outArrays, 5, inArrays );
		Indices = outArrays[0]; weight = outArrays[1]; outArrays[0] = NULL; outArrays[1] = NULL;
//		std::cout<<"ComputeIndicesAndWeights is over!\n";

		inArrays[0] = Indices; inArrays[1] = weight; outArrays[0] = temp;
		DotProductionMap<uint64_t,double,double,double>(1, outArrays, 2, inArrays, Accumulator);
		temp = outArrays[0]; outArrays[0] = NULL;
//		std::cout<<"DotProductionMap is over!\n";

		inArrays[0] = Slice; inArrays[1]=temp; outArrays[0] = Slice;
		mexAddition<double>(1,outArrays,2,inArrays);
		Slice = outArrays[0]; outArrays[0] = NULL;
//		std::cout<<"mexAddition is over!\n";

		doublescalar[0] = 1;
		inArrays[0]=weight; inArrays[1]=doubleScalar;
		mexResetMemory<double>(0,NULL,2,inArrays);

		uint64scalar[0] = 0;
		inArrays[0]=Indices; inArrays[1]=uint64Scalar;
		mexResetMemory<uint64_t>(0,NULL,2,inArrays);
//		std::cout<<"mexResetMemory is over!\n";
	}
//-------------------------------------------------------------------------------------------------
// Free memories:
//-------------------------------------------------------------------------------------------------
//	SplitLine; std::cout<<"Freeing memory...\n";
	for (uint64_t n = 0; n < nClass; n ++)
		Accumulator[n]->clear();
	delete Accumulator;

	for (int i = 0; i < inNum; i ++)
		inArrays[i] = NULL;
	for (int i = 0; i < outNum; i ++)
		outArrays[i] = NULL;

	mxDestroyArray(NPoints);
	mxDestroyArray(weight);
	mxDestroyArray(Indices);
	mxDestroyArray(Grid);
	mxDestroyArray(Index);
 	mxDestroyArray(doubleScalar);
	mxDestroyArray(uint64Scalar);
    	mxDestroyArray(temp);

	delete inArrays;
	delete outArrays;
	delete Weights;
//	SplitLine; std::cout<<"Slice is over!\n";
}

#endif /* SRC_SLICE_ADJACENT_CPP_BACKUP_ */

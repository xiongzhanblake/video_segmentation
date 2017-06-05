/*
 * Splat_Adjacent.hpp
 *
 *  Created on: Feb 7, 2017
 *  Author: xiong
 */

#ifndef SPLAT_ADJACENT_HPP_
#define SPLAT_ADJACENT_HPP_

#include "ComputeIndicesAndWeights.hpp"
#include "../../Common/MapOperation.hpp"

// #define SPLAT_ADJACENT_HPP_DEBUG
#define SPLAT_ADJACENT_HPP_DETAILS
//-------------------------------------------------------------------------------------------------
//Main function:
//-------------------------------------------------------------------------------------------------
void MexSplat_Adjacent( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
//-------------------------------------------------------------------------------------------------
// Check validation:
//-------------------------------------------------------------------------------------------------
	if (nrhs != 6)
		mexErrMsgTxt(Error_INPUT_NUMBER);
	if (nlhs != 1)
		mexErrMsgTxt(Error_OUTPUT_NUMBER);

	#ifdef SPLAT_ADJACENT_HPP_DETAILS
	SplitLine; std::cout<<"Function Splat_Adjacent...";
    EndLine; SplitLine;
	#endif
    // Floor, weights, BilateralVals GridSize are all double type;
    const mxArray **Weights = new const mxArray *[2];
    
	const mxArray *Floor = prhs[0];
	const mxArray *BilateralVals = prhs[1];
	const mxArray *GridSize = prhs[2];
    Weights[0] = prhs[3];
    Weights[1] = prhs[4];
    uint64_t nDims = (uint64_t)((double *)mxGetData(prhs[5]))[0];
    
	mxArray *NPoints = mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
	uint64_t *nPoints = (uint64_t *)mxGetData(NPoints);
	*nPoints = mxGetM(Floor);
	mxArray *NClass = mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
	uint64_t *nClass = (uint64_t *)mxGetData(NClass);
	nClass[0] = mxGetN(BilateralVals);
    mxArray *Grid = NULL;
    mxArray *doubleScalar = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
 	double *doublescalar = (double *)mxGetData(doubleScalar);
//-------------------------------------------------------------------------------------------------
// Preparation procedure:
//-------------------------------------------------------------------------------------------------
// 1. Horizental catanetation Visual cues, Motion cues and their cooresponding grid size;
//-------------------------------------------------------------------------------------------------
	int inNum = 5;
	int outNum = 2;
	const mxArray **inArrays = new const mxArray*[inNum];
	for (int i = 0; i < inNum; i ++)
		inArrays[i] = NULL;
	mxArray **outArrays = new mxArray*[outNum];
	for (int i = 0; i < outNum; i ++)
			outArrays[i] = NULL;

//    #ifdef SPLAT_ADJACENT_HPP_DETAILS
//	SplitLine; std::cout<<"Compute indices and Weights:";
//    EndLine;
//	#endif
//-------------------------------------------------------------------------------------------------
// 2. Compute Indices and Weights required;
//-------------------------------------------------------------------------------------------------
	inArrays[0] = GridSize; outArrays[0] = Grid;
	mexFactorialDotProduct<double>(1,outArrays,1,inArrays);
	Grid = outArrays[0]; outArrays[0] = NULL;

	#ifdef SPLAT_ADJACENT_HPP_DEBUG
	std::cout<<"GridSize:\n"; inArrays[0] = Grid; mexDisplayArray<double>(0,NULL,1,inArrays);
	#endif

	mxArray *Index = mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
	uint64_t *index = (uint64_t *)mxGetData(Index);
	mxArray *uint64Scalar = mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
	uint64_t *uint64scalar = (uint64_t *)mxGetData(uint64Scalar);

	std::map<uint64_t, double> **accumulator = new std::map<uint64_t, double> *[nClass[0]];
	for (uint64_t n = 0; n < nClass[0]; n ++)
		accumulator[n] = NULL;
	mxArray *weight = NULL;
	mxArray *Indices = NULL;
	uint64_t dimensions = pow(2,nDims);

//	#ifdef SPLAT_ADJACENT_HPP_DETAILS
//	std::cout<<"Compute Weight of video pixels for "<<dimensions<<" dimension:\n"; SplitLine;
//	#endif

//-------------------------------------------------------------------------------------------------
// Computing procedure:
//-------------------------------------------------------------------------------------------------
// 1. Accumulate indices and weight of vertices along all dimensions:
//-------------------------------------------------------------------------------------------------
	for (uint64_t i = 0; i < dimensions; i ++)
	{
//		#ifdef SPLAT_ADJACENT_HPP_DETAILS
//		std::cout<<i+1<<"th dimension (out of "<<dimensions<<" dimensions):\n";
//		#endif

		*index = i;
		inArrays[0] = Index; inArrays[1] = Floor; inArrays[2] = Grid; inArrays[3] = Weights[0]; inArrays[4] = Weights[1];
		outArrays[0] = Indices; outArrays[1] = weight;
		ComputeIndicesAndWeights( 2, outArrays, 5, inArrays );
		Indices = outArrays[0]; weight = outArrays[1]; outArrays[0] = NULL; outArrays[1] = NULL;

		inArrays[0] = Indices; inArrays[1] = weight; inArrays[2] = BilateralVals;
		AccMapSum<uint64_t,double,uint64_t,double>( 3, inArrays, accumulator );

		doublescalar[0] = 1;
		inArrays[0]=weight; inArrays[1]=doubleScalar;
		mexResetMemory<double>(0,NULL,2,inArrays);

		uint64scalar[0] = 0;
		inArrays[0]=Indices; inArrays[1]=uint64Scalar;
		mexResetMemory<uint64_t>(0,NULL,2,inArrays);

	}
//-------------------------------------------------------------------------------------------------
// 2. Create a sparse matrix for all occupied vertices:
//-------------------------------------------------------------------------------------------------
	SplitMap<uint64_t,double>(nlhs,plhs,nClass[0],accumulator);

	#ifdef SPLAT_ADJACENT_HPP_DEBUG
	SplitLine; std::cout<<"Display contents of outputs:\n";
	SplitLine; std::cout<<"Format:\n"<<"Class index:\n\t<Iteration, index, weight>\n";
	uint64_t length = mxGetM(plhs[0]);
	uint64_t *ptr1 = (uint64_t *)mxGetData(plhs[0]);
	double *ptr2 = (double *)mxGetData(plhs[1]);
	for (uint64_t n = 0; n < nClass[0]; n ++)
	{
		SplitLine; std::cout<<"Class "<<n<<":\n";
		for( uint64_t index = 0; index < length; index++ )
			std::cout<<"\t"<<index<<"\t"<<"<"<<ptr1[index+n*length]<<","<<ptr2[index+n*length]<<">\t"<<"\n";
		SplitLine;
	}
	#endif
//-------------------------------------------------------------------------------------------------
// Post-process procedure:
//-------------------------------------------------------------------------------------------------
//     Free dynamic memories and delete structures and classes safely;
//-------------------------------------------------------------------------------------------------
	for (uint64_t n = 0; n < nClass[0]; n ++)
	{
		accumulator[n]->clear();
		delete accumulator[n];
	}
	delete accumulator;
	for (int i = 0; i < inNum; i ++)
		inArrays[i] = NULL;
	for (int i = 0; i < outNum; i ++)
		outArrays[i] = NULL;

	mxDestroyArray(weight);
	mxDestroyArray(Grid);
	mxDestroyArray(Indices);
	mxDestroyArray(Index);
	mxDestroyArray(uint64Scalar);
	mxDestroyArray(NPoints);
	mxDestroyArray(doubleScalar);
	mxDestroyArray(NClass);
	delete outArrays;
	delete inArrays;
	delete Weights;
//-------------------------------------------------------------------------------------------------
}
//-------------------------------------------------------------------------------------------------

#endif /* SPLAT_ADJACENT_HPP_ */

/*
 * BasicMexMatrixOperation.hpp
 *
 *  Created on: Feb 13, 2017
 *      Author: xiong
 */
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//	Some important tips before use:
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//  I. Arrangment about matrices:
//-------------------------------------------------------------------------------------------------
//  1. Base index starts from 0 (accordiance with C/C++);
//  2. Row-first Principle (accordiance with MATLAB);
//-------------------------------------------------------------------------------------------------
// II. Type coorespondance (from MATLAB ClassID document):
// Index 	MATLAB		mex-type	 C/C++
//-------------------------------------------------------------------------------------------------
//   1.		Single		float		 float
//   2.		double		double		 double
//   3.		int8		int8_t		 char/byte
//   4.		int16		int16_t		 short
//   5.		int32		int32_t		 int
//   6.		int64		int64_t		 long long
//   7.		uint8		uint8_t		 unsigned char/byte
//   8.		uint16		uint16_t	 unsigned short
//   9.		uint32		uint32_t	 unsigned int
//  10.		uint64		uint64_t	 unsigned long long
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//	Debug mode with debuging information:
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
#ifndef SRC_BASICMEXMATRIXOPERATIONS_HPP_
#define SRC_BASICMEXMATRIXOPERATIONS_HPP_

#define SRC_BASICMEXMATRIXOPERATIONS_HPP_DEBUG
#define SRC_BASICMEXMATRIXOPERATIONS_HPP_DETAIL
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//	User-defined Macros:
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
#define Error_INPUT_NUMBER "Error: The number of input(s) is not correct\n"
#define Error_OUTPUT_NUMBER "Error: The number of output(s) is not correct\n"
#define Error_TYPE "Error: The data type is not matached\n"
#define Error_DIMS "Error: The dimension(s) dismatached\n"
#define SplitLine std::cout<<"----------------------------------------------------------------\n"
#define EndLine std::cout<<"\n"
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//	Head files:
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
#include <matrix.h>
#include <mex.h>
#include "CUtilities.hpp"
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//	MEX Functions:
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//	1. Basic operation:
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//	Function mexReshape: only change the dimension-array with its elements unchanged. (Untested yet)
//-------------------------------------------------------------------------------------------------
//	Inputs:		4
//-------------------------------------------------------------------------------------------------
//	nlsh: output number;
//	plhs: output pointer-array;
//	nrhs: input number;
//	prhs: input	pointer-array;
//-------------------------------------------------------------------------------------------------
//	OutPuts: 	0
//-------------------------------------------------------------------------------------------------
// Notice:
// 1. The number of total element is unchanged all the time!
// 2. The element-sequence is never changed, and give the current matrix a newer explannation!
//-------------------------------------------------------------------------------------------------
void mexReshape(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const int iNum = 2, oNum = 1;
	const char type[] = "uint8";
	const char fun_name[] = "mexReshape";

	if (nrhs != iNum)
	{
		mexPrintf("In function %s:\t",fun_name);
		mexPrintf(Error_INPUT_NUMBER);
		mexPrintf("Expects %d instead of %d!\n", iNum ,nrhs);
		mexErrMsgTxt("");
	}
	if (nlhs != oNum)
	{
		mexPrintf("In function %s:\t",fun_name);
		mexPrintf(Error_OUTPUT_NUMBER);
		mexPrintf("Expects %d instead of %d!\n", oNum, nlhs);
		mexErrMsgTxt("");
	}
	if (mxGetClassID(prhs[1]) != mxUINT8_CLASS)
	{
		mexPrintf("In function %s:\t",fun_name);
		mexPrintf(Error_TYPE);
		mexPrintf("Expects %s instead of %s!\n", type, mxGetClassName(prhs[1]));
		mexErrMsgTxt("");
	}
	mwSize *shape = (mwSize *)mxGetData(prhs[1]);
	if (mxGetNumberOfDimensions(prhs[1]) != 2)
	{
		mexPrintf("In function %s:\t",fun_name);
		mexPrintf("Error: The shape is not a 2-D matrix!\n");
		mexErrMsgTxt("");
	}
	if (mxGetM(prhs[1])!=1 && mxGetN(prhs[1])!=1)
	{
		mexPrintf("In function %s:\t",fun_name);
		mexPrintf("Error: The shape is not a vetor!\n");
		mexErrMsgTxt("");
	}
	mwSize ndim = (mxGetM(prhs[1])>mxGetN(prhs[1]))?mxGetM(prhs[1]):mxGetN(prhs[1]);
	uint64_t length = CFactorialProduct<mwSize,mwSize,uint64_t>(shape,ndim);
	if (length != mxGetNumberOfElements(prhs[0]))
	{
		mexPrintf("In function %s:\t",fun_name);
		mexPrintf("Error: The total element number of new shape differs with that of the original!\n");
		mexErrMsgTxt("");
	}
	if ((plhs[0] != NULL && length != mxGetNumberOfElements(plhs[0])) || plhs[0] == NULL)
	{
		mexPrintf("In function %s:\t",fun_name);
		mexPrintf("The total element number of the reshaped differs from that of the original!\n");
		mexPrintf("Destroy the old reshaped matrix and create a new one!\n");
		mxDestroyArray(plhs[0]);
		plhs[0] = NULL;
		plhs[0] = mxDuplicateArray(prhs[0]);
	}
	mxSetDimensions(plhs[0], shape, ndim);
}
//-------------------------------------------------------------------------------------------------
//	Function mexPermute: Rearrange dimensions of a N-D matrix. (Untested yet)
//-------------------------------------------------------------------------------------------------
//	Inputs:		4
//-------------------------------------------------------------------------------------------------
//	nlsh: output number;
//	plhs: output pointer-array;
//	nrhs: input number;
//	prhs: input	pointer-array;
//-------------------------------------------------------------------------------------------------
//	OutPuts: 	0
//-------------------------------------------------------------------------------------------------
// Notice:
// 1. a permutation operation is an extension of transpose operation!
// 2. order vetocr start index is from 0!
//-------------------------------------------------------------------------------------------------
template <typename T>
void mexPermute(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const int iNum = 2, oNum = 1;
	const char type[] = "uint8";
	const char fun_name[] = "mexPermute";

	if (nrhs != iNum)
	{
		mexPrintf("In function %s:\t",fun_name);
		mexPrintf(Error_INPUT_NUMBER);
		mexPrintf("Expects %d instead of %d!\n", iNum ,nrhs);
		mexErrMsgTxt("");
	}
	if (nlhs != oNum)
	{
		mexPrintf("In function %s:\t",fun_name);
		mexPrintf(Error_OUTPUT_NUMBER);
		mexPrintf("Expects %d instead of %d!\n", oNum, nlhs);
		mexErrMsgTxt("");
	}
	if (mxGetClassID(prhs[1]) != mxUINT8_CLASS)
	{
		mexPrintf("In function %s:\t",fun_name);
		mexPrintf(Error_TYPE);
		mexPrintf("Expects %s instead of %s!\n", type, mxGetClassName(prhs[1]));
		mexErrMsgTxt("");
	}
	mwSize *order = (mwSize *)mxGetData(prhs[1]);
	uint64_t length = mxGetNumberOfElements(prhs[0]);
	if (mxGetNumberOfDimensions(prhs[1]) != 2)
	{
		mexPrintf("In function %s:\t",fun_name);
		mexPrintf("Error: The order is not 2-D matrix!\n");
		mexErrMsgTxt("");
	}
	if (mxGetM(prhs[1])!=1 && mxGetN(prhs[1])!=1)
	{
		mexPrintf("In function %s:\t",fun_name);
		mexPrintf("Error: The order is not a vetor!\n");
		mexErrMsgTxt("");
	}
	mwSize ndim = (mxGetM(prhs[1])>mxGetN(prhs[1]))?mxGetM(prhs[1]):mxGetN(prhs[1]);
	if (ndim != mxGetNumberOfDimensions(prhs[0]))
	{
		mexPrintf("In function %s:\t",fun_name);
		mexPrintf("Error: The dimension is not the same with the original!\n");
		mexErrMsgTxt("");
	}
	mwSize min = CMin<mwSize>(order,ndim);
	if (min > 0)
	{
//		mexPrintf("In function %s:\t",fun_name);
//		mexPrintf("The start index of order vector is from 0\n");
		CMinusAS<mwSize>(order,ndim, min, order);
	}
	mwSize *old_dims = (mwSize *)mxGetDimensions(prhs[0]);
	mwSize *new_dims = new mwSize [ndim];
	mwSize *old_weight = new mwSize [ndim];
	mwSize *new_weight = new mwSize [ndim];
	mwSize *pos = new mwSize [ndim];
	for (mwSize i = 0; i < ndim; i ++)
	{
		new_dims[i] = old_dims[order[i]];
		pos[order[i]] = i;
	}
	CFactorialDotProduct<mwSize,mwSize>(old_dims,ndim,old_weight);
	CFactorialDotProduct<mwSize,mwSize>(new_dims,ndim,new_weight);
	mxArray *temporal = NULL;
	if ((plhs[0] != NULL && length != mxGetNumberOfElements(plhs[0])) || plhs[0] == NULL)
	{
//		mexPrintf("In function %s:\t",fun_name);
//		mexPrintf("The total element number of the permuted differs from that of the original!\n");
//		mexPrintf("Destroy the old reshaped matrix and create a new one!\n");
		mxDestroyArray(plhs[0]);
		plhs[0] = NULL;
		plhs[0] = mxCreateNumericArray(ndim,new_dims,mxGetClassID(prhs[0]),mxREAL);
	}
	else if (prhs[0] == plhs[0])
		temporal = mxDuplicateArray(prhs[0]);
	//Rearrange elements and dimensions
	T *array = NULL;
	if (temporal == NULL)
		array = (T *)mxGetData(prhs[0]);
	else
		array = (T *)mxGetData(temporal);
	T *result = (T *)mxGetData(plhs[0]);
	uint64_t accumulator, temp, offset;
	for (uint64_t index = 0; index < length; index ++)
	{
		accumulator = 0; temp = index;
		for	(int n = ndim - 1; n >= 0; n --)
		{
			offset = temp /old_weight[n];
			temp = temp - offset * old_weight[n];
			accumulator += offset * new_weight[pos[n]];
		}
		result[accumulator] = array[index];
	}
	if (temporal != NULL)
	{
		mxDestroyArray(temporal);
		temporal = NULL;
	}
	mxSetDimensions(plhs[0], new_dims, ndim);
	delete new_dims;
	delete old_weight;
	delete new_weight;
	delete pos;
}
//-------------------------------------------------------------------------------------------------
//	Function mexPermute: Rearrange dimensions of a N-D matrix iversely. (Untested yet)
//-------------------------------------------------------------------------------------------------
//	Inputs:		4
//-------------------------------------------------------------------------------------------------
//	nlsh: output number;
//	plhs: output pointer-array;
//	nrhs: input number;
//	prhs: input	pointer-array;
//-------------------------------------------------------------------------------------------------
//	OutPuts: 	0
//-------------------------------------------------------------------------------------------------
// Notice:
// 1. a permutation operation is an extension of transpose operation!
// 2. order vetocr start index is from 0!
//-------------------------------------------------------------------------------------------------
template <typename T>
void mexiPermute(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const int iNum = 2, oNum = 1;
	const char type[] = "uint8";
	const char fun_name[] = "mexPermute";

	if (nrhs != iNum)
	{
		mexPrintf("In function %s:\t",fun_name);
		mexPrintf(Error_INPUT_NUMBER);
		mexPrintf("Expects %d instead of %d!\n", iNum ,nrhs);
		mexErrMsgTxt("");
	}
	if (nlhs != oNum)
	{
		mexPrintf("In function %s:\t",fun_name);
		mexPrintf(Error_OUTPUT_NUMBER);
		mexPrintf("Expects %d instead of %d!\n", oNum, nlhs);
		mexErrMsgTxt("");
	}
	if (mxGetClassID(prhs[1]) != mxUINT8_CLASS)
	{
		mexPrintf("In function %s:\t",fun_name);
		mexPrintf(Error_TYPE);
		mexPrintf("Expects %s instead of %s!\n", type, mxGetClassName(prhs[1]));
		mexErrMsgTxt("");
	}
	mwSize *order = (mwSize *)mxGetData(prhs[1]);
	uint64_t length = mxGetNumberOfElements(prhs[0]);
	if (mxGetNumberOfDimensions(prhs[1]) != 2)
	{
		mexPrintf("In function %s:\t",fun_name);
		mexPrintf("Error: The order is not 2-D matrix!\n");
		mexErrMsgTxt("");
	}
	if (mxGetM(prhs[1])!=1 && mxGetN(prhs[1])!=1)
	{
		mexPrintf("In function %s:\t",fun_name);
		mexPrintf("Error: The order is not a vetor!\n");
		mexErrMsgTxt("");
	}
	mwSize ndim = (mxGetM(prhs[1])>mxGetN(prhs[1]))?mxGetM(prhs[1]):mxGetN(prhs[1]);
	if (ndim != mxGetNumberOfDimensions(prhs[0]))
	{
		mexPrintf("In function %s:\t",fun_name);
		mexPrintf("Error: The dimension is not the same with the original!\n");
		mexErrMsgTxt("");
	}
	mwSize min = CMin<mwSize>(order,ndim);
	if (min > 0)
	{
//		mexPrintf("In function %s:\t",fun_name);
//		mexPrintf("The start index of order vector is from 0\n");
		CMinusAS<mwSize>(order,ndim, min, order);
	}
	mwSize *old_dims = (mwSize *)mxGetDimensions(prhs[0]);
	mwSize *new_dims = new mwSize [ndim];
	mwSize *old_weight = new mwSize [ndim];
	mwSize *new_weight = new mwSize [ndim];
	mwSize *pos = new mwSize [ndim];
	for (mwSize i = 0; i < ndim; i ++)
	{
		new_dims[order[i]] = old_dims[i];
		pos[i] = order[i];
	}
	CFactorialDotProduct<mwSize,mwSize>(old_dims,ndim,old_weight);
	CFactorialDotProduct<mwSize,mwSize>(new_dims,ndim,new_weight);
	mxArray *temporal = NULL;
	if ((plhs[0] != NULL && length != mxGetNumberOfElements(plhs[0])) || plhs[0] == NULL)
	{
//		mexPrintf("In function %s:\t",fun_name);
//		mexPrintf("The total element number of the permuted differs from that of the original!\n");
//		mexPrintf("Destroy the old reshaped matrix and create a new one!\n");
		mxDestroyArray(plhs[0]);
		plhs[0] = NULL;
		plhs[0] = mxCreateNumericArray(ndim,new_dims,mxGetClassID(prhs[0]),mxREAL);
	}
	else if (prhs[0] == plhs[0])
		temporal = mxDuplicateArray(prhs[0]);
	//Rearrange elements and dimensions
	T *array = NULL;
	if (temporal == NULL)
		array = (T *)mxGetData(prhs[0]);
	else
		array = (T *)mxGetData(temporal);
	T *result = (T *)mxGetData(plhs[0]);
	uint64_t accumulator, temp, offset;
	for (uint64_t ind = 0; ind < length; ind ++)
	{
		accumulator = 0; temp = ind;
		for	(int n = ndim - 1; n >= 0; n --)
		{
			offset = temp /old_weight[n];
			temp = temp - offset * old_weight[n];
			accumulator += offset * new_weight[pos[n]];
		}
		result[accumulator] = array[ind];
	}
	if (temporal != NULL)
	{
		mxDestroyArray(temporal);
		temporal = NULL;
	}
	mxSetDimensions(plhs[0], new_dims, ndim);
	delete new_dims;
	delete old_weight;
	delete new_weight;
	delete pos;
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//	2. Catenation operation:
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//	Function mexhorzcat: horzentally catanenate two matrices. (Untested yet)
//-------------------------------------------------------------------------------------------------
//	Inputs:		4
//-------------------------------------------------------------------------------------------------
//	nlsh: output number;
//	plhs: output pointer-array;
//	nrhs: input number;
//	prhs: input	pointer-array;
//-------------------------------------------------------------------------------------------------
//	OutPuts: 	0
//-------------------------------------------------------------------------------------------------
template<typename T>
void mexHorzcat(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 2)
	{
		mexPrintf("In function mexhorzcat:\t");
		mexPrintf(Error_INPUT_NUMBER);
		mexPrintf("The number of input(s) expects 2 instead of %d!\n", nrhs);
		mexErrMsgTxt("");
	}
	if (nlhs != 1)
	{
		mexPrintf("In function mexhorzcat:\t");
		mexPrintf(Error_OUTPUT_NUMBER);
		mexPrintf("The number of output(s) expects 1 instead of %d!\n", nlhs);
		mexErrMsgTxt("");
	}

	if (prhs[0] == NULL && prhs[1] == NULL) mexErrMsgTxt("No input is valid!\n");

	if (prhs[0] == NULL || prhs[1] == NULL)
	{
		if (prhs[0] == NULL) plhs[0] = mxDuplicateArray(prhs[1]);
		else plhs[0] = mxDuplicateArray(prhs[0]);
	}	
	else
	{
		uint64_t sz1[2],sz2[2];
		if(mxGetNumberOfDimensions(prhs[0]) !=2 || mxGetNumberOfDimensions(prhs[1]) != 2)
			mexErrMsgTxt("This function only concatenates 2D matrices");
		sz1[0] = mxGetM(prhs[0]); sz1[1] = mxGetN(prhs[0]);
		sz2[0] = mxGetM(prhs[1]); sz2[1] = mxGetN(prhs[1]);
		if (sz1[0] != sz2[0])
			mexErrMsgTxt("Dimension mismatch");
		uint64_t nbytes, step;
		step = sz1[1] + sz2[1];
		nbytes = sz1[0];//size of new array
		if (plhs[0] != NULL) mxDestroyArray(plhs[0]);
		plhs[0] = mxCreateNumericMatrix(nbytes,step,mxGetClassID(prhs[0]),mxREAL);
		T *result = (T *)mxGetData(plhs[0]);
		T *array1 = (T *)mxGetData(prhs[0]);
		T *array2 = (T *)mxGetData(prhs[1]);
		for (uint64_t row = 0; row < sz1[0]; row++)
		{
			for (uint64_t col = 0; col < sz1[1]; col++)
				result[row + col * nbytes] = array1[row + nbytes * col];
			for (uint64_t col = 0; col < sz2[1]; col++)
				result[row + (col+sz1[1]) * nbytes]=array2[row + nbytes * col];
		}
	}
}
//-------------------------------------------------------------------------------------------------
//	Function mexVertcat: vertically catanenate two matrices. (Untested yet)
//-------------------------------------------------------------------------------------------------
//	Inputs:		4
//-------------------------------------------------------------------------------------------------
//	nlsh: output number;
//	plhs: output pointer-array;
//	nrhs: input number;
//	prhs: input	pointer-array;
//-------------------------------------------------------------------------------------------------
//	OutPuts: 	0
//-------------------------------------------------------------------------------------------------
template<typename T>
void mexVertcat(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 2)
		mexErrMsgTxt("2 inputs are required!\n");
	if (nlhs != 1)
		mexErrMsgTxt("1 output is required!\n");
	if (prhs[0] == NULL && prhs[1] == NULL)
		mexErrMsgTxt("No valid input!\n");
	if (prhs[0] != NULL && prhs[1] != NULL)
	{
		uint64_t sz1[2],sz2[2];
		if(mxGetNumberOfDimensions(prhs[0]) !=2 || mxGetNumberOfDimensions(prhs[1]) != 2)
			mexErrMsgTxt("This function only concatenates 2D matrices!\n");
		sz1[0] = mxGetM(prhs[0]); sz1[1] = mxGetN(prhs[0]);
		sz2[0] = mxGetM(prhs[1]); sz2[1] = mxGetN(prhs[1]);
		if (sz1[1] != sz2[1])
			mexErrMsgTxt("Dimension mismatch");
		uint64_t nbytes, step;
		step = sz1[1];
		nbytes = sz1[0] + sz2[0];//size of new array

		mxArray *temp = mxCreateNumericMatrix(nbytes,step,mxGetClassID(prhs[0]),mxREAL);

		T *result = (T *)mxGetData(temp);
		T *array1 = (T *)mxGetData(prhs[0]);
		T *array2 = (T *)mxGetData(prhs[1]);
		for (uint64_t col = 0; col < sz1[1]; col++)
		{
			for (uint64_t row = 0; row < sz1[0]; row++)
				result[row + col * nbytes] = array1[row + col * sz1[0]];
			for (uint64_t row = 0; row < sz2[0]; row++)
				result[row + sz1[0] + col * nbytes] = array2[row + col * sz2[0]];
		}
		if (plhs[0] != NULL)
			mxDestroyArray(plhs[0]);
		plhs[0] = temp;
	}
	else
	{
		if (prhs[0] == NULL) plhs[0] = mxDuplicateArray(prhs[1]);
		else plhs[0] = mxDuplicateArray(prhs[0]);
	}
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//	3. Element-wise Operation:
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//	Function mexCeil: Call ceil function element-wisely. (Untested yet)
//-------------------------------------------------------------------------------------------------
//	Inputs:		4
//-------------------------------------------------------------------------------------------------
//	nlsh: output number;
//	plhs: output pointer-array;
//	nrhs: input number;
//	prhs: input	pointer-array;
//-------------------------------------------------------------------------------------------------
//	OutPuts: 	0
//-------------------------------------------------------------------------------------------------
template<typename T>
void mexCeil(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 1)
		mexErrMsgTxt("1 inputs is required!\n");
	if (nlhs != 1)
		mexErrMsgTxt("1 input is required!\n");
	uint64_t sz[2];
	if(mxGetNumberOfDimensions(prhs[0]) !=2)
		mexErrMsgTxt("This function only supports 2D matrices");
	sz[0] = mxGetM(prhs[0]); sz[1] = mxGetN(prhs[0]);
	if (plhs[0] == NULL)
		plhs[0] = mxCreateNumericMatrix(sz[0],sz[1],mxGetClassID(prhs[0]),mxREAL);
	T *result = (T *)mxGetData(plhs[0]);
	T *array = (T *)mxGetData(prhs[0]);
	for (uint64_t i = 0; i < sz[0]*sz[1]; i++)
		result[i] = (T)ceil(array[i]);
}
//-------------------------------------------------------------------------------------------------
//	Function mexRound: Call round function element-wisely. (Untested yet)
//-------------------------------------------------------------------------------------------------
//	Inputs:		4
//-------------------------------------------------------------------------------------------------
//	nlsh: output number;
//	plhs: output pointer-array;
//	nrhs: input number;
//	prhs: input	pointer-array;
//-------------------------------------------------------------------------------------------------
//	OutPuts: 	0
//-------------------------------------------------------------------------------------------------
template<typename T>
void mexRound(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 1)
		mexErrMsgTxt("1 inputs is required!\n");
	if (nlhs != 1)
		mexErrMsgTxt("1 input is required!\n");
	uint64_t sz[2];
	if(mxGetNumberOfDimensions(prhs[0]) !=2)
		mexErrMsgTxt("This function only supports 2D matrices");
	sz[0] = mxGetM(prhs[0]); sz[1] = mxGetN(prhs[0]);
	if (plhs[0] == NULL)
		plhs[0] = mxCreateNumericMatrix(sz[0],sz[1],mxGetClassID(prhs[0]),mxREAL);
	T *result = (T *)mxGetData(plhs[0]);
	T *array = (T *)mxGetData(prhs[0]);
	for (uint64_t i = 0; i < sz[0]*sz[1]; i++)
		result[i] = (T)round(array[i]);
}
//-------------------------------------------------------------------------------------------------
//	Function mexFloor: Call floor function element-wisely. (Untested yet)
//-------------------------------------------------------------------------------------------------
//	Inputs:		4
//-------------------------------------------------------------------------------------------------
//	nlsh: output number;
//	plhs: output pointer-array;
//	nrhs: input number;
//	prhs: input	pointer-array;
//-------------------------------------------------------------------------------------------------
//	OutPuts: 	0
//-------------------------------------------------------------------------------------------------
template<typename T>
void mexFloor(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//    #ifdef SRC_BASICMEXMATRIXOPERATIONS_HPP_DEBUG
//    SplitLine;
//    std::cout<<"Function mexFloor parameter Check:";
//    EndLine;
//    #endif

//    #ifdef SRC_BASICMEXMATRIXOPERATIONS_HPP_DEBUG
//	SplitLine;
//    std::cout<<"Check: nrhs and nlhs\t=>\t";
//	#endif
    
	if (nrhs != 1)
		mexErrMsgTxt("1 inputs is required!\n");
	if (nlhs != 1)
		mexErrMsgTxt("1 output is required!\n");
    
//    #ifdef SRC_BASICMEXMATRIXOPERATIONS_HPP_DEBUG
//	SplitLine;
//    std::cout<<"Passed!";
//    EndLine;
//	#endif
//
//    #ifdef SRC_BASICMEXMATRIXOPERATIONS_HPP_DEBUG
//	SplitLine;
//    std::cout<<"Check: NumberOfDimensions\t=>\t";
//	#endif
    
	if(mxGetNumberOfDimensions(prhs[0]) !=2)
		mexErrMsgTxt("This function only supports 2D matrices");
    
//    #ifdef SRC_BASICMEXMATRIXOPERATIONS_HPP_DEBUG
//	SplitLine;
//    std::cout<<"Passed!";
//    EndLine;
//	#endif

    uint64_t sz[2];
	sz[0] = mxGetM(prhs[0]); sz[1] = mxGetN(prhs[0]);
	if (plhs[0] == NULL)
		plhs[0] = mxCreateNumericMatrix(sz[0],sz[1],mxGetClassID(prhs[0]),mxREAL);
	T *result = (T *)mxGetData(plhs[0]);
	T *array = (T *)mxGetData(prhs[0]);
	for (uint64_t i = 0; i < sz[0]*sz[1]; i++)
		result[i] = (T)floor(array[i]);
}
//-------------------------------------------------------------------------------------------------
//	Function mexAddition: add two matrices element-wisely. (Untested yet)
//-------------------------------------------------------------------------------------------------
//	Inputs:		4
//-------------------------------------------------------------------------------------------------
//	nlsh: output number;
//	plhs: output pointer-array;
//	nrhs: input number;
//	prhs: input	pointer-array;
//-------------------------------------------------------------------------------------------------
//	OutPuts: 	0
//-------------------------------------------------------------------------------------------------
// Notice: supports 4 verisons: scalar-matrix, matrix-scalar, scalar-scalar, matrix-matrix.
//-------------------------------------------------------------------------------------------------
template<typename T>
void mexAddition(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 2)
		mexErrMsgTxt("2 inputs are required!\n");
	if (nlhs != 1)
		mexErrMsgTxt("1 output is required!\n");
	uint64_t sz1[2],sz2[2];
	if (mxGetNumberOfDimensions(prhs[0]) !=2 )
		mexErrMsgTxt("Only 2-D matrix Addition is supported!");
	sz1[0] = mxGetM(prhs[0]); sz1[1] = mxGetN(prhs[0]);
	sz2[0] = mxGetM(prhs[1]); sz2[1] = mxGetN(prhs[1]);
	if (sz1[0] == sz2[0] && sz1[1] == sz2[1])
	{

		#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
		mexPrintf("2-D Matrices Addition!\n");
		#endif

		if (plhs[0] == NULL)
			plhs[0] = mxCreateNumericMatrix(sz1[0],sz1[1], mxGetClassID(prhs[0]),mxREAL);
		T *result = (T *)mxGetData(plhs[0]);
		T *array1 = (T *)mxGetData(prhs[0]);
		T *array2 = (T *)mxGetData(prhs[1]);
		for (uint64_t i = 0; i < sz1[0]*sz1[1]; i++)
			result[i] = array1[i] + array2[i];
	}
	else
	{

		#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
		mexPrintf("Addition:\t");
		#endif

		if (sz1[0] == sz2[0] && sz2[0] > 1 && sz2[1] == 1)
		{

			#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
			mexPrintf("Matrix-ColVector Addition!\n");
			#endif

			if (plhs[0] == NULL)
				plhs[0] = mxCreateNumericMatrix(sz1[0],sz1[1], mxGetClassID(prhs[0]),mxREAL);
			T *result = (T *)mxGetData(plhs[0]);
			T *array1 = (T *)mxGetData(prhs[0]);
			T *array2 = (T *)mxGetData(prhs[1]);
			uint64_t step = sz1[0];
			for (uint64_t col = 0; col < sz1[1]; col++)
				for (uint64_t row = 0; row < sz1[0]; row++)
					result[row + col * step] = array1[row + col * step] + array2[row];
		}
		else if (sz1[1] == sz2[1] && sz2[0] == 1 && sz2[1] > 1)
		{

			#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
			mexPrintf("Matrix-RowVector Addition!\n");
			#endif

			if (plhs[0] == NULL)
				plhs[0] = mxCreateNumericMatrix(sz1[0],sz1[1], mxGetClassID(prhs[0]),mxREAL);
			T *result = (T *)mxGetData(plhs[0]);
			T *array1 = (T *)mxGetData(prhs[0]);
			T *array2 = (T *)mxGetData(prhs[1]);
			uint64_t step = sz1[0];
			for (uint64_t col = 0; col < sz1[1]; col++)
				for (uint64_t row = 0; row < sz1[0]; row++)
					result[col * step + row] = array1[col * step + row] + array2[col];
		}
		else if ( sz2[0] == 1 && sz2[1] == 1 )
		{

			#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
			mexPrintf("Matrix-Scalar Addition!\n");
			#endif

			if (plhs[0] == NULL)
				plhs[0] = mxCreateNumericMatrix(sz1[0],sz1[1], mxGetClassID(prhs[0]),mxREAL);
			T *result = (T *)mxGetData(plhs[0]);
			T *array1 = (T *)mxGetData(prhs[0]);
			T *array2 = (T *)mxGetData(prhs[1]);
			for (uint64_t i = 0; i < sz1[0]*sz1[1]; i++)
				result[i] = array1[i] + array2[0];
		}
		else if ( sz1[0] == 1 && sz1[1] == 1 )
		{

			#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
			mexPrintf("Scalar-Matrix Addition!\n");
			#endif

			if (plhs[0] == NULL)
				plhs[0] = mxCreateNumericMatrix(sz2[0],sz2[1], mxGetClassID(prhs[0]),mxREAL);
			T *result = (T *)mxGetData(plhs[0]);
			T *array1 = (T *)mxGetData(prhs[0]);
			T *array2 = (T *)mxGetData(prhs[1]);
			for (uint64_t i = 0; i < sz2[0]*sz2[1]; i++)
				result[i] = array1[0] + array2[i];
		}
		else
		{
			mexPrintf("Error of dismatched sizes: It is none of below substraction operations: \n");
			mexPrintf("Matrix-Matrix Addition\n");
			mexPrintf("Matrix-Vector Addition\n");
			mexPrintf("Matrix-Scalar Addition\n");
			mexErrMsgTxt("");
		}
	}
}
//-------------------------------------------------------------------------------------------------
//	Function mexSubstract: substract two matrices element-wisely. (Untested yet)
//-------------------------------------------------------------------------------------------------
//	Inputs:		4
//-------------------------------------------------------------------------------------------------
//	nlsh: output number;
//	plhs: output pointer-array;
//	nrhs: input number;
//	prhs: input	pointer-array;
//-------------------------------------------------------------------------------------------------
//	OutPuts: 	0
//-------------------------------------------------------------------------------------------------
// Notice: supports 4 verisons: scalar-matrix, matrix-scalar, scalar-scalar, matrix-matrix.
//-------------------------------------------------------------------------------------------------
template<typename T>
void mexSubstract(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 2)
		mexErrMsgTxt("2 inputs are required!\n");
	if (nlhs != 1)
		mexErrMsgTxt("1 output is required!\n");
	uint64_t sz1[2],sz2[2];
	if (mxGetNumberOfDimensions(prhs[0]) !=2 )
		mexErrMsgTxt("Only 2-D matrix substraction is supported!");
	sz1[0] = mxGetM(prhs[0]); sz1[1] = mxGetN(prhs[0]);
	sz2[0] = mxGetM(prhs[1]); sz2[1] = mxGetN(prhs[1]);
	if (sz1[0] == sz2[0] && sz1[1] == sz2[1])
	{

		#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
		mexPrintf("2-D Matrices Substraction!\n");
		#endif

		if (plhs[0] == NULL)
			plhs[0] = mxCreateNumericMatrix(sz1[0],sz1[1], mxGetClassID(prhs[0]),mxREAL);
		T *result = (T *)mxGetData(plhs[0]);
		T *array1 = (T *)mxGetData(prhs[0]);
		T *array2 = (T *)mxGetData(prhs[1]);
		for (uint64_t i = 0; i < sz1[0]*sz1[1]; i++)
			result[i] = array1[i] - array2[i];
	}
	else
	{

		#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
		mexPrintf("Substraction:\t");
		#endif

		if (sz1[0] == sz2[0] && sz2[0] > 1 && sz2[1] == 1)
		{

			#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
			mexPrintf("Matrix-ColVector substraction!\n");
			#endif

			if (plhs[0] == NULL)
				plhs[0] = mxCreateNumericMatrix(sz1[0],sz1[1], mxGetClassID(prhs[0]),mxREAL);
			T *result = (T *)mxGetData(plhs[0]);
			T *array1 = (T *)mxGetData(prhs[0]);
			T *array2 = (T *)mxGetData(prhs[1]);
			uint64_t step = sz1[0];
			for (uint64_t col = 0; col < sz1[1]; col++)
				for (uint64_t row = 0; row < sz1[0]; row++)
					result[col * step + row] = array1[col * step + row] - array2[row];
		}

		else if (sz1[1] == sz2[1] && sz2[0] == 1 && sz2[1] > 1)
		{

			#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
			mexPrintf("Matrix-RowVector substraction!\n");
			#endif

			if (plhs[0] == NULL)
				plhs[0] = mxCreateNumericMatrix(sz1[0],sz1[1], mxGetClassID(prhs[0]),mxREAL);
			T *result = (T *)mxGetData(plhs[0]);
			T *array1 = (T *)mxGetData(prhs[0]);
			T *array2 = (T *)mxGetData(prhs[1]);
			uint64_t step = sz1[0];
			for (uint64_t col = 0; col < sz1[1]; col++)
				for (uint64_t row = 0; row < sz1[0]; row++)
					result[col * step + row] = array1[col * step + row] - array2[col];
		}
		else if ( sz2[0] == 1 && sz2[1] == 1 )
		{

			#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
			mexPrintf("Matrix-Scalar substraction!\n");
			#endif

			if (plhs[0] == NULL)
				plhs[0] = mxCreateNumericMatrix(sz1[0],sz1[1], mxGetClassID(prhs[0]),mxREAL);
			T *result = (T *)mxGetData(plhs[0]);
			T *array1 = (T *)mxGetData(prhs[0]);
			T *array2 = (T *)mxGetData(prhs[1]);
			for (uint64_t i = 0; i < sz1[0]*sz1[1]; i++)
				result[i] = array1[i] - array2[0];
		}
		else if ( sz1[0] == 1 && sz1[1] == 1 )
		{

			#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
			mexPrintf("Scalar-Matrix substraction!\n");
			#endif

			if (plhs[0] == NULL)
				plhs[0] = mxCreateNumericMatrix(sz2[0],sz2[1], mxGetClassID(prhs[0]),mxREAL);
			T *result = (T *)mxGetData(plhs[0]);
			T *array1 = (T *)mxGetData(prhs[0]);
			T *array2 = (T *)mxGetData(prhs[1]);
			for (uint64_t i = 0; i < sz2[0]*sz2[1]; i++)
				result[i] = array1[0] - array2[i];
		}
		else
		{
			mexPrintf("Error of dismatched sizes: It is none of below substraction operations: \n");
			mexPrintf("Matrix-Matrix Substraction\n");
			mexPrintf("Matrix-Vector Substraction\n");
			mexPrintf("Matrix-Scalar Substraction\n");
			mexErrMsgTxt("");
		}
	}
}
//-------------------------------------------------------------------------------------------------
//	Function mexDotProduction: product two matrices element-wisely. (Untested yet)
//-------------------------------------------------------------------------------------------------
//	Inputs:		4
//-------------------------------------------------------------------------------------------------
//	nlsh: output number;
//	plhs: output pointer-array;
//	nrhs: input number;
//	prhs: input	pointer-array;
//-------------------------------------------------------------------------------------------------
//	OutPuts: 	0
//-------------------------------------------------------------------------------------------------
// Notice: supports 4 verisons: scalar-matrix, matrix-scalar, scalar-scalar, matrix-matrix.
//-------------------------------------------------------------------------------------------------
template<typename T>
void mexDotProduction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 2)
		mexErrMsgTxt("2 inputs are required!\n");
	if (nlhs != 1)
		mexErrMsgTxt("1 output is required!\n");
	uint64_t sz1[2],sz2[2];
	if(mxGetNumberOfDimensions(prhs[0]) != mxGetNumberOfDimensions(prhs[1]))
		mexErrMsgTxt("the dimensions of 2 matrces is not equial!\n ");
	sz1[0] = mxGetM(prhs[0]); sz1[1] = mxGetN(prhs[0]);
	sz2[0] = mxGetM(prhs[1]); sz2[1] = mxGetN(prhs[1]);
	if (sz1[0] == sz2[0] && sz1[1] == sz2[1])
	{

		#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
		mexPrintf("Matrix-Matrix element-wise production!\n");
		#endif

		if (plhs[0] == NULL)
			plhs[0] = mxCreateNumericMatrix(sz1[0],sz1[1], mxGetClassID(prhs[0]),mxREAL);
		T *result = (T *)mxGetData(plhs[0]);
		T *array1 = (T *)mxGetData(prhs[0]);
		T *array2 = (T *)mxGetData(prhs[1]);
		for (uint64_t i = 0; i < sz1[0]*sz1[1]; i++)
			result[i] = array1[i] * array2[i];
	}
	else if(sz2[0] == 1 && sz2[1] == 1 )
	{

		#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
		mexPrintf("Matrix-scalar element-wise production!\n");
		#endif

		if (plhs[0] == NULL)
			plhs[0] = mxCreateNumericMatrix(sz1[0],sz1[1], mxGetClassID(prhs[0]),mxREAL);
		T *result = (T *)mxGetData(plhs[0]);
		T *array1 = (T *)mxGetData(prhs[0]);
		T *array2 = (T *)mxGetData(prhs[1]);
		for (uint64_t i = 0; i < sz1[0]*sz1[1]; i++)
			result[i] = array1[i] * array2[0];
	}
	else if(sz1[0] == 1 && sz1[1] == 1 )
	{

		#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
		mexPrintf("Scalar-matrix element-wise production!\n");
		#endif

		if (plhs[0] == NULL)
			plhs[0] = mxCreateNumericMatrix(sz2[0],sz2[1], mxGetClassID(prhs[0]),mxREAL);
		T *result = (T *)mxGetData(plhs[0]);
		T *array1 = (T *)mxGetData(prhs[0]);
		T *array2 = (T *)mxGetData(prhs[1]);
		for (uint64_t i = 0; i < sz1[0]*sz1[1]; i++)
			result[i] = array1[0] * array2[i];
	}
	else
	{
		mexPrintf("Error of dismatched sizes: It is none of below substraction operations: \n");
		mexPrintf("Matrix-Matrix element-wise production\n");
		mexPrintf("Matrix-Scalar element-wise production\n");
		mexPrintf("Scalar-Matrix element-wise production\n");
		mexErrMsgTxt("");
	}
}
//-------------------------------------------------------------------------------------------------
//	Function mexDotDivision: divide two matrices element-wisely. (Untested yet)
//-------------------------------------------------------------------------------------------------
//	Inputs:		4
//-------------------------------------------------------------------------------------------------
//	nlsh: output number;
//	plhs: output pointer-array;
//	nrhs: input number;
//	prhs: input	pointer-array;
//-------------------------------------------------------------------------------------------------
//	OutPuts: 	0
//-------------------------------------------------------------------------------------------------
// Notices:
// 1. supports 4 verisons: scalar-matrix, matrix-scalar, scalar-scalar, matrix-matrix.
// 2. zero-denominator check is not included!
//-------------------------------------------------------------------------------------------------
template<typename T>
void mexDotDivision(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 2)
		mexErrMsgTxt("2 inputs are required!\n");
	if (nlhs != 1)
		mexErrMsgTxt("1 output is required!\n");
	uint64_t sz1[2],sz2[2];
	if(mxGetNumberOfDimensions(prhs[0]) != mxGetNumberOfDimensions(prhs[1]))
		mexErrMsgTxt("the dimensions of 2 matrces is not equial!\n ");
	sz1[0] = mxGetM(prhs[0]); sz1[1] = mxGetN(prhs[0]);
	sz2[0] = mxGetM(prhs[1]); sz2[1] = mxGetN(prhs[1]);
	if (sz1[0] == sz2[0] && sz1[1] == sz2[1])
	{

		#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
		mexPrintf("Matrix-Matrix element-wise division!\n");
		#endif

		if (plhs[0] == NULL)
			plhs[0] = mxCreateNumericMatrix(sz1[0],sz1[1], mxGetClassID(prhs[0]),mxREAL);
		T *result = (T *)mxGetData(plhs[0]);
		T *array1 = (T *)mxGetData(prhs[0]);
		T *array2 = (T *)mxGetData(prhs[1]);
		for (uint64_t i = 0; i < sz1[0]*sz1[1]; i++)
			result[i] = max(array1[i], array2[i]);
	}
	else if(sz2[0] == 1 && sz2[1] == 1 )
	{

		#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
		mexPrintf("Matrix-scalar element-wise division!\n");
		#endif

		if (plhs[0] == NULL)
			plhs[0] = mxCreateNumericMatrix(sz1[0],sz1[1], mxGetClassID(prhs[0]),mxREAL);
		T *result = (T *)mxGetData(plhs[0]);
		T *array1 = (T *)mxGetData(prhs[0]);
		T *array2 = (T *)mxGetData(prhs[1]);
		if (array2[0] != 0)
		{
			for (uint64_t i = 0; i < sz1[0]*sz1[1]; i++)
				result[i] = array1[i] / array2[0];
		}
		else
		{
			for (uint64_t i = 0; i < sz1[0]*sz1[1]; i++)
				result[i] = mxGetNaN();
		}
	}
	else if(sz1[0] == 1 && sz1[1] == 1 )
	{

		#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
		mexPrintf("Scalar-matrix element-wise division!\n");
		#endif

		if (plhs[0] == NULL)
			plhs[0] = mxCreateNumericMatrix(sz2[0],sz2[1], mxGetClassID(prhs[0]),mxREAL);
		T *result = (T *)mxGetData(plhs[0]);
		T *array1 = (T *)mxGetData(prhs[0]);
		T *array2 = (T *)mxGetData(prhs[1]);
		for (uint64_t i = 0; i < sz1[0]*sz1[1]; i++)
		{
			if (array2[i] == 0)
				result[i] = mxGetNaN();
			else
				result[i] = array1[0] / array2[i];
		}
	}
	else
	{
		mexPrintf("Error of dismatched sizes: It is none of below substraction operations: \n");
		mexPrintf("Matrix-Matrix element-wise division\n");
		mexPrintf("Matrix-Scalar element-wise division\n");
		mexPrintf("Scalar-Matrix element-wise division\n");
		mexErrMsgTxt("");
	}
}
//-------------------------------------------------------------------------------------------------
//	Function mexDotMax: choose the max between matrices element-wisely. (Untested yet)
//-------------------------------------------------------------------------------------------------
//	Inputs:		4
//-------------------------------------------------------------------------------------------------
//	nlsh: output number;
//	plhs: output pointer-array;
//	nrhs: input number;
//	prhs: input	pointer-array;
//-------------------------------------------------------------------------------------------------
//	OutPuts: 	0
//-------------------------------------------------------------------------------------------------
// Notices:
// 1. supports 4 verisons: scalar-matrix, matrix-scalar, scalar-scalar, matrix-matrix.
// 2. zero-denominator check is not included!
//-------------------------------------------------------------------------------------------------
template<typename T>
void mexDotMax(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 2)
		mexErrMsgTxt("2 inputs are required!\n");
	if (nlhs != 1)
		mexErrMsgTxt("1 output is required!\n");
	uint64_t sz1[2],sz2[2];
	if(mxGetNumberOfDimensions(prhs[0]) != mxGetNumberOfDimensions(prhs[1]))
		mexErrMsgTxt("the dimensions of 2 matrces is not equial!\n ");
	sz1[0] = mxGetM(prhs[0]); sz1[1] = mxGetN(prhs[0]);
	sz2[0] = mxGetM(prhs[1]); sz2[1] = mxGetN(prhs[1]);
	if (sz1[0] == sz2[0] && sz1[1] == sz2[1])
	{

		#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
		mexPrintf("Matrix-Matrix element-wise max!\n");
		#endif

		if (plhs[0] == NULL)
			plhs[0] = mxCreateNumericMatrix(sz1[0],sz1[1], mxGetClassID(prhs[0]),mxREAL);
		T *result = (T *)mxGetData(plhs[0]);
		T *array1 = (T *)mxGetData(prhs[0]);
		T *array2 = (T *)mxGetData(prhs[1]);
		for (uint64_t i = 0; i < sz1[0]*sz1[1]; i++)
			result[i] = max(array1[i],array2[i]);
	}
	else if(sz2[0] == 1 && sz2[1] == 1 )
	{

		#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
		mexPrintf("Matrix-scalar element-wise max!\n");
		#endif

		if (plhs[0] == NULL)
			plhs[0] = mxCreateNumericMatrix(sz1[0],sz1[1], mxGetClassID(prhs[0]),mxREAL);
		T *result = (T *)mxGetData(plhs[0]);
		T *array1 = (T *)mxGetData(prhs[0]);
		T *array2 = (T *)mxGetData(prhs[1]);
		for (uint64_t i = 0; i < sz1[0]*sz1[1]; i++)
			result[i] = max(array1[i],array2[0]);
	}
	else if(sz1[0] == 1 && sz1[1] == 1 )
	{

		#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
		mexPrintf("Scalar-matrix element-wise max!\n");
		#endif

		if (plhs[0] == NULL)
			plhs[0] = mxCreateNumericMatrix(sz2[0],sz2[1], mxGetClassID(prhs[0]),mxREAL);
		T *result = (T *)mxGetData(plhs[0]);
		T *array1 = (T *)mxGetData(prhs[0]);
		T *array2 = (T *)mxGetData(prhs[1]);
		for (uint64_t i = 0; i < sz1[0]*sz1[1]; i++)
			result[i] = max(array1[0],array2[i]);
	}
	else
	{
		mexPrintf("Error of dismatched sizes: It is none of below substraction operations: \n");
		mexPrintf("Matrix-Matrix element-wise production\n");
		mexPrintf("Matrix-Scalar element-wise production\n");
		mexPrintf("Scalar-Matrix element-wise production\n");
		mexErrMsgTxt("");
	}
}
//-------------------------------------------------------------------------------------------------
//	Function mexDotMin: choose the max between matrices element-wisely. (Untested yet)
//-------------------------------------------------------------------------------------------------
//	Inputs:		4
//-------------------------------------------------------------------------------------------------
//	nlsh: output number;
//	plhs: output pointer-array;
//	nrhs: input number;
//	prhs: input	pointer-array;
//-------------------------------------------------------------------------------------------------
//	OutPuts: 	0
//-------------------------------------------------------------------------------------------------
// Notices:
// 1. supports 4 verisons: scalar-matrix, matrix-scalar, scalar-scalar, matrix-matrix.
// 2. zero-denominator check is not included!
//-------------------------------------------------------------------------------------------------
template<typename T>
void mexDotMin(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 2)
		mexErrMsgTxt("2 inputs are required!\n");
	if (nlhs != 1)
		mexErrMsgTxt("1 output is required!\n");
	uint64_t sz1[2],sz2[2];
	if(mxGetNumberOfDimensions(prhs[0]) != mxGetNumberOfDimensions(prhs[1]))
		mexErrMsgTxt("the dimensions of 2 matrces is not equial!\n ");
	sz1[0] = mxGetM(prhs[0]); sz1[1] = mxGetN(prhs[0]);
	sz2[0] = mxGetM(prhs[1]); sz2[1] = mxGetN(prhs[1]);
	if (sz1[0] == sz2[0] && sz1[1] == sz2[1])
	{

		#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
		mexPrintf("Matrix-Matrix element-wise min!\n");
		#endif

		if (plhs[0] == NULL)
			plhs[0] = mxCreateNumericMatrix(sz1[0],sz1[1], mxGetClassID(prhs[0]),mxREAL);
		T *result = (T *)mxGetData(plhs[0]);
		T *array1 = (T *)mxGetData(prhs[0]);
		T *array2 = (T *)mxGetData(prhs[1]);
		for (uint64_t i = 0; i < sz1[0]*sz1[1]; i++)
			result[i] = min(array1[i],array2[i]);
	}
	else if(sz2[0] == 1 && sz2[1] == 1 )
	{

		#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
		mexPrintf("Matrix-scalar element-wise min!\n");
		#endif

		if (plhs[0] == NULL)
			plhs[0] = mxCreateNumericMatrix(sz1[0],sz1[1], mxGetClassID(prhs[0]),mxREAL);
		T *result = (T *)mxGetData(plhs[0]);
		T *array1 = (T *)mxGetData(prhs[0]);
		T *array2 = (T *)mxGetData(prhs[1]);
		for (uint64_t i = 0; i < sz1[0]*sz1[1]; i++)
			result[i] = min(array1[i],array2[0]);
	}
	else if(sz1[0] == 1 && sz1[1] == 1 )
	{

		#ifdef SRC_USEFULMEXFUNCTIONS_HPP_DEBUG
		mexPrintf("Scalar-matrix element-wise min!\n");
		#endif

		if (plhs[0] == NULL)
			plhs[0] = mxCreateNumericMatrix(sz2[0],sz2[1], mxGetClassID(prhs[0]),mxREAL);
		T *result = (T *)mxGetData(plhs[0]);
		T *array1 = (T *)mxGetData(prhs[0]);
		T *array2 = (T *)mxGetData(prhs[1]);
		for (uint64_t i = 0; i < sz1[0]*sz1[1]; i++)
			result[i] = min(array1[0],array2[i]);
	}
	else
	{
		mexPrintf("Error of dismatched sizes: It is none of below substraction operations: \n");
		mexPrintf("Matrix-Matrix element-wise production\n");
		mexPrintf("Matrix-Scalar element-wise production\n");
		mexPrintf("Scalar-Matrix element-wise production\n");
		mexErrMsgTxt("");
	}
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// 4. Matrix index operaton:
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//	Function mexExtractRow: Extract a row/rows from vector/matrix. (Untested yet)
//-------------------------------------------------------------------------------------------------
//	Inputs:		4
//-------------------------------------------------------------------------------------------------
//	nlsh: output number;
//	plhs: output pointer-array;
//	nrhs: input number;
//	prhs: input	pointer-array;
//-------------------------------------------------------------------------------------------------
//	OutPuts: 	0
//-------------------------------------------------------------------------------------------------
// Notices: multiple rows extraction has not been supported yet!
//-------------------------------------------------------------------------------------------------
template<typename T, typename Pos>
void mexExtractRow(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 2)
		mexErrMsgTxt("2 inputs are required!\n");
	if (nlhs != 1)
		mexErrMsgTxt("1 output is required!\n");
	if (mxGetNumberOfDimensions(prhs[0]) != 2)
		mexErrMsgTxt("2-D matrices are expacted!\n");
	uint64_t rows, cols;
	rows = mxGetM(prhs[0]); cols = mxGetN(prhs[0]);
	if (rows < 1 || cols < 1)
		mexErrMsgTxt("The input matrix is empty!\n");
	T *row = (T *)mxGetData(prhs[1]);
	Pos nDims = mxGetM(prhs[1]) < mxGetN(prhs[1]) ? mxGetN(prhs[1]) : mxGetM(prhs[1]);
	if (rows < nDims) mexErrMsgTxt("Indices exceeds dimensions of matrix!\n");
	if (rows > 1)
	{
		if (plhs[0] == NULL)
			plhs[0] = mxCreateNumericMatrix(nDims,cols, mxGetClassID(prhs[0]), mxREAL);
		T *result = (T *)mxGetData(plhs[0]);
		T *array = (T *)mxGetData(prhs[0]);
		for (uint64_t i = 0; i < nDims; i++)
			for (uint64_t col = 0; col < cols; col++)
				result[i + col * rows] = array[(uint64_t)row[i] + col * rows];
	}
	else
	{
		if (plhs[0] != NULL)
			mxDestroyArray(plhs[0]);
		plhs[0] = mxDuplicateArray(prhs[0]);
	}
}
//-------------------------------------------------------------------------------------------------
//	Function mexExtractColumn: Extract a column/columns from vector/matrix. (Untested yet)
//-------------------------------------------------------------------------------------------------
//	Inputs:		4
//-------------------------------------------------------------------------------------------------
//	nlsh: output number;
//	plhs: output pointer-array;
//	nrhs: input number;
//	prhs: input	pointer-array;
//-------------------------------------------------------------------------------------------------
//	OutPuts: 	0
//-------------------------------------------------------------------------------------------------
// Notices: multiple cols extraction has not been supported yet!
//-------------------------------------------------------------------------------------------------
template<typename T, typename Pos>
void mexExtractColumn(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 2)
		mexErrMsgTxt("2 inputs are required!\n");
	if (nlhs != 1)
		mexErrMsgTxt("1 output is required!\n");
	if (mxGetNumberOfDimensions(prhs[0]) != 2)
		mexErrMsgTxt("2-D matrices are expacted!\n");
	uint64_t rows, cols;
	rows = mxGetM(prhs[0]); cols = mxGetN(prhs[0]);
	if (rows < 1 || cols < 1)
		mexErrMsgTxt("The input matrix is empty!\n");
	Pos *col = (Pos *)mxGetData(prhs[1]);
	Pos nDims = mxGetM(prhs[1]) < mxGetN(prhs[1]) ? mxGetN(prhs[1]) : mxGetM(prhs[1]);
	if (cols < nDims) mexErrMsgTxt("Indices exceeds dimensions of matrix!\n");
	if (cols > 1)
	{
		if (plhs[0] == NULL)
			plhs[0] = mxCreateNumericMatrix(rows,1, mxGetClassID(prhs[0]),mxREAL);
		T *result = (T *)mxGetData(plhs[0]);
		T *array = (T *)mxGetData(prhs[0]);
		for (uint64_t i = 0; i < nDims; i++)
			for (uint64_t row = 0; row < rows; row++)
				result[row] = array[row + (uint64_t)col[i] * rows];
	}
	else
	{
		if (plhs[0] != NULL)
			mxDestroyArray(plhs[0]);
		plhs[0] = mxDuplicateArray(prhs[0]);
	}
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// 5. Memory management functions:
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//	Function mexAllocateMemoryWithInitials: denote memory with an intial value for a matrix. (Untested yet)
//-------------------------------------------------------------------------------------------------
//	Inputs:		4
//-------------------------------------------------------------------------------------------------
//	nlsh: output number;
//	plhs: output pointer-array;
//	nrhs: input number;
//	prhs: input	pointer-array;
//-------------------------------------------------------------------------------------------------
//	OutPuts: 	0
//-------------------------------------------------------------------------------------------------
template<typename T>
void mexAllocateMemoryWithInitials(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs !=2 )
		mexErrMsgTxt("2 inputs are required!\n");
	if (nlhs != 1)
		mexErrMsgTxt("1 input is required!\n");
	uint64_t *sz = (uint64_t *)mxGetData(prhs[0]);
	uint64_t rows = sz[0];
	uint64_t cols = 0;
	if (mxGetN(prhs[0]) == 1)
		cols = 1;
	else
		cols = sz[1];
	T *val = (T *)mxGetData(prhs[1]);
	if (plhs[0] != NULL)
		mxDestroyArray(plhs[0]);
	plhs[0] = mxCreateNumericMatrix(rows,cols,mxGetClassID(prhs[1]),mxREAL);
	T *result = (T *)mxGetData(plhs[0]);
	for (uint64_t i = 0; i < rows*cols; i ++)
		result[i] = val[0];
}
//-------------------------------------------------------------------------------------------------
//	Function mexResetMemory: reset memory with a specific value for a matrix. (Untested yet)
//-------------------------------------------------------------------------------------------------
//	Inputs:		4
//-------------------------------------------------------------------------------------------------
//	nlsh: output number;
//	plhs: output pointer-array;
//	nrhs: input number;
//	prhs: input	pointer-array;
//-------------------------------------------------------------------------------------------------
//	OutPuts: 	0
//-------------------------------------------------------------------------------------------------
template<typename T>
void mexResetMemory(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs !=2 )
		mexErrMsgTxt("2 inputs are required!\n");
	if (nlhs !=0 )
		mexErrMsgTxt("No input is required!\n");
	uint64_t sz[2];
	sz[0] = mxGetM(prhs[0]); sz[1] = mxGetN(prhs[0]);
	T *array = (T *)mxGetData(prhs[0]);
	T *val = (T *)mxGetData(prhs[1]);
	for (uint64_t i = 0; i < sz[0]*sz[1]; i ++)
		array[i] = val[0];
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// 4. Miscellaneous utility functions:
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
#endif /* SRC_BASICMEXMATRIXOPERATIONS_HPP_ */

/*
 * Utilities.hpp
 *
 *  Created on: Feb 14, 2017
 *      Author: xiong
 */

#ifndef SRC_MEXUTILITIES_HPP_
#define SRC_MEXUTILITIES_HPP_

#include "BasicMexMatrixOperations.hpp"


//-------------------------------------------------------------------------------------------------
// template<typename T>
// void mexSparseMatrix(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// {
// 	if (nrhs != 5)
// 		mexErrMsgTxt(Error_INPUT_NUMBER);
// 	if (nlhs != 1)
// 		mexErrMsgTxt(Error_OUTPUT_NUMBER);
// 	const mxArray *sp_i = prhs[0];
// 	const mxArray *sp_j = prhs[1];
// 	const mxArray *sp_v = prhs[2];
// 	mwSize rows = (mwSize)((double *)mxGetData(prhs[3]))[0];
// 	mwSize cols = (mwSize)((double *)mxGetData(prhs[4]))[0];
// 	mwSize nzmax = mxGetM(sp_i);
// 	if (plhs[0] != NULL) mxDestroyArray(plhs[0]);
// 	plhs[0] = mxCreateSparse(rows,cols,nzmax,mxREAL);
// 	mwIndex *ir = CDataCast<T,mwIndex,mwSize>((T *)mxGetData(sp_i),nzmax);
// 	mxSetIr(plhs[0],ir);
// 	mwIndex *ic = CDataCast<T,mwIndex,mwSize>((T *)mxGetData(sp_j),nzmax);
// 	mxSetJc(plhs[0],ic);
// 	double *pr = CDataCast<T,double,mwSize>((T *)mxGetData(sp_v),nzmax);
// 	mxSetPr(plhs[0],pr);
// }

//-------------------------------------------------------------------------------------------------
template<typename T>
void mexSplitSum(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 3)
		mexErrMsgTxt(Error_INPUT_NUMBER);
	if (nlhs != 1)
		mexErrMsgTxt(Error_OUTPUT_NUMBER);
	const mxArray *Labels = prhs[0];
	const mxArray *Indices = prhs[1];
	const mxArray *Weights = prhs[2];
	double *labels = (double *)mxGetData(Labels);
	uint64_t *indices = (uint64_t *)mxGetData(Indices);
	double *weights = (double *)mxGetData(Weights);
	uint64_t cols = mxGetN(Labels);
	uint64_t rows = mxGetM(Labels);

	if (plhs[0] == NULL)
		plhs[0] = mxCreateNumericMatrix(rows,cols,mxDOUBLE_CLASS,mxREAL);
	uint64_t *slice = (uint64_t *)mxGetData(plhs[0]);
//	for (uint64_t col = 0; col < cols; col++)
//		for (uint64_t row = 0; row < rows; row++)
//			slice[row + col * rows] += labels + (uint64_t)bin[col] - 1) * (uint64_t)gridSize[col];
}
//-------------------------------------------------------------------------------------------------
void mexComputeIndices(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 3)
		mexErrMsgTxt(Error_INPUT_NUMBER);
	if (nlhs != 1)
		mexErrMsgTxt(Error_OUTPUT_NUMBER);
	const mxArray *Floors = prhs[0];
	const mxArray *Bin = prhs[1];
	const mxArray *GridSize = prhs[2];
	double *floor = (double *)mxGetData(Floors);
	uint64_t *bin = (uint64_t *)mxGetData(Bin);
	double *gridSize = (double *)mxGetData(GridSize);
	uint64_t cols = mxGetN(Floors);
	uint64_t rows = mxGetM(Floors);
	if (plhs[0] == NULL)
		plhs[0] = mxCreateNumericMatrix(rows,1,mxUINT64_CLASS,mxREAL);
	uint64_t *indices = (uint64_t *)mxGetData(plhs[0]);
	for (uint64_t row = 0; row < rows; row++)
		for (uint64_t col = 0; col < cols; col++)
			indices[row] += (uint64_t)((floor[row + col*rows] + bin[col] - 1) * gridSize[col]);
}
//-------------------------------------------------------------------------------------------------
template<typename T>
void mexFactorialDotProduct(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 1)
		mexErrMsgTxt(Error_INPUT_NUMBER);
	if (nlhs != 1)
		mexErrMsgTxt(Error_OUTPUT_NUMBER);
	uint64_t rows,cols, length;
	rows = mxGetM(prhs[0]);
	cols = mxGetN(prhs[0]);
	if (rows != 1 && cols != 1)
		mexErrMsgTxt("Input is neither row-vector nor col-vector!\n");
	length = cols;
	if (rows > 1)
		length = rows;
	if (plhs[0] == NULL)
		plhs[0] = mxCreateNumericMatrix(rows,cols,mxGetClassID(prhs[0]),mxREAL);
	T *Size = (T *)mxGetData(prhs[0]);
	T *ProdSize = (T *)mxGetData(plhs[0]);;
	CFactorialDotProduct<T,uint64_t>(Size,length,ProdSize);
}
//-------------------------------------------------------------------------------------------------
//	Function mxDec2Bin: transfer a decimal integer into a binary array. (Untested yet)
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
// Notes: The output binary array is not a string!
//-------------------------------------------------------------------------------------------------
template<typename T>
void mexDec2Bin(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 2)
		mexErrMsgTxt("2 inputs are required!\n");
	if (nlhs != 1)
		mexErrMsgTxt("1 output is required!\n");
	T *dec = (T *)mxGetData(prhs[0]);
	uint64_t *ndims = (uint64_t *)mxGetData(prhs[1]);
//	uint8_t *format = (uint8_t *)mxGetData(prhs[2]);
	if (plhs[0] == NULL)
		plhs[0] = mxCreateNumericMatrix(1,ndims[0],mxUINT64_CLASS,mxREAL);
	uint64_t *container = (uint64_t* )mxGetData(plhs[0]);
	if ( container == NULL )
		mexErrMsgTxt("Memory allocation is failed!\n");
	int offset = 0;
	for ( int64_t c = ndims[0]-1; c >= 0; c-- )
   {
	  int d = *dec >> (T)c;
	  if ( d & 1 )
		 *(container+offset) = 1;
	  else
		 *(container+offset) = 0;
	  offset++;
   }
}
//-------------------------------------------------------------------------------------------------
//	Function mexReverse: transfer a vector into its reverse order. (Untested yet)
//-------------------------------------------------------------------------------------------------
//	Inputs:		1
//-------------------------------------------------------------------------------------------------
//	nlsh: output number;
//	plhs: output pointer-array;
//	nrhs: input number;
//	prhs: input	pointer-array;
//-------------------------------------------------------------------------------------------------
//	OutPuts: 	1
//-------------------------------------------------------------------------------------------------
// Notes: The output binary array is not a string!
//-------------------------------------------------------------------------------------------------
template<typename T>
void mexReverse(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 1)
		mexErrMsgTxt("1 inputs are required!\n");
	if (nlhs != 1)
		mexErrMsgTxt("1 output is required!\n");
	T *array = (T *)mxGetData(prhs[0]);
	if (mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1)
		mexErrMsgTxt("The input is not a vector!\n");
	mwSize rows = mxGetM(prhs[0]);
	mwSize cols = mxGetN(prhs[0]);
	mwSize length = rows*cols;
//	uint8_t *format = (uint8_t *)mxGetData(prhs[2]);
	if (plhs[0] == prhs[0])
	{
		T swap;
		mwSize offset;
		for ( mwSize c = 0; c < length/2; c++ )
		{
			offset = length - 1 - c;
			swap = array[c];
			array[c] = array[offset];
			array[offset] = swap;
		}
	}
	else if (plhs[0] == NULL)
	{
		plhs[0] = mxCreateNumericMatrix(rows,cols, mxGetClassID(prhs[0]), mxREAL);
		T *result = (T* )mxGetData(plhs[0]);
		for ( mwSize c = 0; c < length; c++ )
			result[c] = array[length - 1 - c];
	}
}
//-------------------------------------------------------------------------------------------------
//	Function mexDisplayArray: display a matrix on screen. (Untested yet)
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
void mexDisplayArray(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//	mexPrintf("Display matrix which is row-first:\n");
	if (nrhs != 1 )
		mexErrMsgTxt("1 inputs are required!\n");
	if (nlhs != 0)
		mexErrMsgTxt("No input is required!\n");
	uint64_t rows, cols;
	rows = mxGetM(prhs[0]); cols = mxGetN(prhs[0]);
	T *array = (T *)mxGetData(prhs[0]);
	for (uint64_t row = 0; row < rows; row ++)
	{
		for (uint64_t col = 0; col < cols; col ++)
			std::cout<<array[row + col * rows]<<"\t";
		std::cout<<"\n";
	}
	SplitLine;
}
//-------------------------------------------------------------------------------------------------
#endif /* SRC_MEXUTILITIES_HPP_ */

/*
 * MapOperation.hpp
 *
 *  Created on: Feb 7, 2017
 *      Author: xiong
 */

#ifndef ACCUMULATESUM_HPP_
#define ACCUMULATESUM_HPP_
//-------------------------------------------------------------------------------------------------
// User-defined Macros
//-------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------
//Head files
//-------------------------------------------------------------------------------------------------
#include "BasicMexMatrixOperations.hpp"
#include <map>
//-------------------------------------------------------------------------------------------------
//Functions:
//-------------------------------------------------------------------------------------------------
template<typename T1, typename T2>
void CreateMapByKeysAndValues(mwSize nClass, const mxArray *prhs[],std::map<T1, T2>*Accumulator[])
{
	for (mwSize n = 0; n < nClass; n++)
	{
		if (prhs[n] == NULL)
		{
			mexPrintf("%dth Input exists!\n",n);
			mexErrMsgTxt("");
		}
		uint64_t rows, cols;
		rows = mxGetM(prhs[n]); cols = mxGetN(prhs[n]);
		if (cols != 2)
			mexErrMsgTxt("Input is not a Mx2 matrix which contains indices and values on each on columns!\n");
		T1 *Data = (T1 *)mxGetData(prhs[n]);
		if (Accumulator[n] == NULL)
			Accumulator[n] = new std::map<T1, T2>();
		for( uint64_t row = 0; row < rows; row ++ )
			Accumulator[n][0][Data[row]] = Data[row + rows];
	}
}
//-------------------------------------------------------------------------------------------------
template<typename T1,typename T2>
void CreateMapByKeysAndInitials(int nrhs, const mxArray *prhs[],std::map<T1, T2>*Accumulator[])
{
	if (nrhs != 2)
		mexErrMsgTxt(Error_INPUT_NUMBER);
	T1 *indices = (T1 *)mxGetData(prhs[0]);
	T2 *values = (T2 *)mxGetData(prhs[1]);
	uint64_t rows, cols;
	rows = mxGetM(prhs[0]); cols = mxGetN(prhs[0]);
	for (uint64_t col = 0; col < cols; col ++)
	{
		if (Accumulator[col] == NULL)
				Accumulator[col] = new std::map<T1, T2>();
		for( uint64_t row = 0; row < rows; row ++ )
			Accumulator[col][0][indices[row + col * rows]] = values[0];
	}
}
//-------------------------------------------------------------------------------------------------
template<typename T1, typename T3, typename T4>
void AssignmentByIndicesMap(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[],std::map<T3, T4>*Accumulator[])
{
	if (nrhs != 1)
		mexErrMsgTxt(Error_INPUT_NUMBER);
	if (nlhs != 1)
		mexErrMsgTxt(Error_OUTPUT_NUMBER);
	if (plhs[0] == NULL)
		mexErrMsgTxt("Plhs[0] array is empty\n");
	T1 *indices = (T1 *)mxGetData(prhs[0]);
	uint64_t rows, cols;
	rows= mxGetM(plhs[0]); cols= mxGetN(plhs[0]);
	T1 *Weight = (T1 *)mxGetData(plhs[0]);
	typename std::map<T3, T4>::iterator iter;
	for (uint64_t col = 0; col < cols; col ++)
	{
		if (Accumulator[col] == NULL)
			mexErrMsgTxt("No reference-map exist!\n");
		for( uint64_t row = 0; row < rows; row ++ )
		{
			iter = Accumulator[col]->find(indices[row]);
			if (iter != Accumulator[col]->end())
				Weight[row + col *rows] = iter->second;
			else
				continue;
		}
	}
}
//-------------------------------------------------------------------------------------------------
template<typename T1, typename T2, typename T3, typename T4>
void DotProductionMap(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[],std::map<T3, T4>*Accumulator[])
{
	if (nrhs != 2)
		mexErrMsgTxt(Error_INPUT_NUMBER);
	T1 *indices = (T1 *)mxGetData(prhs[0]);
	T2 *values = (T2 *)mxGetData(prhs[1]);
	uint64_t rows, cols;
	rows= mxGetM(prhs[0]); cols= mxGetN(prhs[0]);
	if (nlhs != 1)
		mexErrMsgTxt(Error_OUTPUT_NUMBER);
	if (mxGetM(prhs[0]) != mxGetM(prhs[1]) || mxGetN(prhs[0]) != mxGetN(prhs[1]))
		mexErrMsgTxt("Index array is not match with value array\n");
	if (plhs[0] == NULL)
		plhs[0] = mxCreateNumericMatrix(rows,cols,mxGetClassID(prhs[1]),mxREAL);
	T2 *result = (T2 *)mxGetData(plhs[0]);
	typename std::map<T3, T4>::iterator iter;
	for (uint64_t col = 0; col < cols; col ++)
	{
		if (Accumulator[col] == NULL)
				mexErrMsgTxt("No reference-map exist!\n");
		for( uint64_t row = 0; row < rows; row ++ )
		{
			iter = Accumulator[col]->find(indices[row + col * rows]);
			if (iter != Accumulator[col]->end())
				result[row + col * rows] = values[row + col * rows]*(iter->second);
			else
				continue;
//				mexErrMsgTxt("No resource-index found in the reference map!\n");
		}
	}
}
//-------------------------------------------------------------------------------------------------
template<typename T1, typename T2>
void SplitMap(int nlhs, mxArray *plhs[], uint64_t nClass, std::map<T1, T2>*Accumulator[])
{
	if (nlhs != 1)
		mexErrMsgTxt(Error_OUTPUT_NUMBER);
	if (plhs[0] == NULL)
		plhs[0] = mxCreateCellMatrix( 1, nClass );
	for (uint64_t n = 0; n < nClass; n ++)
	{
		uint64_t length = Accumulator[n]->size();
		mxArray *cell = mxGetCell(plhs[0],n);
		if (cell == NULL)
		{
			cell = mxCreateNumericMatrix(length,2,mxDOUBLE_CLASS,mxREAL);
			mxSetCell(plhs[0],n,cell);
		}
		double *Data = ( double * )mxGetData( cell );
		uint64_t count = 0;
		for( typename std::map<T1, T2>::iterator iter = Accumulator[n]->begin(); iter != Accumulator[n]->end(); iter++ )
		{
			Data[ count] = (double) iter->first;
			Data[ count + length ] = (double) iter->second;
			count++;
		}
	}
}
//-------------------------------------------------------------------------------------------------
template<typename T1, typename T2, typename T3, typename T4>
void AccMapSum( int nrhs, const mxArray *prhs[], std::map<T3, T4>*Accumulator[] )
{
	if (nrhs != 3)
		mexErrMsgTxt(Error_INPUT_NUMBER);
	const mxArray *Indices = prhs[0];
	const mxArray *Weight = prhs[1];
	const mxArray *Data = prhs[2];
	T1 *indices = (T1 *)mxGetData(Indices);
	T2 *weight = (T2 *)mxGetData(Weight);
	T2 *data = (T2 *)mxGetData(Data);
	mwSize rows = mxGetM(Data);
	mwSize cols = mxGetN(Data);
	typename std::map<T3, T4>::iterator iter;
	for (mwSize col = 0; col < cols; col ++)
	{
		if (Accumulator[col] == NULL)
			Accumulator[col] = new std::map<T3, T4>();
		for (uint64_t row = 0; row < rows; row ++)
		{
			if (weight[row]*data[row + col * rows] > 0)
			{
				Accumulator[col][0][indices[row]] += weight[row]*data[row + col * rows];
//				iter = Accumulator[col][0].find(indices[row]);
//				std::cout<<"\t<"<<col+1<<","<<row+1<<","<<iter->first<<","<<weight[row]*data[row + col * rows]<<">\t"<<"=>";
//				std::cout<<"\t<"<<col+1<<","<<row+1<<","<<iter->first<<","<<iter->second<<">\t"<<"\n";
			}
		}
	}
}
//-------------------------------------------------------------------------------------------------
template<typename T1, typename T2>
void DisplayMap(int nrhs, const mxArray *prhs[],std::map<T1, T2>*Accumulator[])
{
	if (nrhs != 1)
		mexErrMsgTxt(Error_INPUT_NUMBER);
	uint64_t *nClass = (uint64_t*)mxGetData(prhs[0]);
	SplitLine; std::cout<<"Display contents of accumulator:\n";
	for (uint64_t n = 0; n < nClass[0]; n ++)
	{
		std::cout<<"Class "<<n<<":\n";
		uint64_t offset = 0;
		for( typename std::map<T1, T2>::iterator iter = Accumulator[n]->begin(); iter != Accumulator[n]->end(); iter++ )
		{
			std::cout<<"\t<"<<offset<<","<<iter->first<<","<<iter->second<<">\t"<<"\n";
			offset++;
		}
	}
	SplitLine;
}
//-------------------------------------------------------------------------------------------------
#endif /* ACCUMULATESUM_HPP_ */

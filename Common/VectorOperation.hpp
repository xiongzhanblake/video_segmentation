/*
 * VectorOperation.hpp
 *
 *  Created on: May 15, 2017
 *      Author: xiong
 */

#ifndef VECTOROPERATION_HPP_
#define VECTOROPERATION_HPP_

//-------------------------------------------------------------------------------------------------
// User-defined Macros
//-------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------
//Head files
//-------------------------------------------------------------------------------------------------
#include "BasicMexMatrixOperations.hpp"
#include <vector>
//-------------------------------------------------------------------------------------------------
//Functions:
//-------------------------------------------------------------------------------------------------

template<typename T, typename L>
std::vector<T> *InitVector(T *data, L length)
{
	std::vector<T> *ret = new std::vector<T>(data,data + length);
	return ret;
}

//-------------------------------------------------------------------------------------------------
template<typename T, typename Pos>
void mexIntersectVector(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	typedef typename std::vector<T> VectorT;
	typedef typename std::vector<Pos> VectorP;
	typedef typename VectorT::iterator IteratorT;
//	typedef typename VectorP::iterator IteratorP;

	if (nrhs != 3)
		mexErrMsgTxt("3 inputs are required!\n");
	if (nlhs != 3)
		mexErrMsgTxt("3 output is required!\n");
	if (mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1)
		mexErrMsgTxt("1-D vector is expected!\n");
	uint64_t rows1, rows2;
	rows1 = mxGetN(prhs[0]) < mxGetM(prhs[0]) ? mxGetM(prhs[0]) : mxGetN(prhs[0]);
	rows2 = mxGetN(prhs[1]) < mxGetM(prhs[1]) ? mxGetM(prhs[1]) : mxGetN(prhs[1]);;
	if (rows1 < 1 || rows2 < 1)
		mexErrMsgTxt("The input vector is empty!\n");
	VectorT v1 ((T*)mxGetData(prhs[0]),(T*)mxGetData(prhs[0]) + rows1);
	VectorT v2 ((T*)mxGetData(prhs[1]),(T*)mxGetData(prhs[1]) + rows2);
	T lowerBound = ((T*)mxGetData(prhs[2]))[0];
	VectorT r;
	VectorP i1;
	VectorP i2;

	IteratorT first1 = v1.begin();
	IteratorT first2 = v2.begin();
	IteratorT last1 = v1.end();
	IteratorT last2 = v2.end();
//	IteratorT result = r.begin();
//	uint64_t circles = 0;
//	IteratorP index1 = i1.begin();
//	IteratorP index2 = i2.begin();
	while (first1 != last1 && first2 != last2)
	{
		if (*first1 < lowerBound){ first1++; continue; }
		if (*first2 < lowerBound) {first2++;continue;}
		if (*first1 < *first2) ++first1;
		else if (*first2 < *first1) ++first2;
		else 
		{
		  r.push_back(*first1);
		  i1.push_back(Pos(first1 - v1.begin()));
		  i2.push_back(Pos(first2 - v2.begin()));
//		  Pos d1 = (Pos)(first1 - v1.begin());
//		  Pos d2 = (Pos)(first2 - v2.begin());
		  ++first1; ++first2;
//		  double d1 = i1.at(circles); double d2 = i2.at(circles);
//		  circles++;
		}
//		circles++;
	}
	v1.clear();
	v2.clear();
	if(plhs[0] != NULL) mxDestroyArray(plhs[0]);
	plhs[0] = mxCreateNumericMatrix(i2.size(),1,mxGetClassID(prhs[0]),mxREAL);
	for (typename std::vector<T>::size_type i = 0; i < i2.size(); i++)
		((Pos *)mxGetData(plhs[0]))[i] = (Pos)i2.at(i);
	if(plhs[1] != NULL) mxDestroyArray(plhs[1]);
	plhs[1] = mxCreateNumericMatrix(i1.size(),1,mxGetClassID(prhs[1]),mxREAL);
	for (typename std::vector<T>::size_type i = 0; i < i1.size(); i++)
		((Pos *)mxGetData(plhs[1]))[i] = (Pos)i1.at(i);
	if(plhs[2] != NULL) mxDestroyArray(plhs[2]);
	plhs[2] = mxCreateNumericMatrix(r.size(),1,mxGetClassID(prhs[2]),mxREAL);
	for (typename std::vector<T>::size_type i = 0; i < r.size(); i++)
		((T *)mxGetData(plhs[2]))[i] = (T)r.at(i);
	r.clear();
	i1.clear();
	i2.clear();
}
////-------------------------------------------------------------------------------------------------
//// a function to merge two 2-D matrix by element order of rows
////-------------------------------------------------------------------------------------------------
//template<typename T, typename Pos>
//void mexOrderedMergeByRow(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//{
//	typedef typename std::vector<T> VectorT;
//	typedef typename std::vector<Pos> VectorP;
//	typedef typename VectorT::iterator IteratorT;
////	typedef typename VectorP::iterator IteratorP;
//
//	if (nrhs != 3)
//		mexErrMsgTxt("3 inputs are required!\n");
//	if (nlhs != 3)
//		mexErrMsgTxt("3 output is required!\n");
//	if (mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1)
//		mexErrMsgTxt("1-D vector is expected!\n");
//	uint64_t rows1, rows2;
//	rows1 = mxGetN(prhs[0]) < mxGetM(prhs[0]) ? mxGetM(prhs[0]) : mxGetN(prhs[0]);
//	rows2 = mxGetN(prhs[1]) < mxGetM(prhs[1]) ? mxGetM(prhs[1]) : mxGetN(prhs[1]);;
//	if (rows1 < 1 || rows2 < 1)
//		mexErrMsgTxt("The input vector is empty!\n");
//	VectorT v1 ((T*)mxGetData(prhs[0]),(T*)mxGetData(prhs[0]) + rows1);
//	VectorT v2 ((T*)mxGetData(prhs[1]),(T*)mxGetData(prhs[1]) + rows2);
//	T lowerBound = ((T*)mxGetData(prhs[2]))[0];
//	VectorT r;
//	VectorP i1;
//	VectorP i2;
//
//	IteratorT first1 = v1.begin();
//	IteratorT first2 = v2.begin();
//	IteratorT last1 = v1.end();
//	IteratorT last2 = v2.end();
////	IteratorT result = r.begin();
////	uint64_t circles = 0;
////	IteratorP index1 = i1.begin();
////	IteratorP index2 = i2.begin();
//	while (first1 != last1 && first2 != last2)
//	{
//		if (*first1 < lowerBound){ first1++; continue; }
//		if (*first2 < lowerBound) {first2++;continue;}
//		if (*first1 < *first2) ++first1;
//		else if (*first2 < *first1) ++first2;
//		else
//		{
//		  r.push_back(*first1);
//		  i1.push_back(Pos(first1 - v1.begin()));
//		  i2.push_back(Pos(first2 - v2.begin()));
////		  Pos d1 = (Pos)(first1 - v1.begin());
////		  Pos d2 = (Pos)(first2 - v2.begin());
//		  ++first1; ++first2;
////		  double d1 = i1.at(circles); double d2 = i2.at(circles);
////		  circles++;
//		}
////		circles++;
//	}
//	v1.clear();
//	v2.clear();
//	if(plhs[0] != NULL) mxDestroyArray(plhs[0]);
//	plhs[0] = mxCreateNumericMatrix(i2.size(),1,mxGetClassID(prhs[0]),mxREAL);
//	for (typename std::vector<T>::size_type i = 0; i < i2.size(); i++)
//		((Pos *)mxGetData(plhs[0]))[i] = (Pos)i2.at(i);
//	if(plhs[1] != NULL) mxDestroyArray(plhs[1]);
//	plhs[1] = mxCreateNumericMatrix(i1.size(),1,mxGetClassID(prhs[1]),mxREAL);
//	for (typename std::vector<T>::size_type i = 0; i < i1.size(); i++)
//		((Pos *)mxGetData(plhs[1]))[i] = (Pos)i1.at(i);
//	if(plhs[2] != NULL) mxDestroyArray(plhs[2]);
//	plhs[2] = mxCreateNumericMatrix(r.size(),1,mxGetClassID(prhs[2]),mxREAL);
//	for (typename std::vector<T>::size_type i = 0; i < r.size(); i++)
//		((T *)mxGetData(plhs[2]))[i] = (T)r.at(i);
//	r.clear();
//	i1.clear();
//	i2.clear();
//}
//-------------------------------------------------------------------------------------------------
#endif

/*
 * CUtilities.hpp
 *
 *  Created on: Feb 15, 2017
 *      Author: xiong
 */

#ifndef SRC_CUTILITIES_HPP_
#define SRC_CUTILITIES_HPP_

#include <stdint.h>
#include <iostream>
#include <math.h>

using namespace std;

template<typename T, typename T2, typename T3>
T3 CFactorialProduct(const T *array, T2 length)
{
	if (array == NULL)
	{
		printf("Input is empty!\n");
		return -1;
	}
	if (length <= 0)
	{
		printf("length is zero or negative!\n");
		return -1;
	}
	T3 ret = 1;
	for (T2 i = 0; i < length; i++)
		ret *= array[i];
	return ret;
}

template<typename T, typename T2>
T CMin(T *array, T2 length)
{
	if (array == NULL)
	{
		printf("Input is empty!\n");
		return -1;
	}
	if (length <= 0)
	{
		printf("length is zero or negative!\n");
		return -1;
	}
	T ret = array[0];
	for (int i = 0; i < length; i++)
		if (ret > array[i])
			ret = array[i];
	return ret;
}

template<typename T, typename T2>
T* CMin(T *array1, T *array2, T2 length)
{
	if (array1 == NULL || array2 == NULL)
	{
		printf("Input is empty!\n");
		return -1;
	}
	if (length <= 0)
	{
		printf("length is zero or negative!\n");
		return -1;
	}
	T *ret = new T[length];
	for (int i = 0; i < length; i++)
	{
		if (array1[i] < array2[i])
			ret[i] = array1[i];
		else
			ret[i] = array2[i];
	}
	return ret;
}

template<typename T, typename T2>
T CMax(T *array, T2 length)
{
	if (array == NULL)
	{
		printf("Input is empty!\n");
		return -1;
	}
	if (length <= 0)
	{
		printf("length is zero or negative!\n");
		return -1;
	}
	T ret = array[0];
	for (int i = 0; i < length; i++)
		if (ret < array[i])
			ret = array[i];
	return ret;
}

template<typename T, typename T2>
T* CMax(T *array1, T *array2, T2 length)
{
	if (array1 == NULL || array2 == NULL)
	{
		printf("Input is empty!\n");
		return -1;
	}
	if (length <= 0)
	{
		printf("length is zero or negative!\n");
		return -1;
	}
	T *ret = new T[length];
	for (int i = 0; i < length; i++)
	{
		if (array1[i] > array2[i])
			ret[i] = array1[i];
		else
			ret[i] = array2[i];
	}
	return ret;
}

template<typename T, typename T2>
void CMinusAS(T *array, T2 length, T scalar, T *result)
{
	if (array == NULL)
	{
		printf("Input is empty!\n");
		exit(-1);
	}
	if (length <= 0)
	{
		printf("length is zero or negative!\n");
		exit(-1);
	}
	if (result == NULL)
		result = (T *)malloc(length*sizeof(T));
	for (int i = 0; i < length; i++)
	{
		result[i] = array[i] - scalar;
	}
}

template<typename T, typename T2>
void CPLUSAS(T *array, T2 length, T scalar, T *result)
{
	if (array == NULL)
	{
		printf("Input is empty!\n");
		exit(-1);
	}
	if (length <= 0)
	{
		printf("length is zero or negative!\n");
		exit(-1);
	}
	if (result == NULL)
		result = (T *)malloc(length*sizeof(T));
	for (int i = 0; i < length; i++)
	{
		result[i] = array[i] + scalar;
	}
}

template<typename T, typename T2>
void CProductAS(T *array, T2 length, T scalar, T *result)
{
	if (array == NULL)
	{
		printf("Input is empty!\n");
		exit(-1);
	}
	if (length <= 0)
	{
		printf("length is zero or negative!\n");
		exit(-1);
	}
	if (result == NULL)
		result = (T *)malloc(length*sizeof(T));
	for (int i = 0; i < length; i++)
	{
		result[i] = array[i] * scalar;
	}
}

template<typename T, typename T2>
void CDivideAS(T *array, T2 length, T scalar, T *result)
{
	if (array == NULL)
	{
		printf("Input is empty!\n");
		exit(-1);
	}
	if (length <= 0)
	{
		printf("length is zero or negative!\n");
		exit(-1);
	}
	if (scalar == 0)
	{
		printf("Denominator is zero!\n");
		exit(-1);
	}
	if (result == NULL)
		result = (T *)malloc(length*sizeof(T));
	for (int i = 0; i < length; i++)
	{
		result[i] = array[i] / scalar;
	}
}
template<typename T, typename T2>
void CFactorialDotProduct(T *array, T2 length, T *result)
{
	if (array == NULL)
	{
		printf("Input is empty!\n");
		exit(-1);
	}
	if (length <= 0)
	{
		printf("length is zero or negative!\n");
		exit(-1);
	}
	if (result == NULL)
		result = (T *)malloc(length*sizeof(T));
	T ret, swap;
	ret = 1;
	for (uint64_t i = 0; i < length; i ++)
	{
		swap = array[i];
		result[i] = ret;
		ret = ret * swap;
	}
}

template<typename T, typename L>
void CMemoryCopy(T *det, T *src, L length)
{
	if (src == NULL) {printf("Memory allocation failed!\n");exit(-1);}
	if (det == NULL) det = (T *)malloc(length*sizeof(T));
	if (det == NULL) {printf("Memory allocation failed!\n");exit(-1);}
	for (L i = 0; i < length; i++)
		det[i] = src[i];
}

template<typename T, typename L>
T *CMemoryInit(T val, L length)
{
	T *det = (T*) malloc(length*sizeof(T));
	if (det == NULL) {printf("Memory allocation failed!\n");exit(-1);}
	if (length < 0) {printf("Memory size can be negative!\n");exit(-1);}
	for (L i = 0; i < length; i++)
		det[i] = val;
	return det;
}

template<typename T>
void CMemoryRelease(T **det)
{
	if (det != NULL)
	{
		free(*det);
		*det = NULL;
	}
}

template<typename T1, typename T2, typename L>
T2* CDataCast(T1 *src, L length)
{
	if (src == NULL) {printf("Input is empty!\n");exit(-1);}
	T2 *ret = (T2 *)malloc(length*sizeof(T2));
	if (ret == NULL) {printf("Memory allocation failed!\n");exit(-1);}
	for (L i = 0; i < length; i++)
			ret[i] = (T2)src[i];
		return ret;
}

#endif /* SRC_CUTILITIES_HPP_ */

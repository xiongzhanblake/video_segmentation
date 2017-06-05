/*
 * CreateAdjacentMatrix.hpp
 *
 *  Created on: Feb 7, 2017
 *  Author: xiong
 */

#ifndef CREATADJACENTMATRIX_HPP_
#define CREATADJACENTMATRIX_HPP_

#define MINGRAPHWEIGHT 0.01

//#include "ComputeIndicesAndWeights.hpp"
#include "../../Common/VectorOperation.hpp"
#include "../../Common/ClassUtilities.hpp"
#include "../../Common/MexUtilities.hpp"
#include <Eigen/Sparse>

//typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SPMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> Triplet;

// #define CREATADJACENTMATRIX_HPP_DEBUG
// #define CREATADJACENTMATRIX_HPP_DETAILS

template<typename T, typename L>
void ReadOffSparseMatrix(int nlhs, mxArray *plhs[],Eigen::SparseMatrix<T, Eigen::ColMajor>* SPMat)
{
	if (nlhs != 3) {std::cout<<"output number is not 0!\n", exit(-1);}
	L Rows = SPMat->nonZeros();
//	L Cols = SPMat->outerSize();
	if (plhs[0] != NULL) mxDestroyArray(plhs[0]);
	if (plhs[1] != NULL) mxDestroyArray(plhs[1]);
	if (plhs[2] != NULL) mxDestroyArray(plhs[2]);
	plhs[0] = mxCreateNumericMatrix(Rows,1,mxDOUBLE_CLASS,mxREAL);
	plhs[1] = mxCreateNumericMatrix(Rows,1,mxDOUBLE_CLASS,mxREAL);
	plhs[2] = mxCreateNumericMatrix(Rows,1,mxDOUBLE_CLASS,mxREAL);
	double *Ic = (double *)mxGetData(plhs[0]);
	double *Ir = (double *)mxGetData(plhs[1]);
	double *Pv = (double *)mxGetData(plhs[2]);
	L index = 0;
	for (int k=0; k<SPMat->outerSize(); ++k)
		for (Eigen::SparseMatrix<double,Eigen::ColMajor>::InnerIterator it(*SPMat,k); it; ++it)
		{
			Pv[index] = it.value();
			Ir[index] = it.row();   // row index starts at 0;
			Ic[index] = it.col();   // col index starts at 0;
//			it.index(); // inner index, here it is equal to it.row()
			index++;
		}
}

template<typename T, typename L>
Eigen::SparseMatrix<T, Eigen::ColMajor>* ConstructSparseMatrix(T *Ir, T *Ic, T *Vp, L length, L Rows, L Cols)
{
	Eigen::SparseMatrix<T, Eigen::ColMajor> *mat = new Eigen::SparseMatrix<T, Eigen::ColMajor>(Rows,Cols);
	std::vector<Triplet> elements;
	for (L i = 0; i < length; i++)
		elements.push_back(Triplet(Ir[i],Ic[i],Vp[i]));
	mat->setFromTriplets(elements.begin(),elements.end());
	elements.clear();
	return mat;
}

template<typename T>
void checkBoundValidation(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
	if (nlhs != 0) {std::cout<<"output number is not 0!\n", exit(-1);}
	if (nrhs != 6) {std::cout<<"input number is not 6!\n", exit(-1);}
	T *centerIndices  = (T*)mxGetData(prhs[0]);
	T *leftIndices  = (T*)mxGetData(prhs[1]);
	T *rightIndices = (T*)mxGetData(prhs[2]);
	T leftBound  = ((T*)mxGetData(prhs[3]))[0];
	T rightBound = ((T*)mxGetData(prhs[4]))[0];
	T flag = ((T*)mxGetData(prhs[5]))[0];
	if (centerIndices == NULL || leftIndices == NULL || rightIndices == NULL) {std::cout<<"Empty inputs exist!"; exit(-1);}
	uint64_t length = mxGetM(prhs[0]);

	for (uint64_t i = 0; i < length; i++)
	{
		if (flag == 0)
		{
			leftIndices[i]  = (leftIndices[i]  - centerIndices[i]) < leftBound  ? 0 : leftIndices[i];
			rightIndices[i] = (rightIndices[i] - centerIndices[i]) > rightBound ? 0 : rightIndices[i];
		}
		else
		{
			leftIndices[i]  = (leftIndices[i]  - centerIndices[i]) < leftBound  ? (centerIndices[i] + rightBound) : leftIndices[i];
			rightIndices[i] = (rightIndices[i] - centerIndices[i]) > rightBound ? centerIndices[i] : rightIndices[i];
		}
	}
}

void CreateAdjacentMatrix( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
//-------------------------------------------------------------------------------------------------
// Check validation:
//-------------------------------------------------------------------------------------------------
	if (nrhs != 4)
		mexErrMsgTxt(Error_INPUT_NUMBER);
	if (nlhs != 3)
		mexErrMsgTxt(Error_OUTPUT_NUMBER);

	#ifdef SPLAT_ADJACENT_HPP_DETAILS
	SplitLine; std::cout<<"Function CreateAdjacentMatrix:";
    EndLine;
	#endif
    // Floor, weights, BilateralVals GridSize are all double type;
    const mxArray *DataIndices = prhs[0];
	const mxArray *DataWeights = prhs[1];
	const mxArray *GridSize = prhs[2];
	const mxArray *FactorWeight = prhs[3];
	mxArray *CenterModulo = NULL;
	mxArray *LeftIndices = NULL;
	mxArray *RightIndices = NULL;
	mxArray *DoubleScalar = NULL;
	mxArray *LeftBound = NULL;
	mxArray *RightBound = NULL;
	mxArray *wLeft = NULL;
	mxArray *wRight = NULL;

    mwSize nDims = mxGetN(FactorWeight) > mxGetM(FactorWeight) ? mxGetN(FactorWeight) : mxGetM(FactorWeight);
    mwSize nLength = mxGetN(DataIndices) > mxGetM(DataIndices) ? mxGetN(DataIndices) : mxGetM(DataIndices);
    double minGraphWeight = MINGRAPHWEIGHT;
    double *factorWeight = (double*)mxGetData(FactorWeight);
    double *gridsize = NULL;
    double offset = 0;
    double maxid = 0;
    double *doublescalar = NULL;
//-------------------------------------------------------------------------------------------------
// Preparation procedure:
//-------------------------------------------------------------------------------------------------
    int inNum = 5;
	int outNum = 2;
	const mxArray **inArrays = (const mxArray**)mxCalloc(inNum,sizeof(const mxArray*));
	for (int i = 0; i < inNum; i ++)
		inArrays[i] = NULL;
	mxArray **outArrays = (mxArray**)mxCalloc(outNum,sizeof(mxArray*));
	for (int i = 0; i < outNum; i ++)
		outArrays[i] = NULL;
	gridsize = (double *)mxCalloc((nDims+1),sizeof(double));
	gridsize[0] = 1;
	for (mwSize i = 1; i < nDims + 1; i ++)
		gridsize[i] = gridsize[i-1]*((double *)mxGetData(GridSize))[i-1];
	DoubleScalar = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
	doublescalar = (double *)mxGetData(DoubleScalar);
	LeftBound  = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
	RightBound = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);

	mxArray *sp_i = NULL;
	mxArray *sp_j = NULL;
	mxArray *sp_v = NULL;
	mxArray *tmp1, *tmp2;

	Eigen::SparseMatrix<double, Eigen::ColMajor> *B = NULL;
//-------------------------------------------------------------------------------------------------
// Computing procedure:
//-------------------------------------------------------------------------------------------------
	for (mwSize i = 0; i < nDims; i++)
	{
		Eigen::SparseMatrix<double, Eigen::ColMajor> *BD = NULL;

		if (factorWeight[i] == 0) continue;
		offset = gridsize[i];

		*doublescalar = offset;
		inArrays[0] = DataIndices; inArrays[1] = DoubleScalar; outArrays[0] = LeftIndices;
		mexSubstract<double>(1,outArrays,2,inArrays);
		LeftIndices = outArrays[0]; outArrays[0] = NULL;

//		inArrays[0] = LeftIndices; mexDisplayArray<double>(0,NULL,1,inArrays);

		inArrays[0] = DataIndices; inArrays[1] = DoubleScalar; outArrays[0] = RightIndices;
		mexAddition<double>(1,outArrays,2,inArrays);
		RightIndices = outArrays[0]; outArrays[0] = NULL;

//		inArrays[0] = DataIndices; mexDisplayArray<double>(0,NULL,1,inArrays);
//		inArrays[0] = RightIndices; mexDisplayArray<double>(0,NULL,1,inArrays);

		maxid = gridsize[i+1];

		// project onto the front dimension plane;
		// compute centerModulo = floor((centerModulo - 1) / maxid) * maxid;
		*doublescalar = 1;
		inArrays[0] = DataIndices; inArrays[1] = DoubleScalar; outArrays[0] = CenterModulo;
		mexSubstract<double>(1,outArrays,2,inArrays);
		CenterModulo = outArrays[0]; outArrays[0] = NULL;

//		inArrays[0] = CenterModulo; mexDisplayArray<double>(0,NULL,1,inArrays);

		*doublescalar = maxid;
		inArrays[0] = CenterModulo; inArrays[1] = DoubleScalar; outArrays[0] = CenterModulo;
		mexDotDivision<double>(1,outArrays,2,inArrays);
		CenterModulo = outArrays[0]; outArrays[0] = NULL;

//		inArrays[0] = CenterModulo; mexDisplayArray<double>(0,NULL,1,inArrays);

		inArrays[0] = CenterModulo; outArrays[0] = CenterModulo;
		mexFloor<double>(1,outArrays,1,inArrays);
		CenterModulo = outArrays[0]; outArrays[0] = NULL;

//		inArrays[0] = CenterModulo; mexDisplayArray<double>(0,NULL,1,inArrays);

		inArrays[0] = CenterModulo; inArrays[1] = DoubleScalar; outArrays[0] = CenterModulo;
		mexDotProduction<double>(1,outArrays,2,inArrays);
		CenterModulo = outArrays[0]; outArrays[0] = NULL;
		
//		inArrays[0] = CenterModulo; mexDisplayArray<double>(0,NULL,1,inArrays);

		// check if out of bounds
		*doublescalar = 0;
		if (i == nDims) *doublescalar = 1;
		((double*)mxGetData(LeftBound))[0] = 1; ((double*)mxGetData(RightBound))[0] = maxid;
		inArrays[0] = CenterModulo; inArrays[1] = LeftIndices, inArrays[2] = RightIndices;
		inArrays[3] = LeftBound; inArrays[4] = RightBound; inArrays[5] = DoubleScalar;
		checkBoundValidation<double>(0,NULL,6,(mxArray **)inArrays);

		// convert the indices into the occupied Vertex space
		// weight for an edge is the sum of the vertex weights
		*doublescalar = 1;
		inArrays[0] = LeftIndices, inArrays[1] = DataIndices; inArrays[2] = DoubleScalar; outArrays[0] = NULL; outArrays[1] = NULL; outArrays[2] = NULL;
		mexIntersectVector<double,double>(3,outArrays,3,inArrays);
		mxDestroyArray(outArrays[2]); mxDestroyArray(LeftIndices); mxDestroyArray(CenterModulo);
		CenterModulo = outArrays[0]; LeftIndices = outArrays[1];
		outArrays[0] = outArrays[1] = outArrays[2] = NULL;

//		inArrays[0] = LeftIndices; printf("LeftIndices:\n"); mexDisplayArray<double>(0,NULL,1,inArrays);
//		inArrays[0] = CenterModulo; printf("CenterModulo:\n");mexDisplayArray<double>(0,NULL,1,inArrays);

		inArrays[0] = sp_i; inArrays[1] = CenterModulo; outArrays[0] = sp_i;
		mexVertcat<double>(1,outArrays,2,inArrays);
		mxDestroyArray(sp_i); sp_i = outArrays[0]; outArrays[0] = NULL;

//		inArrays[0] = sp_i; printf("sp_i:\n"); mexDisplayArray<double>(0,NULL,1,inArrays);

		inArrays[0] = sp_j; inArrays[1] = LeftIndices; outArrays[0] = sp_j;
		mexVertcat<double>(1,outArrays,2,inArrays);
		mxDestroyArray(sp_j); sp_j = outArrays[0]; outArrays[0] = NULL;

//		inArrays[0] = sp_j; printf("sp_j:\n"); mexDisplayArray<double>(0,NULL,1,inArrays);

		inArrays[0] = DataWeights; inArrays[1] = CenterModulo; outArrays[0] = NULL;
		mexExtractRow<double,uint64_t>(1,outArrays,2,inArrays);
		tmp1 = outArrays[0]; outArrays[0] = NULL;

		inArrays[0] = DataWeights; inArrays[1] = LeftIndices; outArrays[0] = NULL;
		mexExtractRow<double,uint64_t>(1,outArrays,2,inArrays);
		tmp2 = outArrays[0]; outArrays[0] = NULL;
		
		inArrays[0] = tmp2; inArrays[1] = tmp1; outArrays[0] = wLeft;
		mexAddition<double>(1,outArrays,2,inArrays);
		wLeft = outArrays[0]; outArrays[0] = NULL;
		mxDestroyArray(tmp1); tmp1 = NULL;
		mxDestroyArray(tmp2); tmp2 = NULL;

		*doublescalar = factorWeight[i];
		inArrays[0] = wLeft; inArrays[1] = DoubleScalar; outArrays[0] = wLeft;
		mexDotProduction<double>(1,outArrays,2,inArrays);
		wLeft = outArrays[0]; outArrays[0] = NULL;

		inArrays[0] = sp_v; inArrays[1] = wLeft; outArrays[0] = NULL;
		mexVertcat<double>(1,outArrays,2,inArrays);
		mxDestroyArray(sp_v); sp_v = outArrays[0]; outArrays[0] = NULL;
		mxDestroyArray(wLeft); wLeft = NULL;
		mxDestroyArray(LeftIndices); LeftIndices = NULL;
		mxDestroyArray(CenterModulo); CenterModulo = NULL;
//		inArrays[0] = sp_v; printf("sp_v:\n"); mexDisplayArray<double>(0,NULL,1,inArrays);

		*doublescalar = minGraphWeight;
		inArrays[0] = sp_v; inArrays[1] = DoubleScalar; outArrays[0] = sp_v;
		mexDotMax<double>(1,outArrays,2,inArrays);
		sp_v = outArrays[0]; outArrays[0] = NULL;

//		BD1 = new SparseMatrix<double, double>((double *)mxGetData(sp_i), (double *)mxGetData(sp_j), (double *)mxGetData(sp_v), (double)mxGetM(sp_v), nLength, nLength);
//		mxDestroyArray(sp_i); mxDestroyArray(sp_j); mxDestroyArray(sp_v);
//		sp_i = sp_j = sp_v = NULL;

		// compute right indices
		*doublescalar = 1;
		inArrays[0] = RightIndices, inArrays[1] = DataIndices; inArrays[2] = DoubleScalar; outArrays[0] = NULL; outArrays[1] = NULL; outArrays[2] = NULL;
		mexIntersectVector<double,double>(3,outArrays,3,inArrays);
		mxDestroyArray(outArrays[2]); mxDestroyArray(RightIndices); mxDestroyArray(CenterModulo);
		CenterModulo = outArrays[0]; RightIndices = outArrays[1];
		outArrays[0] = outArrays[1] = outArrays[2] = NULL;

//		inArrays[0] = RightIndices; printf("RightIndices:\n"); mexDisplayArray<double>(0,NULL,1,inArrays);
//		inArrays[0] = CenterModulo; printf("CenterModulo:\n");mexDisplayArray<double>(0,NULL,1,inArrays);

		inArrays[0] = sp_i; inArrays[1] = CenterModulo; outArrays[0] = NULL;
		mexVertcat<double>(1,outArrays,2,inArrays);
		mxDestroyArray(sp_i); sp_i = outArrays[0]; outArrays[0] = NULL;

//		inArrays[0] = sp_i; printf("sp_i:\n"); mexDisplayArray<double>(0,NULL,1,inArrays);

		inArrays[0] = sp_j; inArrays[1] = RightIndices; outArrays[0] = NULL;
		mexVertcat<double>(1,outArrays,2,inArrays);
		mxDestroyArray(sp_j); sp_j = outArrays[0]; outArrays[0] = NULL;

//		inArrays[0] = sp_j; printf("sp_j:\n"); mexDisplayArray<double>(0,NULL,1,inArrays);

		inArrays[0] = DataWeights; inArrays[1] = CenterModulo; outArrays[0] = NULL;
		mexExtractRow<double,double>(1,outArrays,2,inArrays);
		tmp1 = outArrays[0]; outArrays[0] = NULL;

		inArrays[0] = DataWeights; inArrays[1] = RightIndices; outArrays[0] = NULL;
		mexExtractRow<double,double>(1,outArrays,2,inArrays);
		tmp2 = outArrays[0]; outArrays[0] = NULL;
		
		inArrays[0] = tmp2; inArrays[1] = tmp1; outArrays[0] = wRight;
		mexAddition<double>(1,outArrays,2,inArrays);
		wRight = outArrays[0]; outArrays[0] = NULL;
		mxDestroyArray(tmp1); tmp1 = NULL;
		mxDestroyArray(tmp2); tmp2 = NULL;

		*doublescalar = factorWeight[i];
		inArrays[0] = wRight; inArrays[1] = DoubleScalar; outArrays[0] = wRight;
		mexDotProduction<double>(1,outArrays,2,inArrays);
		wRight = outArrays[0]; outArrays[0] = NULL;

		inArrays[0] = sp_v; inArrays[1] = wRight; outArrays[0] = NULL;
		mexVertcat<double>(1,outArrays,2,inArrays);
		mxDestroyArray(sp_v); sp_v = outArrays[0]; outArrays[0] = NULL;
		mxDestroyArray(wRight); wRight = NULL;
		mxDestroyArray(RightIndices); RightIndices = NULL;
		mxDestroyArray(CenterModulo); CenterModulo = NULL;

//		inArrays[0] = sp_v; printf("sp_v:\n"); mexDisplayArray<double>(0,NULL,1,inArrays);

		*doublescalar = minGraphWeight;
		inArrays[0] = sp_v; inArrays[1] = DoubleScalar; outArrays[0] = sp_v;
		mexDotMax<double>(1,outArrays,2,inArrays);
		sp_v = outArrays[0]; outArrays[0] = NULL;

//		inArrays[0] = sp_v; printf("sp_v:\n"); mexDisplayArray<double>(0,NULL,1,inArrays);

		BD = ConstructSparseMatrix<double,mwSize>((double *)mxGetData(sp_i), (double *)mxGetData(sp_j), (double *)mxGetData(sp_v), (mwSize)mxGetM(sp_v), nLength, nLength );
		mxDestroyArray(sp_i); mxDestroyArray(sp_j); mxDestroyArray(sp_v);
		sp_i = sp_j = sp_v = NULL;

		if (i == 0)
		{
			B = BD;
			BD = NULL;
		}
		else
		{
			*B = *B + *BD;
			delete BD;
			BD = NULL;
		}
	}

	for (int i = 0; i < inNum; i ++)
			inArrays[i] = NULL;
	mxFree(inArrays); inArrays = NULL;
	for (int i = 0; i < outNum; i ++)
			outArrays[i] = NULL;
	mxFree(outArrays); outArrays = NULL;

	mxDestroyArray(LeftBound); LeftBound = NULL;
	mxDestroyArray(RightBound); RightBound = NULL;
	mxDestroyArray(DoubleScalar); DoubleScalar = NULL;
	mxFree(gridsize); gridsize = NULL;

	if (plhs[0] != NULL) mxDestroyArray(plhs[0]); plhs[0] = NULL;
	if (plhs[1] != NULL) mxDestroyArray(plhs[1]); plhs[1] = NULL;
	if (plhs[2] != NULL) mxDestroyArray(plhs[2]); plhs[2] = NULL;
	ReadOffSparseMatrix<double,uint64_t>(3,plhs,B);
	delete B; B = NULL;
//	delete inArrays;
//	delete outArrays;
}

#endif


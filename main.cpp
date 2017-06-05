/*
 * main.cpp
 *
 *  Created on: May 16, 2017
 *      Author: xiong
 */

#include "mat.h"
#include "mex.h"
#include "matrix.h"

#include "Common/VectorOperation.hpp"
#include "Core/Graph_Cut/CreateAdjacentMatrix.hpp"
#include "Interface/Interface.hpp"
//#include "GCoptimization.h"
//#include "ClassUtilities.hpp"

int main()
{
	MATFile *DataMat;
	MATFile *ResultMat;
	const char *file_data = "data.mat";
	const char *file_para = "para.mat";
	int nrhs = 10;
	int nlhs = 5;
	const mxArray **prhs = (const mxArray**)mxCalloc(nrhs,sizeof(const mxArray*));
	mxArray **plhs = (mxArray**)mxCalloc(nlhs,sizeof(mxArray*));
	for (int i = 0; i < nlhs; i++)
		plhs[i] = NULL;
	for (int i = 0; i < nrhs; i++)
		prhs[i] = NULL;

	printf("-------------------------------------------------------------------------\n");
	printf("Function Main...\n");
	printf("-------------------------------------------------------------------------\n");

	DataMat = matOpen(file_data,"r");
	mxArray * MaskVisualCues = matGetVariable(DataMat,"MaskVisualCues");
	mxArray * MaskMotionCues = matGetVariable(DataMat,"MaskMotionCues");
	mxArray * DataVisualCues = matGetVariable(DataMat,"DataVisualCues");
	mxArray * DataMotionCues = matGetVariable(DataMat,"DataMotionCues");
	mxArray * Mask = matGetVariable(DataMat,"Mask");
	mxArray * MotionGrid = matGetVariable(DataMat,"MotionGrid");
	mxArray * VisualGrid = matGetVariable(DataMat,"VisualGrid");
	mxArray * nLabels = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS, mxREAL);
	((double *)mxGetData(nLabels))[0] = 2;
	matClose(DataMat);

	prhs[0] = MaskVisualCues;
	prhs[1] = MaskMotionCues;
	prhs[2] = DataVisualCues;
	prhs[3] = DataMotionCues;
	prhs[4] = Mask;
	prhs[5] = MotionGrid;
	prhs[6] = VisualGrid;
	prhs[7] = nLabels;
	plhs[0] = plhs[1] = plhs[2] = plhs[3] = NULL;

	LT_SP_Interface<double>( 4, plhs, 8, prhs);
	mxDestroyArray(MaskVisualCues); MaskVisualCues = NULL;
	mxDestroyArray(MaskMotionCues); MaskMotionCues = NULL;
	mxDestroyArray(DataVisualCues); DataVisualCues = NULL;
	mxDestroyArray(DataMotionCues); DataMotionCues = NULL;
	mxDestroyArray(Mask); Mask = NULL;
	mxDestroyArray(MotionGrid); MotionGrid = NULL;
	mxDestroyArray(VisualGrid); VisualGrid = NULL;
	mxDestroyArray(nLabels); nLabels = NULL;

	mxArray *DataIndices = plhs[0];
	mxArray *DataWeights = plhs[1];
	mxArray *MaskWeights = plhs[2];
	mxArray *GridSize = plhs[3];

	DataMat = matOpen(file_para,"r");
	mxArray *VisualWeights = matGetVariable(DataMat,"VisualWeights");
	mxArray *MotionWeights = matGetVariable(DataMat,"MotionWeights");
	matClose(DataMat);
	prhs[0] = DataIndices; prhs[1] = DataWeights; prhs[2] = GridSize;
	prhs[3] = VisualWeights; prhs[4] = MotionWeights;
	plhs[0] = plhs[1] = plhs[2] = plhs[3] = NULL;
	SP_AM_Interface<double>( 3, plhs, 5, prhs);

	mxArray *Sp_i = plhs[0];
	mxArray *Sp_j = plhs[1];
	mxArray *Sp_v = plhs[2];

	mxDestroyArray(VisualWeights); VisualWeights = NULL;
	mxDestroyArray(MotionWeights); MotionWeights = NULL;

	mxArray *nLength = mxCreateNumericMatrix(1,1,mxINT64_CLASS,mxREAL);
	nLabels = mxCreateNumericMatrix(1,1,mxINT64_CLASS,mxREAL);
	((int64_t *)mxGetData(nLength))[0] = mxGetM(DataIndices);
	((int64_t *)mxGetData(nLabels))[0] = mxGetN(MaskWeights);
//
//	printf("reading data...\n");
	DataMat = matOpen(file_para,"r");
	mxArray *unaryWeight = matGetVariable(DataMat,"unaryWeight");
	mxArray *pairwiseWeight = matGetVariable(DataMat,"pairwiseWeight");
	matClose(DataMat);
	prhs[0] = MaskWeights;
	prhs[1] = Sp_i;
	prhs[2] = Sp_j;
	prhs[3] = Sp_v;
	prhs[4] = nLength;
	prhs[5] = nLabels;
	prhs[6] = unaryWeight;
	prhs[7] = pairwiseWeight;
	plhs[0] = NULL;

	// An interface from Create_Adjacent_Matrix to GC_optimization;
	AM_GC_Interface<double>(1, plhs, 8, prhs);
	mxArray *Labels = plhs[0];
	plhs[0] = NULL;

	mxDestroyArray(Sp_i); Sp_i = NULL;
	mxDestroyArray(Sp_j); Sp_j = NULL;
	mxDestroyArray(Sp_v); Sp_v = NULL;
	mxDestroyArray(MaskWeights); MaskWeights = NULL;
	mxDestroyArray(DataWeights); DataWeights = NULL;
	mxDestroyArray(unaryWeight); unaryWeight = NULL;
	mxDestroyArray(pairwiseWeight); pairwiseWeight = NULL;
	mxDestroyArray(nLabels); nLabels = NULL;
	mxDestroyArray(nLength); nLength = NULL;

	DataMat = matOpen(file_para,"r");
	mxArray *VideoSize = matGetVariable(DataMat,"VideoSize");
	matClose(DataMat);

	DataMat = matOpen(file_data,"r");
	DataVisualCues = matGetVariable(DataMat,"DataVisualCues");
	DataMotionCues = matGetVariable(DataMat,"DataMotionCues");
	matClose(DataMat);
	prhs[0] = DataIndices;
	prhs[1] = Labels;
	prhs[2] = DataVisualCues;
	prhs[3] = DataMotionCues;
	prhs[4] = GridSize;
	prhs[5] = VideoSize;
	plhs[0] = NULL;

	GC_SL_Interface<double,uint64_t>(1,plhs,6,prhs);
	mxArray *segmentation = plhs[0]; plhs[0] = NULL;

	mxDestroyArray(DataIndices); DataIndices = NULL;
	mxDestroyArray(Labels); Labels = NULL;
	mxDestroyArray(DataVisualCues); DataVisualCues = NULL;
	mxDestroyArray(DataMotionCues); DataMotionCues = NULL;
	mxDestroyArray(GridSize); GridSize = NULL;
	mxDestroyArray(VideoSize); VideoSize = NULL;

	printf("-------------------------------------------------------------------------\n");
	printf("Function Main...\n");
	printf("-------------------------------------------------------------------------\n");
	printf("-------------------------------------------------------------------------\n");
	printf("Save Result...\n");
	const char *file_Segmentation = "segmentation.mat";
	ResultMat = matOpen(file_Segmentation,"w");
	matPutVariable(ResultMat, "segmentation", segmentation);
	matClose(ResultMat);
	printf("-------------------------------------------------------------------------\n");
	printf("Releasing Memory...\n");
	for (int i = 0; i < nrhs; i++) 	prhs[i] = NULL;
	mxFree(prhs); prhs = NULL;
	for (int i = 0; i < nlhs; i++) plhs[i] = NULL;
	mxFree(plhs); plhs = NULL;
	printf("-------------------------------------------------------------------------\n");
	printf("Exit with success!\n");
	printf("-------------------------------------------------------------------------\n");
	return 0;
}


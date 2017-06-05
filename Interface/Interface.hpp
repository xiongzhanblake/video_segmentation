/*
 * GCInit.cpp
 *
 *  Created on: May 28, 2017
 *      Author: xiong
 */
#include "mat.h"
#include "mex.h"
#include "matrix.h"
#include <Eigen/Dense>
#include "../Core/Splat_Slice/MexSplat_Adjacent.hpp"
#include "../Core/Splat_Slice/MexSlice_Adjacent.hpp"
#include "../Core/Splat_Slice/MexIntersection.hpp"
#include "../Core/Graph_Cut/GCoptimization.h"

template <typename ctype> struct ctype2mx { };
template <> struct ctype2mx<bool>            { static const mxClassID classid = mxLOGICAL_CLASS;   static const char* classname() { return "logical"; } };
template <> struct ctype2mx<float>           { static const mxClassID classid = mxSINGLE_CLASS;    static const char* classname() { return "single"; }  };
template <> struct ctype2mx<double>          { static const mxClassID classid = mxDOUBLE_CLASS;    static const char* classname() { return "double"; }  };
template <> struct ctype2mx<char>            { static const mxClassID classid = mxINT8_CLASS;      static const char* classname() { return "int8"; }    };
template <> struct ctype2mx<unsigned char>   { static const mxClassID classid = mxUINT8_CLASS;     static const char* classname() { return "uint8"; }   };
template <> struct ctype2mx<short>           { static const mxClassID classid = mxINT16_CLASS;     static const char* classname() { return "int16"; }   };
template <> struct ctype2mx<unsigned short>  { static const mxClassID classid = mxUINT16_CLASS;    static const char* classname() { return "uint16"; }  };
template <> struct ctype2mx<int>             { static const mxClassID classid = mxINT32_CLASS;     static const char* classname() { return "int32"; }   };
template <> struct ctype2mx<unsigned int>    { static const mxClassID classid = mxUINT32_CLASS;    static const char* classname() { return "uint32"; }  };
template <> struct ctype2mx<long long>       { static const mxClassID classid = mxINT64_CLASS;     static const char* classname() { return "int64"; }   };
template <> struct ctype2mx<unsigned long long> { static const mxClassID classid = mxUINT64_CLASS; static const char* classname() { return "uint64"; }  };

template<typename T>
void AM_GC_Interface(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	typedef long long LabelID;
	typedef long long SiteID;
	typedef float EnergyType;
	typedef float EnergyTermType;
	mxClassID cLabelClassID      = ctype2mx<LabelID>::classid;
//	mxClassID cSiteClassID       = ctype2mx<SiteID>::classid;
	mxClassID cEnergyTermClassID = ctype2mx<EnergyTermType>::classid;
//	mxClassID cEnergyClassID     = ctype2mx<EnergyType>::classid;
	typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> MatrixXTC;
	typedef Eigen::Matrix<EnergyTermType, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> MatrixXEC;
	typedef Eigen::Map<MatrixXTC> MapMatrixXTC;
	typedef Eigen::Map<MatrixXEC> MapMatrixXEC;

	printf("-------------------------------------------------------------------------\n");
	printf("Function Graph-Cut optimization...\n");

	if (nlhs != 1) {printf("Input is not 1!");exit(-1);}
	if (nrhs != 8) {printf("Input is not 8!");exit(-1);}

	const mxArray *MaskWeight = prhs[0];
	const mxArray *Ir = prhs[1];
	const mxArray *Jc = prhs[2];
	const mxArray *Pr = prhs[3];
	const mxArray *nLength = prhs[4];
	const mxArray *nLabels = prhs[5];
	const mxArray *UnaryWeigths = prhs[6];
	const mxArray *PairWiseWeights = prhs[7];

//	mxArray *TransposedMaskWeight = mxCreateNumericMatrix(mxGetN(MaskWeight),mxGetM(MaskWeight),mxGetClassID(MaskWeight),mxREAL);
	mxArray *TransposedMaskWeight = mxCreateNumericMatrix(mxGetN(MaskWeight),mxGetM(MaskWeight),cEnergyTermClassID,mxREAL);

	// Computing SPM as datacost in GC;
	MapMatrixXTC *map = new MapMatrixXTC((T *)mxGetData(MaskWeight),mxGetM(MaskWeight),mxGetN(MaskWeight));
	MapMatrixXEC *tmp = new MapMatrixXEC((EnergyTermType *)mxGetData(TransposedMaskWeight),mxGetM(TransposedMaskWeight),mxGetN(TransposedMaskWeight));
//	printf("%f\n",((T*)mxGetData(TransposedMaskWeight))[31129]);
	*tmp = (map->transpose() * ((T *)mxGetData(UnaryWeigths))[0]).template cast<EnergyTermType>();
//	printf("%f\n",((T*)mxGetData(TransposedMaskWeight))[31129]);
	delete map; map = NULL;
	delete tmp; tmp = NULL;

	// Computing SMC as smoothcost in GC;
	mxArray *SMC = mxCreateNumericMatrix(2,2,cEnergyTermClassID,mxREAL);
	tmp = new MapMatrixXEC((EnergyTermType *)mxGetData(SMC),mxGetM(SMC),mxGetN(SMC));
//	printf("%f\n",((T*)mxGetData(SMC))[1]);
	*tmp<<0,1,1,0;
//	printf("%f\n",((T*)mxGetData(SMC))[1]);
	*tmp = *tmp * (EnergyTermType)((T *)mxGetData(PairWiseWeights))[0];
//	std::cout<<*tmp<<std::endl;
	delete tmp; tmp = NULL;

	SiteID num_pixels = (SiteID)((SiteID *)mxGetData(nLength))[0];
	LabelID num_labels = (LabelID)((LabelID *)mxGetData(nLabels))[0];
	printf("-------------------------------------------------------------------------\n");
	printf("Site -> %lld, Labels -> %lld\n",num_pixels, num_labels);
	GCoptimizationGeneralGraph<EnergyTermType, EnergyType, SiteID, LabelID> *gc = new GCoptimizationGeneralGraph<EnergyTermType, EnergyType, SiteID, LabelID>(num_pixels,num_labels);

	/* set the nieghbors and weights according to sparse matrix */
//	printf("Read is done...\n");
//	printf("Setting data...\n");

	// Setting neighborhood system in GC;
	/* Get the starting positions of all four data arrays. */
	double *ir = (double *)mxGetData(Ir);
	double *jc = (double *)mxGetData(Jc);
	double *pr = (double *)mxGetData(Pr);

	for (mwSize total = 0; total < mxGetM(Ir); total ++)
		if (ir[total] >= jc[total])
			gc->setNeighbors((SiteID)ir[total], (SiteID)jc[total], (EnergyTermType)pr[total]);

	EnergyTermType *DataCost = (EnergyTermType *)mxGetData(TransposedMaskWeight);
	DataCostFn<EnergyTermType, SiteID, LabelID> *dc = new DataCostFnFromArray<EnergyTermType, SiteID, LabelID>(DataCost,num_labels);
	EnergyTermType *SmoothnessCost = (EnergyTermType *)mxGetData(SMC);
	SmoothCostFn<EnergyTermType, SiteID, LabelID> *sc = new SmoothCostFnFromArray<EnergyTermType, SiteID, LabelID>(SmoothnessCost,num_labels);
	/* set data term */
	gc->setDataCost(dc);
	/* set the smoothness term */
	gc->setSmoothCost(sc);
//	printf("Setting is done...\n");

	// Computing labels;
//	printf("Label expansion...\n");
	printf("-------------------------------------------------------------------------\n");
	printf("Before minimization: smoothcost -> %1.0f\tdatacost -> %1.0f\n", gc->giveSmoothEnergy(),gc->giveDataEnergy());
	printf("-------------------------------------------------------------------------\n");
	gc->expansion(-1);

	if (plhs[0] != NULL) mxDestroyArray(plhs[0]);
	plhs[0] = mxCreateNumericMatrix(gc->numSites(), 1, cLabelClassID, mxREAL);
	gc->whatLabel(0, mxGetNumberOfElements(plhs[0]),(LabelID*) mxGetData(plhs[0]));
	printf("After minimization: smoothcost -> %1.0f\tdatacost -> %1.0f\n", gc->giveSmoothEnergy(),gc->giveDataEnergy());
	printf("-------------------------------------------------------------------------\n");
//	printf("label expansion is done...\n");

//	printf("deleting Graph-Cut optimizer...\n");
	delete gc; gc = NULL;
	mxDestroyArray(TransposedMaskWeight); TransposedMaskWeight = NULL;
	mxDestroyArray(SMC); SMC = NULL;
//	printf("deletion is done...\n");
}
// Input list: 1. DataIndices; 2. DataWeights; 3. Gridsize; 4. FactorWeights;
template<typename T>
void SP_AM_Interface(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	printf("-------------------------------------------------------------------------\n");
	printf("Function Create Adjacent Matrix:\n");
	printf("-------------------------------------------------------------------------\n");
//	printf("reading data...\n");
//	printf("Creating adjacent matrix...\n");

	if (nlhs != 3) {printf("Input is not 3!");exit(-1);}
	if (nrhs != 5) {printf("Input is not 5!");exit(-1);}

	// Computing FactorWeights
	int inNum = 4; int outNum = 3;
	const mxArray **parameters = (const mxArray **)mxCalloc(inNum,sizeof(const mxArray *));
	mxArray **results = (mxArray **)mxCalloc(outNum,sizeof(mxArray *));
	parameters[0] = prhs[3]; parameters[1] = prhs[4]; results[0] = NULL;
	mexHorzcat<T>(1,results,2,parameters);
	mxArray *FactorWeights = results[0]; results[0] = NULL;

	// Computing adjacent matrix
	parameters[0] = prhs[0]; parameters[1] = prhs[1]; parameters[2] = prhs[2]; parameters[3] = FactorWeights;
	results[0] = results[1] = results[2] = NULL;
	CreateAdjacentMatrix( 3, plhs, 4, parameters);
	mxDestroyArray(FactorWeights); FactorWeights = NULL;
//	printf("Adjacent matrix is done...\n");
	for (int i = 0; i < inNum; i++) parameters[i] = NULL;
	mxFree(parameters); parameters = NULL;
	for (int i = 0; i < outNum; i++) results[i] = NULL;
	mxFree(results); results = NULL;
}
template<typename T>
void LT_SP_Interface(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	typedef Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> ArrayXT;
	typedef Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> ArrayXB;
	typedef Eigen::Map<ArrayXT> MapArrayXT;
//	typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MatrixXT;
//	typedef Eigen::Map<MatrixXT> MapMatrixXT;
//	typedef Eigen::Map<ArrayXB> MapArrayXB;

	if (nlhs != 4) {printf("Input is not 3!");exit(-1);}
	if (nrhs != 8) {printf("Input is not 4!");exit(-1);}

	const mxArray *MaskVisualCues = prhs[0];
	const mxArray *MaskMotionCues = prhs[1];
	const mxArray *DataVisualCues = prhs[2];
	const mxArray *DataMotionCues = prhs[3];
	const mxArray *Mask = prhs[4];
	const mxArray *MotionGridSize = prhs[5];
	const mxArray *VisualGridSize = prhs[6];
	mwSize nLabels = (mwSize)((T*)mxGetData(prhs[7]))[0];

	if (nLabels != 2) {printf("Current only support 2 labels!");exit(-1);}

	mxArray *MaskLabels = mxCreateNumericMatrix(mxGetM(MaskVisualCues),nLabels,mxGetClassID(prhs[0]),mxREAL);
	mxArray *DataLabels = mxCreateNumericMatrix(mxGetM(DataVisualCues),1,mxGetClassID(prhs[0]),mxREAL);
	mxArray *GridSize = NULL;
	mxArray *MotionCues = NULL;
	mxArray *VisualCues = NULL;
	mxArray *Cues = NULL;
	mxArray *Round = NULL;
	mxArray *WeightCeil = NULL;
	mxArray *WeightFloor = NULL;
	mxArray *DataInfor = NULL;
	mxArray *MaskInfor = NULL;

	mxArray *DoubleScalar = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);

	MapArrayXT *map = new MapArrayXT((T *)mxGetData(DataLabels),mxGetM(DataLabels),mxGetN(DataLabels));
	map->setOnes();
	delete map; map = NULL;

	map = new MapArrayXT((T *)mxGetData(MaskLabels),mxGetM(MaskLabels),mxGetN(MaskLabels));
	MapArrayXT *mask = new MapArrayXT((T *)mxGetData(Mask),mxGetNumberOfElements(Mask),1);
	ArrayXB *ret = new ArrayXB(mask->rows(),nLabels);
	ret->col(0) = (*mask!=0); ret->col(1) = (*mask==0);
	*map = ret->template cast<T>();
	delete ret; ret = NULL; delete mask; mask = NULL; delete map; map = NULL;

	map = new MapArrayXT((T *)mxGetData(MaskMotionCues),mxGetM(MaskMotionCues),mxGetN(MaskMotionCues));
	map->col(0).swap(map->col(2));
	delete map; map = NULL;

	map = new MapArrayXT((T *)mxGetData(DataMotionCues),mxGetM(DataMotionCues),mxGetN(DataMotionCues));
	map->col(0).swap(map->col(2));
	delete map; map = NULL;

	map = new MapArrayXT((T *)mxGetData(MotionGridSize),mxGetM(MotionGridSize),mxGetN(MotionGridSize));
	map->col(0).swap(map->col(2));
	delete map; map = NULL;

	int inNum = 10; int outNum = 5;
	const mxArray **parameters = (const mxArray **)mxCalloc(inNum,sizeof(const mxArray *));
	mxArray **results = (mxArray **)mxCalloc(outNum,sizeof(mxArray *));

	parameters[0] = VisualGridSize; parameters[1] = MotionGridSize; results[0] = GridSize;
	mexHorzcat<T>(1,results,2,parameters);
	GridSize = results[0]; results[0] = NULL;

	// Splat DataCues;
	parameters[0] = DataMotionCues; results[0] = MotionCues;
	mexRound<T>(1,results,1,parameters);
	MotionCues = results[0]; results[0] = NULL;

	parameters[0] = DataVisualCues; results[0] = VisualCues;
	mexFloor<T>(1,results,1,parameters);
	VisualCues = results[0]; results[0] = NULL;

	parameters[0] = VisualCues; parameters[1] = MotionCues; results[0] = Round;
	mexHorzcat<T>(1,results,2,parameters);
	Round = results[0]; results[0] = NULL;
	mxDestroyArray(MotionCues); mxDestroyArray(VisualCues);
	MotionCues = VisualCues = NULL;

	parameters[0] = DataVisualCues; parameters[1] = DataMotionCues; results[0] = Cues;
	mexHorzcat<T>(1,results,2,parameters);
	Cues = results[0]; results[0] = NULL;

	parameters[0] = Cues; parameters[1] = Round; results[0] = WeightCeil;
	mexSubstract<T>(1,results,2,parameters);
	WeightCeil = results[0]; results[0] = NULL;
	mxDestroyArray(Cues);
	Cues = NULL;

	((double *)mxGetData(DoubleScalar))[0] = 1;
	parameters[0] = DoubleScalar; parameters[1] = WeightCeil; results[0] = WeightFloor;
	mexSubstract<T>(1,results,2,parameters);
	WeightFloor = results[0]; results[0] = NULL;

	((double *)mxGetData(DoubleScalar))[0] = (double)mxGetN(DataVisualCues);
	parameters[0] = Round; parameters[1] = DataLabels; parameters[2] = GridSize; parameters[3] = WeightFloor; parameters[4] = WeightCeil;
	parameters[5] = DoubleScalar; results[0] = DataInfor;
	MexSplat_Adjacent(1,results,6,parameters);
	DataInfor = results[0]; results[0] = NULL;
	mxDestroyArray(Round); Round = NULL;
	mxDestroyArray(DataLabels); DataLabels = NULL;
	mxDestroyArray(WeightFloor); WeightFloor = NULL;
	mxDestroyArray(WeightCeil); WeightCeil = NULL;
	// Splat MaskCues;
	parameters[0] = MaskMotionCues; results[0] = MotionCues;
	mexRound<T>(1,results,1,parameters);
	MotionCues = results[0]; results[0] = NULL;

	parameters[0] = MaskVisualCues; results[0] = VisualCues;
	mexFloor<T>(1,results,1,parameters);
	VisualCues = results[0]; results[0] = NULL;

	parameters[0] = VisualCues; parameters[1] = MotionCues; results[0] = Round;
	mexHorzcat<T>(1,results,2,parameters);
	Round = results[0]; results[0] = NULL;
	mxDestroyArray(MotionCues); mxDestroyArray(VisualCues);
	MotionCues = VisualCues = NULL;

	parameters[0] = MaskVisualCues; parameters[1] = MaskMotionCues; results[0] = Cues;
	mexHorzcat<T>(1,results,2,parameters);
	Cues = results[0]; results[0] = NULL;

	parameters[0] = Cues; parameters[1] = Round; results[0] = WeightCeil;
	mexSubstract<T>(1,results,2,parameters);
	WeightCeil = results[0]; results[0] = NULL;
	mxDestroyArray(Cues);
	Cues = NULL;

	((double *)mxGetData(DoubleScalar))[0] = 1;
	parameters[0] = DoubleScalar; parameters[1] = WeightCeil; results[0] = WeightFloor;
	mexSubstract<T>(1,results,2,parameters);
	WeightFloor = results[0]; results[0] = NULL;

	((double *)mxGetData(DoubleScalar))[0] = (double)mxGetN(MaskVisualCues);
	parameters[0] = Round; parameters[1] = MaskLabels; parameters[2] = GridSize; parameters[3] = WeightFloor; parameters[4] = WeightCeil;
	parameters[5] = DoubleScalar; results[0] = MaskInfor;
	MexSplat_Adjacent(1,results,6,parameters);
	MaskInfor = results[0]; results[0] = NULL;
	mxDestroyArray(Round); Round = NULL;
	mxDestroyArray(MaskLabels); MaskLabels = NULL;
	mxDestroyArray(WeightFloor); WeightFloor = NULL;
	mxDestroyArray(WeightCeil); WeightCeil = NULL;

	// Post-process DataInfor and MaskInfor
	mxArray *cell = mxGetCell(DataInfor,0);
	map = new MapArrayXT((T *)mxGetData(cell),mxGetM(cell),mxGetN(cell));
	mxArray *DataIndices = mxCreateNumericMatrix(mxGetM(cell),1,mxGetClassID(cell),mxREAL);
	mxArray *DataWeights = mxCreateNumericMatrix(mxGetM(cell),1,mxGetClassID(cell),mxREAL);
	mxArray *MaskWeights = NULL;
	MapArrayXT *col1 = new MapArrayXT((T *)mxGetData(DataIndices),mxGetM(DataIndices),mxGetN(DataIndices));
	MapArrayXT *col2 = new MapArrayXT((T *)mxGetData(DataWeights),mxGetM(DataWeights),mxGetN(DataWeights));
	*col1 = map->col(0); *col2 = map->col(1);
	delete map; map = NULL;
	delete col1; col1 = NULL;
	delete col2; col2 = NULL;
	mxDestroyArray(DataInfor); DataInfor = NULL;
	cell = NULL;

	parameters[0] = MaskInfor; parameters[1] = DataIndices; results[0] = MaskWeights;
	MexIntersection(1,results,2,parameters);
	MaskWeights = results[0]; results[0] = NULL;
	mxDestroyArray(MaskInfor); MaskInfor = NULL;
//	printf("%f\n",((T *)mxGetData(MaskWeights))[31130]);

	//Free memory:
	for (int i = 0; i < inNum; i++) parameters[i] = NULL;
	mxFree(parameters); parameters = NULL;
	for (int i = 0; i < outNum; i++) results[i] = NULL;
	mxFree(results); results = NULL;
	mxDestroyArray(DoubleScalar); DoubleScalar = NULL;
	// OutPut:
	if (plhs[0] != NULL) mxDestroyArray(plhs[0]);
	if (plhs[1] != NULL) mxDestroyArray(plhs[1]);
	if (plhs[2] != NULL) mxDestroyArray(plhs[2]);
	if (plhs[3] != NULL) mxDestroyArray(plhs[3]);
	plhs[0] = DataIndices;
	plhs[1] = DataWeights;
	plhs[2] = MaskWeights;
	plhs[3] = GridSize;
}
template<typename T,typename L>
void GC_SL_Interface(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	typedef Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> ArrayXT;
	typedef Eigen::Array<L,Eigen::Dynamic,Eigen::Dynamic> ArrayXL;
	typedef Eigen::Map<ArrayXT> MapArrayXT;
	typedef Eigen::Map<ArrayXL> MapArrayXL;

	if (nlhs != 1) {printf("Input is not 1!");exit(-1);}
	if (nrhs != 6) {printf("Input is not 7!");exit(-1);}

	const mxArray *DataIndices = prhs[0];
	const mxArray *Labels = prhs[1];
	const mxArray *DataVisualCues = prhs[2];
	const mxArray *DataMotionCues = prhs[3];
	const mxArray *GridSize = prhs[4];
	const mxArray *VideoSize = prhs[5];

	mxArray *MotionCues = NULL;
	mxArray *VisualCues = NULL;
	mxArray *Cues = NULL;
	mxArray *Round = NULL;
	mxArray *WeightCeil = NULL;
	mxArray *WeightFloor = NULL;
	mxArray *SliceData = NULL;

	mxArray *DoubleScalar = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);

	int inNum = 6; int outNum = 5;
	const mxArray **parameters = (const mxArray **)mxCalloc(inNum,sizeof(const mxArray *));
	mxArray **results = (mxArray **)mxCalloc(outNum,sizeof(mxArray *));

	MapArrayXT *map = new MapArrayXT((T *)mxGetData(DataMotionCues),mxGetM(DataMotionCues),mxGetN(DataMotionCues));
	map->col(0).swap(map->col(2));
	delete map; map = NULL;

	parameters[0] = DataMotionCues; results[0] = MotionCues;
	mexRound<T>(1,results,1,parameters);
	MotionCues = results[0]; results[0] = NULL;

	parameters[0] = DataVisualCues; results[0] = VisualCues;
	mexFloor<T>(1,results,1,parameters);
	VisualCues = results[0]; results[0] = NULL;

	parameters[0] = VisualCues; parameters[1] = MotionCues; results[0] = Round;
	mexHorzcat<T>(1,results,2,parameters);
	Round = results[0]; results[0] = NULL;
	mxDestroyArray(MotionCues); mxDestroyArray(VisualCues);
	MotionCues = VisualCues = NULL;

	parameters[0] = DataVisualCues; parameters[1] = DataMotionCues; results[0] = Cues;
	mexHorzcat<T>(1,results,2,parameters);
	Cues = results[0]; results[0] = NULL;

	parameters[0] = Cues; parameters[1] = Round; results[0] = WeightCeil;
	mexSubstract<T>(1,results,2,parameters);
	WeightCeil = results[0]; results[0] = NULL;
	mxDestroyArray(Cues);
	Cues = NULL;

	((double *)mxGetData(DoubleScalar))[0] = 1;
	parameters[0] = DoubleScalar; parameters[1] = WeightCeil; results[0] = WeightFloor;
	mexSubstract<T>(1,results,2,parameters);
	WeightFloor = results[0]; results[0] = NULL;

	mxArray *Cell = mxCreateCellMatrix( 1, 1 );
	mxArray *data = mxCreateNumericMatrix(mxGetM(DataIndices),2,mxDOUBLE_CLASS,mxREAL);
	mxSetCell(Cell,0,data);

	map = new MapArrayXT((T*)mxGetData(data),mxGetM(DataIndices),2);
	MapArrayXT *col1 = new MapArrayXT((T*)mxGetData(DataIndices),mxGetM(DataIndices),1);
	MapArrayXL *col2 = new MapArrayXL((L*)mxGetData(Labels),mxGetM(Labels),1);
	map->col(0) = *col1;
	map->col(1) = col2->template cast<T>();
	delete map; delete col1; delete col2;
	map = col1 = NULL;
	col2 = NULL;

	((double *)mxGetData(DoubleScalar))[0] = (double)mxGetN(DataVisualCues);
	parameters[0] = Cell; parameters[1] = Round; parameters[2] = GridSize; parameters[3] = WeightFloor; parameters[4] = WeightCeil;
	parameters[5] = DoubleScalar; results[0] = SliceData;

	MexSlice_Adjacent(1,results,6,parameters);
	SliceData = results[0]; results[0] = NULL;
	mxDestroyArray(Cell); Cell = data = NULL;
	mxDestroyArray(Round); Round = NULL;
	mxDestroyArray(WeightFloor); WeightFloor = NULL;
	mxDestroyArray(WeightCeil); WeightCeil = NULL;

	// Reshape SliceData back into video
	if (plhs[0] != NULL) mxDestroyArray(plhs[0]);
	mwSize nDims = mxGetN(VideoSize);
	mwSize *sz = new mwSize[nDims];
	for (mwSize i = 0; i < nDims; i ++) sz[i] = (mwSize)((T *)mxGetData(VideoSize))[i];
	plhs[0] = mxCreateNumericArray(nDims,(const mwSize *)sz,mxGetClassID(DataIndices),mxREAL);
	CMemoryCopy<T,mwSize>((T *)mxGetData(plhs[0]), (T *)mxGetData(SliceData), mxGetNumberOfElements(SliceData));
	mxDestroyArray(SliceData); SliceData = NULL;
	delete sz; sz = NULL;

	// Free memory
	for (int i = 0; i < inNum; i++) parameters[i] = NULL;
	mxFree(parameters); parameters = NULL;
	for (int i = 0; i < outNum; i++) results[i] = NULL;
	mxFree(results); results = NULL;
	mxDestroyArray(DoubleScalar); DoubleScalar = NULL;
}

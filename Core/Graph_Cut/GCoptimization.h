#ifndef __GCOPTIMIZATION_H__
#define __GCOPTIMIZATION_H__

/////////////////////////////////////////////////////////////////////
// Inclusion of files and macros
/////////////////////////////////////////////////////////////////////

// Due to quiet bugs in function template specialization, it is not 
// safe to use earlier MS compilers. 
#if defined(_MSC_VER) && _MSC_VER < 1400
#error Requires Visual C++ 2005 (VC8) compiler or later.
#endif
#ifdef MATLAB_MEX_FILE
#include <mex.h>
#endif

#include <cstddef>
#include "energy.h"
#include "graph.h"
#include "maxflow.cpp"
#include "GCLinkedBlockList.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <string.h>


/////////////////////////////////////////////////////////////////////
// Utility functions, classes, and macros
/////////////////////////////////////////////////////////////////////

#ifdef _WIN32
typedef __int64 gcoclock_t;
#else
#include <ctime>
typedef clock_t gcoclock_t;
#endif
extern "C" gcoclock_t gcoclock(); // fairly high-resolution timer... better than clock() when available
extern "C" gcoclock_t GCO_CLOCKS_PER_SEC; // this variable will stay 0 until gcoclock() is called for the first time

#ifdef _MSC_EXTENSIONS
#define OLGA_INLINE __forceinline
#else
#define OLGA_INLINE inline
#endif

#ifndef GCO_MAX_ENERGYTERM
#define GCO_MAX_ENERGYTERM 10000000000  // maximum safe coefficient to avoid integer overflow
                                     	// if a data/smooth/label cost term is larger than this,
                                     	// the library will raise an exception
#define MAX_INTT 1000000
#endif

#if defined(GCO_ENERGYTYPE) && !defined(GCO_ENERGYTERMTYPE)
#define GCO_ENERGYTERMTYPE GCO_ENERGYTYPE
#endif
#if !defined(GCO_ENERGYTYPE) && defined(GCO_ENERGYTERMTYPE)
#define GCO_ENERGYTYPE GCO_ENERGYTERMTYPE
#endif

class GCException {
public:
	const char* message;
	GCException( const char* m ): message(m) { }
	void Report();
};

void GCException::Report()
{
	printf("\n%s\n",message);
	exit(0);
}

#define ErrMessage NULL

OLGA_INLINE static void handleError(const char *message) {throw GCException(message);};

//////////////////////////////////////////////////////////////////////////////////////////////////
// Define DataCostFn class and its derives
//////////////////////////////////////////////////////////////////////////////////////////////////

template <typename EnergyTermType, typename SiteID, typename LabelID>
class DataCostFn
{
public:
	typedef EnergyTermType (*DataCostFun)(SiteID s, LabelID l);
	typedef EnergyTermType (*DataCostFunExtra)(SiteID s, LabelID l,void *);

	virtual ~DataCostFn(){} ;
	virtual EnergyTermType compute(SiteID s, LabelID l) = 0;
};

template <typename EnergyTermType, typename SiteID, typename LabelID>
class DataCostFnFromArray : public DataCostFn<EnergyTermType, SiteID, LabelID>
{
public:
	DataCostFnFromArray(EnergyTermType* theArray, LabelID num_labels) : m_array(theArray), m_num_labels(num_labels){}
	~DataCostFnFromArray(){};
	OLGA_INLINE EnergyTermType compute(SiteID s, LabelID l){return m_array[s*m_num_labels+l];}
private:
	const EnergyTermType* const m_array;
	const LabelID m_num_labels;
};

//template <typename EnergyTermType, typename SiteID, typename LabelID>
//DataCostFnFromArray<EnergyTermType, SiteID, LabelID>::~DataCostFnFromArray() { if (m_array) /*delete m_array;*/ m_array = NULL; }

template <typename DataCostFun, typename EnergyTermType, typename SiteID, typename LabelID>
class DataCostFnFromFunction : public DataCostFn<EnergyTermType, SiteID, LabelID>
{
public:
	DataCostFnFromFunction(DataCostFun fn) : m_fn(fn){};
	OLGA_INLINE EnergyTermType compute(SiteID s, LabelID l){return m_fn(s,l);}
private:
	const DataCostFun m_fn;
};

template <typename DataCostFunExtra, typename EnergyTermType, typename SiteID, typename LabelID>
class DataCostFnFromFunctionExtra : public DataCostFn<EnergyTermType, SiteID, LabelID>
{
public:
	DataCostFnFromFunctionExtra(DataCostFunExtra fn,void *extraData) : m_fn(fn),m_extraData(extraData){};
	~DataCostFnFromFunctionExtra(){if (m_extraData) free(m_extraData);};
	OLGA_INLINE EnergyTermType compute(SiteID s, LabelID l){return m_fn(s,l,m_extraData);};
private:
	const DataCostFunExtra m_fn;
	void *m_extraData;
};
/////////////////////////////////////////////////////////////////////
// DataCostFnSparse
//   This data cost functor maintains a simple sparse structure
//   to quickly find the cost associated with any (site,label) pair.
/////////////////////////////////////////////////////////////////////

// Set cost of assigning 'l' to a specific subset of sites.
// The sites are listed as (SiteID,cost) pairs.

template <typename EnergyTermType, typename SiteID, typename LabelID>
class DataCostFnSparse : public DataCostFn<EnergyTermType, SiteID, LabelID>
{
public:
	typedef struct SparseDataCost {
		SiteID site;
		EnergyTermType cost;
	}SparseDataCost;
	typedef struct DataCostBucket {
		const SparseDataCost* begin;
		const SparseDataCost* end;     // one-past-the-last item in the range
		const SparseDataCost* predict; // predicts the next cost to be needed
	}DataCostBucket;
	// cLogSitesPerBucket basically controls the compression ratio
	// of the sparse structure: 1 => a dense array, num_sites => a single sparse list.
	// The amount (cLogSitesPerBucket - cLinearSearchSize) determines the maximum
	// number of binary search steps taken for a cost lookup for specific (site,label).

	static const int cLogSitesPerBucket = 9;
	static const int cSitesPerBucket = (1 << cLogSitesPerBucket);
	static const size_t    cDataCostPtrMask = ~(sizeof(SparseDataCost)-1);
	static const ptrdiff_t cLinearSearchSize = 64/sizeof(SparseDataCost);

public:
	DataCostFnSparse(SiteID num_sites, LabelID num_labels): m_num_sites(num_sites)
														  , m_num_labels(num_labels)
														  , m_buckets_per_label((m_num_sites + cSitesPerBucket-1)/cSitesPerBucket)
														  , m_buckets(0){};
	DataCostFnSparse(const DataCostFnSparse& src);
	~DataCostFnSparse();

	void           set(LabelID l, const SparseDataCost* costs, SiteID count);
	EnergyTermType compute(SiteID s, LabelID l);
	SiteID         queryActiveSitesExpansion(LabelID alpha_label, const LabelID* labeling, SiteID* activeSites);

	class iterator {
	public:
		OLGA_INLINE iterator(): m_ptr(0) { }
		OLGA_INLINE iterator& operator++() { m_ptr++; return *this; }
		OLGA_INLINE SiteID         site() const { return m_ptr->site; }
		OLGA_INLINE EnergyTermType cost() const { return m_ptr->cost; }
		OLGA_INLINE bool      operator==(const iterator& b) const { return m_ptr == b.m_ptr; }
		OLGA_INLINE bool      operator!=(const iterator& b) const { return m_ptr != b.m_ptr; }
		OLGA_INLINE ptrdiff_t operator- (const iterator& b) const { return m_ptr  - b.m_ptr; }
	private:
		OLGA_INLINE iterator(const SparseDataCost* ptr): m_ptr(ptr) { }
		const SparseDataCost* m_ptr;
		friend class DataCostFnSparse;
	};

	OLGA_INLINE iterator begin(LabelID label) const { return m_buckets[label*m_buckets_per_label].begin; }
	OLGA_INLINE iterator end(LabelID label)   const { return m_buckets[label*m_buckets_per_label + m_buckets_per_label-1].end; }

private:
	EnergyTermType search(DataCostBucket& b, SiteID s);
	const SiteID  m_num_sites;
	const LabelID m_num_labels;
	const int m_buckets_per_label;
	mutable DataCostBucket* m_buckets;
};

//////////////////////////////////////////////////////////////////////////////////////////////////
// Define SmoothCostFn class and its derives
//////////////////////////////////////////////////////////////////////////////////////////////////
template <typename EnergyTermType, typename SiteID, typename LabelID>
class SmoothCostFn
{
public:
	typedef EnergyTermType (*SmoothCostFun)(SiteID s1, SiteID s2, LabelID l1, LabelID l2);
	typedef EnergyTermType (*SmoothCostFunExtra)(SiteID s1, SiteID s2, LabelID l1, LabelID l2,void *);
	virtual ~SmoothCostFn(){};
	virtual EnergyTermType compute(SiteID s1, SiteID s2, LabelID l1, LabelID l2) = 0;
};

template <typename EnergyTermType, typename SiteID, typename LabelID>
class SmoothCostFnFromArray : public SmoothCostFn<EnergyTermType, SiteID, LabelID>
{
public:
	SmoothCostFnFromArray(EnergyTermType* theArray, LabelID num_labels) : m_array(theArray), m_num_labels(num_labels){}
	~SmoothCostFnFromArray(){};
	OLGA_INLINE EnergyTermType compute(SiteID s1, SiteID s2, LabelID l1, LabelID l2){return m_array[l1*m_num_labels+l2];}
private:
	const EnergyTermType* const m_array;
	const LabelID m_num_labels;
};

//template <typename EnergyTermType, typename SiteID, typename LabelID>
//SmoothCostFnFromArray<EnergyTermType, SiteID, LabelID>::~SmoothCostFnFromArray() {if (m_array) /*delete m_array;*/ m_array = NULL; }

template <typename SmoothCostFun, typename EnergyTermType, typename SiteID, typename LabelID>
class SmoothCostFnFromFunction : public SmoothCostFn<EnergyTermType, SiteID, LabelID>
{
public:
	SmoothCostFnFromFunction(SmoothCostFun fn) : m_fn(fn){}
	OLGA_INLINE EnergyTermType compute(SiteID s1, SiteID s2, LabelID l1, LabelID l2){return m_fn(s1,s2,l1,l2);}
private:
	const SmoothCostFun m_fn;
};

template <typename SmoothCostFunExtra, typename EnergyTermType, typename SiteID, typename LabelID>
class SmoothCostFnFromFunctionExtra : public SmoothCostFn<EnergyTermType, SiteID, LabelID>
{
public:
	SmoothCostFnFromFunctionExtra(SmoothCostFunExtra fn,void *extraData)	: m_fn(fn),m_extraData(extraData){}
	~SmoothCostFnFromFunctionExtra(){if (m_extraData) free(m_extraData);}
	OLGA_INLINE EnergyTermType compute(SiteID s1, SiteID s2, LabelID l1, LabelID l2){return m_fn(s1,s2,l1,l2,m_extraData);}
private:
	const SmoothCostFunExtra m_fn;
	void *m_extraData;
};

template <typename EnergyTermType, typename SiteID, typename LabelID>
class SmoothCostFnPotts : public SmoothCostFn<EnergyTermType, SiteID, LabelID>
{
public:
	OLGA_INLINE EnergyTermType compute(SiteID, SiteID, LabelID l1, LabelID l2){return l1 != l2 ? (EnergyTermType)1 : (EnergyTermType)0;}
};

/////////////////////////////////////////////////////////////////////
// GCoptimization class
/////////////////////////////////////////////////////////////////////

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
class GCoptimization
{
public:

	typedef DataCostFn<EnergyTermType, SiteID, LabelID> DataCostFnT;
	typedef SmoothCostFn<EnergyTermType, SiteID, LabelID> SmoothCostFnT;
	typedef Energy<EnergyTermType,EnergyTermType,EnergyType> EnergyT;
	typedef typename EnergyT::Var VarID;

	GCoptimization(SiteID num_sites, LabelID num_labels);
	virtual ~GCoptimization();

	/* Peforms expansion algorithm. Runs the number of iterations specified by max_num_iterations */
	/* Returns the total energy of labeling   */
	EnergyType expansion(int max_num_iterations = -1);

	/* Peforms one iteration (one pass over all labels)  of expansion algorithm.*/
	EnergyType oneExpansionIteration();

	/* Peforms  expansion on one label, specified by the input parameter alpha_label */
	void alpha_expansion(LabelID alpha_label);

	// /* Peforms  expansion on label alpha_label only for pixels specified by *pixels.  */
	//  num is the size of array pixels                                                
	// EnergyType alpha_expansion(LabelID alpha_label, SiteID *pixels, int num);

	/* Peforms swap algorithm. Runs it the specified number of iterations */
	EnergyType swap(int max_num_iterations = -1);

	/* Peforms one iteration (one pass over all labels)  of swap algorithm.*/
	EnergyType oneSwapIteration();

	/* Peforms  swap on a pair of labels, specified by the input parameters alpha_label, beta_label */
	void alpha_beta_swap(LabelID alpha_label, LabelID beta_label);

	void setupDataCostsExpansion(SiteID size,LabelID alpha_label,EnergyT *e,VarID *variables);

	void setupSmoothCostsExpansion(SiteID size,LabelID alpha_label,EnergyT *e,VarID *variables);

	void setupDataCostsSwap(SiteID size, LabelID alpha_label, LabelID beta_label, EnergyT *e,VarID *variables, SiteID *sites);

	void setupSmoothCostsSwap(SiteID size, LabelID alpha_label,LabelID beta_label, EnergyT *e,VarID *variables, SiteID *sites);

	// Check for overflow and submodularity issues when setting up binary graph cut
	void addterm1_checked(EnergyT *e,VarID i,EnergyTermType e0,EnergyTermType e1,EnergyTermType w);
	void addterm2_checked(EnergyT *e,VarID i,VarID j,EnergyTermType e00,EnergyTermType e01,EnergyTermType e10,EnergyTermType e11,EnergyTermType w);

	void permuteLabelTable();

	// Set cost for all (SiteID,LabelID) pairs. Default data cost is all zeros.
	void setDataCost(DataCostFnT *fn);
	
//	 void setDataCost(LabelID l, SparseDataCost *costs, SiteID count);

	// Set cost for all (LabelID,LabelID) pairs; the actual smooth cost is then weighted
	// at each pair of on neighbors. Defaults to Potts model (0 if l1==l2, 1 otherwise)
	void setSmoothCost(SmoothCostFnT *fn);

	/* Returns current label assigned to input pixel */
	OLGA_INLINE LabelID whatLabel(SiteID site){assert(site >= 0 && site < m_num_sites);return(m_labeling[site]);};
	/* copy the internal labels array into a given external array                  */
    OLGA_INLINE void whatLabel(SiteID start, SiteID count, LabelID* labeling) {assert(start >= 0 && start+count <= m_num_sites); memcpy(labeling, m_labeling+start, count*sizeof(LabelID));};

	// This function can be used to change the label of any site at any time      
	OLGA_INLINE void setLabel(SiteID site, LabelID label){assert(label >= 0 && label < m_num_labels && site >= 0 && site < m_num_sites);
	m_labeling[site] = label;};

	// setLabelOrder(false) sets the order to be not random; setLabelOrder(true) 
	//	sets the order to random. By default, the labels are visited in non-random order 
	//	for both the swap and alpha-expansion moves                         
	//	Note that srand() must be initialized with an appropriate seed in order for 
	//	random order to take effect!
	void setLabelOrder(bool isRandom);

	// Returns total energy for the current labeling
	EnergyType compute_energy();

	// Returns separate Data, Smooth, and Label energy of current labeling 
	EnergyType giveDataEnergy();
	EnergyType giveSmoothEnergy();
	// EnergyType giveLabelEnergy();

	// Returns number of sites/labels as specified in the constructor
	OLGA_INLINE SiteID  numSites() const {return m_num_sites;};
	OLGA_INLINE LabelID numLabels() const {return m_num_labels;};
	void readOffGraph(EnergyT *e, VarID *varibles, int num_sites);
	// Prints output to stdout during exansion/swap execution.
	//   0 => no output
	//   1 => cycle-level output (cycle number, current energy)
	//   2 => expansion-/swap-level output (label(s), current energy)
	// void setVerbosity(int level) { m_verbosity = level; }

	/* validate class function */
    // bool IsClassValid() { return class_sig == VALID_CLASS_SIGNITURE; }; 
    bool GetTruncateState() { return truncate_flag; };
    void SetTruncateState(bool tf) { truncate_flag = tf; };

protected:

	LabelID m_num_labels;
	SiteID  m_num_sites;
	SmoothCostFnT *m_smoothcostFn;	// base class pointer to smooth cost
	DataCostFnT *m_datacostFn; 		// base class pointer to data cost
	SiteID *m_numNeighbors;   		// holds num of neighbors for each site
	SiteID  m_numNeighborsTotal;    // holds total num of neighbor relationships
	SiteID  *m_lookupSiteVar; 		// holds index of variable corresponding to site participating in a move,
			                          		// -1 for nonparticipating site
	LabelID *m_labeling;
	LabelID *m_labelTable;    		// to figure out label order in which to do expansion/swaps


	// returns a pointer to the neighbors of a site and the weights
	virtual void giveNeighborInfo(SiteID site, SiteID *numSites, SiteID **neighbors, EnergyTermType **weights)=0;
	virtual void finalizeNeighbors() = 0;

private:

// 	int class_sig;  bagon: signiture value to verify class is ok 
	bool truncate_flag;
	bool m_random_label_order;
};

//////////////////////////////////////////////////////////////////////////////////////////////////
// GCoptimizationGridGraph classes derived from basic graphs
//////////////////////////////////////////////////////////////////////////////////////////////////
template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
class GCoptimizationGridGraph: public GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>
{
public:
	GCoptimizationGridGraph(SiteID width,SiteID height,LabelID num_labels);
	virtual ~GCoptimizationGridGraph();

	void setSmoothCostVH(EnergyTermType *smoothArray, EnergyTermType *vCosts, EnergyTermType *hCosts);

protected:
	virtual void giveNeighborInfo(SiteID site, SiteID *numSites, SiteID **neighbors, EnergyTermType **weights);
	virtual void finalizeNeighbors();
	EnergyTermType m_unityWeights[4];
	int m_weightedGraph;  // true if spatially varying w_pq's are present. False otherwise.

private:
	SiteID m_width;
	SiteID m_height;
	SiteID *m_neighbors;                 // holds neighbor indexes
	EnergyTermType *m_neighborsWeights;    // holds neighbor weights

	void setupNeighbData(SiteID startY,SiteID endY,SiteID startX,SiteID endX,SiteID maxInd,SiteID *indexes);
	void computeNeighborWeights(EnergyTermType *vCosts,EnergyTermType *hCosts);
};

//////////////////////////////////////////////////////////////////////////////////////////////////
// Functions for the GCoptimizationGridGraph, derived from GCoptimization
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
GCoptimizationGridGraph<EnergyTermType, EnergyType, SiteID, LabelID>::GCoptimizationGridGraph(SiteID width, SiteID height,LabelID num_labels) : GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>(width*height,num_labels)
{
	assert( (width > 1) && (height > 1) && (num_labels > 1 ));

	m_weightedGraph = 0;
	for (int  i = 0; i < 4; i ++ )	m_unityWeights[i] = 1;

	m_width  = width;
	m_height = height;

	this->m_numNeighbors = new SiteID[this->m_num_sites];
	this->m_neighbors = new SiteID[4*this->m_num_sites];

	SiteID indexes[4] = {-1,1,-m_width,m_width};

	SiteID indexesL[3] = {1,-m_width,m_width};
	SiteID indexesR[3] = {-1,-m_width,m_width};
	SiteID indexesU[3] = {1,-1,m_width};
	SiteID indexesD[3] = {1,-1,-m_width};

	SiteID indexesUL[2] = {1,m_width};
	SiteID indexesUR[2] = {-1,m_width};
	SiteID indexesDL[2] = {1,-m_width};
	SiteID indexesDR[2] = {-1,-m_width};

	setupNeighbData(1,m_height-1,1,m_width-1,4,indexes);

	setupNeighbData(1,m_height-1,0,1,3,indexesL);
	setupNeighbData(1,m_height-1,m_width-1,m_width,3,indexesR);
	setupNeighbData(0,1,1,width-1,3,indexesU);
	setupNeighbData(m_height-1,m_height,1,m_width-1,3,indexesD);

	setupNeighbData(0,1,0,1,2,indexesUL);
	setupNeighbData(0,1,m_width-1,m_width,2,indexesUR);
	setupNeighbData(m_height-1,m_height,0,1,2,indexesDL);
	setupNeighbData(m_height-1,m_height,m_width-1,m_width,2,indexesDR);
}

//-------------------------------------------------------------------

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
GCoptimizationGridGraph<EnergyTermType, EnergyType, SiteID, LabelID>::~GCoptimizationGridGraph()
{
	delete [] this->m_numNeighbors;
	if ( m_neighbors )
		delete [] m_neighbors;
	if (m_weightedGraph) delete [] m_neighborsWeights;
	delete this->m_datacostFn; this->m_datacostFn = NULL;
	delete this->m_smoothcostFn; this->m_smoothcostFn = NULL;
}


//-------------------------------------------------------------------

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
void GCoptimizationGridGraph<EnergyTermType, EnergyType, SiteID, LabelID>::setupNeighbData(SiteID startY,SiteID endY,SiteID startX,
											  SiteID endX,SiteID maxInd,SiteID *indexes)
{
	SiteID x,y,pix;
	SiteID n;

	for ( y = startY; y < endY; y++ )
		for ( x = startX; x < endX; x++ )
		{
			pix = x+y*m_width;
			this->m_numNeighbors[pix] = maxInd;
			this->m_numNeighborsTotal += maxInd;

			for (n = 0; n < maxInd; n++ )
				m_neighbors[pix*4+n] = pix+indexes[n];
		}
}

//-------------------------------------------------------------------

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
void GCoptimizationGridGraph<EnergyTermType, EnergyType, SiteID, LabelID>::finalizeNeighbors()
{
}

//-------------------------------------------------------------------

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
void GCoptimizationGridGraph<EnergyTermType, EnergyType, SiteID, LabelID>::setSmoothCostVH(EnergyTermType *smoothArray, EnergyTermType *vCosts, EnergyTermType *hCosts)
{
	setSmoothCost(smoothArray);
	m_weightedGraph = 1;
	computeNeighborWeights(vCosts,hCosts);
}

//-------------------------------------------------------------------

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
void GCoptimizationGridGraph<EnergyTermType, EnergyType, SiteID, LabelID>::giveNeighborInfo(SiteID site, SiteID *numSites, SiteID **neighbors, EnergyTermType **weights)
{
	*numSites  = this->m_numNeighbors[site];
	*neighbors = &m_neighbors[site*4];
	
	if (m_weightedGraph) *weights  = &m_neighborsWeights[site*4];
	else *weights = m_unityWeights;
}

//-------------------------------------------------------------------

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
void GCoptimizationGridGraph<EnergyTermType, EnergyType, SiteID, LabelID>::computeNeighborWeights(EnergyTermType *vCosts,EnergyTermType *hCosts)
{
	SiteID i,n,nSite;
	EnergyTermType weight;
	
	this->m_neighborsWeights = new EnergyTermType[this->m_num_sites*4];

	for ( i = 0; i < this->m_num_sites; i++ )
	{
		for ( n = 0; n < this->m_numNeighbors[i]; n++ )
		{
			nSite = m_neighbors[4*i+n];
			if ( i-nSite == 1 )            weight = hCosts[nSite];
			else if (i-nSite == -1 )       weight = hCosts[i];
			else if ( i-nSite == m_width ) weight = vCosts[nSite];
			else if (i-nSite == -m_width ) weight = vCosts[i];
	
			m_neighborsWeights[i*4+n] = weight;
		}
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////
// GCoptimizationGeneralGraph classes derived from basic graphs
//////////////////////////////////////////////////////////////////////////////////////////////////

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
class GCoptimizationGeneralGraph:public GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>
{
public:
	// This is the constructor for non-grid graphs. Neighborhood structure must  be specified by 
	// setNeighbors()  function
	GCoptimizationGeneralGraph(SiteID num_sites,LabelID num_labels);
	virtual ~GCoptimizationGeneralGraph();

	// Makes site1 and site2 neighbors of each other. Can be called only 1 time for each      
	// unordered pair of sites. Parameter weight can be used to set spacially varying terms     
	// If the desired penalty for neighboring sites site1 and site2 is                        
	// V(label1,label2) = weight*SmoothnessPenalty(label1,label2), then                        
	// member function setLabel should be called as: setLabel(site1,site2,weight)             
	void setNeighbors(SiteID site1, SiteID site2, EnergyTermType weight=1);

	// passes pointers to arrays storing neighbor information
	// numNeighbors[i] is the number of neighbors for site i
	// neighborsIndexes[i] is a pointer to the array storing the sites which are neighbors to site i
	// neighborWeights[i] is a pointer to array storing the weights between site i and its neighbors
	// in the same order as neighborIndexes[i] stores the indexes
	void setAllNeighbors(SiteID *numNeighbors,SiteID **neighborsIndexes,EnergyTermType **neighborsWeights);

protected: 
	virtual void giveNeighborInfo(SiteID site, SiteID *numSites, SiteID **neighbors, EnergyTermType **weights);
	virtual void finalizeNeighbors();

private:

	typedef struct NeighborStruct{
		SiteID  to_node;
		EnergyTermType weight;
	} Neighbor;

	GCLinkedBlockList *m_neighbors;
	bool m_needToFinishSettingNeighbors;
	SiteID **m_neighborsIndexes;
	EnergyTermType **m_neighborsWeights;
	bool m_needTodeleteNeighbors;
};

////////////////////////////////////////////////////////////////////////////////////////////////
// Functions for the GCoptimizationGeneralGraph, derived from GCoptimization
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
GCoptimizationGeneralGraph<EnergyTermType, EnergyType, SiteID, LabelID>::GCoptimizationGeneralGraph(SiteID num_sites,LabelID num_labels): GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>(num_sites,num_labels)
{
	assert( num_sites > 1 && num_labels > 1 );

	m_neighborsIndexes 		 = 0;
	m_neighborsWeights		 = 0;
	this->m_numNeighbors     = 0;
	m_neighbors       		 = 0;

	m_needTodeleteNeighbors        = true;
	m_needToFinishSettingNeighbors = true;
}

//------------------------------------------------------------------

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
GCoptimizationGeneralGraph<EnergyTermType, EnergyType, SiteID, LabelID>::~GCoptimizationGeneralGraph()
{
		
	if ( m_neighbors )
		delete [] m_neighbors;

	if ( this->m_numNeighbors && m_needTodeleteNeighbors )
	{
		for ( SiteID i = 0; i < this->m_num_sites; i++ )
		{
			if (this->m_numNeighbors[i] != 0 )
			{
				delete [] m_neighborsIndexes[i];
				delete [] m_neighborsWeights[i];
			}
		}

		delete [] this->m_numNeighbors;
		delete [] m_neighborsIndexes;
		delete [] m_neighborsWeights;
		delete this->m_datacostFn; this->m_datacostFn = NULL;
		delete this->m_smoothcostFn; this->m_smoothcostFn = NULL;
	}
}

//------------------------------------------------------------------

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
void GCoptimizationGeneralGraph<EnergyTermType, EnergyType, SiteID, LabelID>::finalizeNeighbors()
{
	if ( !m_needToFinishSettingNeighbors )
		return;
	m_needToFinishSettingNeighbors = false;

	Neighbor *tmp;
	SiteID i,site,count;

	EnergyTermType *tempWeights = new EnergyTermType[this->m_num_sites];
	SiteID *tempIndexes         = new SiteID[this->m_num_sites];
	
	if ( !tempWeights || !tempIndexes ) handleError("Not enough memory");

	this->m_numNeighbors     = new SiteID[this->m_num_sites];
	m_neighborsIndexes = new SiteID*[this->m_num_sites];
	m_neighborsWeights = new EnergyTermType*[this->m_num_sites];
	
	if ( !this->m_numNeighbors || !m_neighborsIndexes || !m_neighborsWeights ) handleError("Not enough memory.");

	for ( site = 0; site < this->m_num_sites; site++ )
	{
		if ( m_neighbors && !m_neighbors[site].isEmpty() )
		{
			m_neighbors[site].setCursorFront();
			count = 0;
			
			while ( m_neighbors[site].hasNext() )
			{
				tmp = (Neighbor *) (m_neighbors[site].next());
				tempIndexes[count] =  tmp->to_node;
				tempWeights[count] =  tmp->weight;
				delete tmp;
				count++;
			}
			this->m_numNeighbors[site]     = count;
			this->m_numNeighborsTotal     += count;
			m_neighborsIndexes[site] = new SiteID[count];
			m_neighborsWeights[site] = new EnergyTermType[count];
			
			if ( !m_neighborsIndexes[site] || !m_neighborsWeights[site] ) handleError("Not enough memory.");
			
			for ( i = 0; i < count; i++ )
			{
				m_neighborsIndexes[site][i] = tempIndexes[i];
				m_neighborsWeights[site][i] = tempWeights[i];
			}
		}
		else this->m_numNeighbors[site] = 0;

	}

	delete [] tempIndexes;
	delete [] tempWeights;
	if (m_neighbors) {
		delete [] m_neighbors;
		m_neighbors = 0;
	}
}
//------------------------------------------------------------------------------

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
void GCoptimizationGeneralGraph<EnergyTermType, EnergyType, SiteID, LabelID>::giveNeighborInfo(SiteID site, SiteID *numSites, 
												  SiteID **neighbors, EnergyTermType **weights)
{
	if (this->m_numNeighbors) {
		(*numSites)  =  this->m_numNeighbors[site];
		(*neighbors) = m_neighborsIndexes[site];
		(*weights)   = m_neighborsWeights[site];
	} else {
		*numSites = 0;
		*neighbors = 0;
		*weights = 0;
	}
}


//------------------------------------------------------------------

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
void GCoptimizationGeneralGraph<EnergyTermType, EnergyType, SiteID, LabelID>::setNeighbors(SiteID site1, SiteID site2, EnergyTermType weight)
{

	assert( site1 < this->m_num_sites && site1 >= 0 && site2 < this->m_num_sites && site2 >= 0);
	if ( m_needToFinishSettingNeighbors == false )
		handleError("Already set up neighborhood system.");

	if ( !m_neighbors )
	{
		this->m_neighbors = (GCLinkedBlockList *) new GCLinkedBlockList[this->m_num_sites];
		if ( !m_neighbors ) handleError("Not enough memory.");
	}

	Neighbor *temp1 = (Neighbor *) new Neighbor;
	Neighbor *temp2 = (Neighbor *) new Neighbor;

	temp1->weight  = weight;
	temp1->to_node = site2;

	temp2->weight  = weight;
	temp2->to_node = site1;

	m_neighbors[site1].addFront(temp1);
	m_neighbors[site2].addFront(temp2);
	
}
//------------------------------------------------------------------

//template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
//void GCoptimizationGeneralGraph<EnergyTermType, EnergyType, SiteID, LabelID>::setAllNeighbors(SiteID *numNeighbors,SiteID **neighborsIndexes,
//												 EnergyTermType **neighborsWeights)
//{
//	m_needTodeleteNeighbors = false;
//	m_needToFinishSettingNeighbors = false;
//	if ( m_numNeighborsTotal > 0 )
//		handleError("Already set up neighborhood system.");
//	m_numNeighbors     = numNeighbors;
//	m_numNeighborsTotal = 0;
//	for (int site = 0; site < m_num_sites; site++ ) m_numNeighborsTotal += m_numNeighbors[site];
//	m_neighborsIndexes = neighborsIndexes;
//	m_neighborsWeights = neighborsWeights;
//}

//-------------------------------------------------------------------
// GCoptimization methods
//-------------------------------------------------------------------

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>::GCoptimization(SiteID nSites, LabelID nLabels)
: m_num_labels(nLabels)
, m_num_sites(nSites)
, m_smoothcostFn(0)
, m_datacostFn(0)
, m_numNeighborsTotal(0)
, m_lookupSiteVar(new SiteID[nSites])
, m_labeling(new LabelID[nSites])
, m_labelTable(new LabelID[nLabels])
, truncate_flag(false)
{
	if ( nLabels <= 1 ) handleError("Number of labels must be >= 2");
	if ( nSites <= 0 )  handleError("Number of sites must be >= 1");

	if ( !m_lookupSiteVar || !m_labelTable || !m_labeling ){
		if (m_lookupSiteVar) delete [] m_lookupSiteVar;
		if (m_labelTable) delete [] m_labelTable;
		if (m_labeling) delete [] m_labeling;
		handleError("Not enough memory.");
	}
	for ( int i = 0; i < m_num_labels; i++ )
		m_labelTable[i] = i;
	memset(m_labeling, 0, m_num_sites*sizeof(LabelID));
	memset(m_lookupSiteVar, -1, m_num_sites*sizeof(SiteID));
	setLabelOrder(false);
}

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>::~GCoptimization()
{
	delete [] m_labelTable;
	delete [] m_lookupSiteVar;
	delete [] m_labeling;

	delete m_datacostFn; m_datacostFn = NULL;
	delete m_smoothcostFn; m_smoothcostFn = NULL;
}

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
void GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>::permuteLabelTable()
{
	if ( !m_random_label_order )
		return;
	for ( LabelID i = 0; i < m_num_labels; i++ )
	{
		LabelID j = i + (rand() % (m_num_labels-i));
		LabelID temp    = m_labelTable[i];
		m_labelTable[i] = m_labelTable[j];
		m_labelTable[j] = temp;
	}
}

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
void GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>::setLabelOrder(bool RANDOM_LABEL_ORDER)
{
	m_random_label_order = RANDOM_LABEL_ORDER;
}

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
EnergyType GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>::giveDataEnergy()
{
	EnergyType energy = 0;
	for ( SiteID i = 0; i < m_num_sites; i++ )
		energy += m_datacostFn->compute(i,m_labeling[i]);
	return energy;
}

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
EnergyType GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>::giveSmoothEnergy()
{
	EnergyType eng = (EnergyType) 0;
	SiteID i,numN,*nPointer,nSite,n;
	EnergyTermType *weights;
	for ( i = 0; i < m_num_sites; i++ )
	{
		giveNeighborInfo(i,&numN,&nPointer,&weights);
		for ( n = 0; n < numN; n++ )
		{
			nSite = nPointer[n];
			if ( nSite < i ) 
				eng += weights[n]*(m_smoothcostFn->compute(i,nSite,m_labeling[i],m_labeling[nSite]));
		}
	}

	return eng;
}

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
EnergyType GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>::compute_energy()
{
	// return giveDataEnergy() + giveSmoothEnergy() + giveLabelEnergy();
	return giveDataEnergy() + giveSmoothEnergy();
}

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
OLGA_INLINE void GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>::addterm1_checked(EnergyT* e, VarID i, EnergyTermType e0, EnergyTermType e1, EnergyTermType w = 1)
{
	if ( e0 > GCO_MAX_ENERGYTERM || e1 > GCO_MAX_ENERGYTERM )
		handleError("Smooth cost term was larger than GCO_MAX_ENERGYTERM; danger of integer overflow.");
	if ( w > GCO_MAX_ENERGYTERM )
		handleError("Smoothness weight was larger than GCO_MAX_ENERGYTERM; danger of integer overflow.");
	e->add_term1(i,e0*w,e1*w);
}

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
OLGA_INLINE void GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>::addterm2_checked(EnergyT* e, VarID i, VarID j, EnergyTermType e00, EnergyTermType e01, EnergyTermType e10, EnergyTermType e11, EnergyTermType w)
{
	if ( e00 > GCO_MAX_ENERGYTERM || e11 > GCO_MAX_ENERGYTERM || e01 > GCO_MAX_ENERGYTERM || e10 > GCO_MAX_ENERGYTERM )
		handleError("Smooth cost term was larger than GCO_MAX_ENERGYTERM; danger of integer overflow.");
	if ( w > GCO_MAX_ENERGYTERM )
		handleError("Smoothness weight was larger than GCO_MAX_ENERGYTERM; danger of integer overflow.");
	// Inside energy/maxflow code the submodularity check is performed as an assertion,
	// but is optimized out. We check it in release builds as well.
	if ( e00+e11 > e01+e10 )
		handleError("Non-submodular expansion term detected; smooth costs must be a metric for expansion");
	e->add_term2(i,j,e00*w,e01*w,e10*w,e11*w);
}

//-------------------------------------------------------------------//
//                  METHODS for EXPANSION MOVES                      //  
//-------------------------------------------------------------------//

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
OLGA_INLINE void GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>:: setDataCost(DataCostFnT *dc)
{
	if (!m_datacostFn) { delete m_datacostFn; m_datacostFn = NULL;}
	m_datacostFn = dc;
}

//template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
//OLGA_INLINE void GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>:: setDataCost(LabelID l, SparseDataCost *costs, SiteID count)
//{
//	m_datacostFn = dc;
//	m_datacostFn->set(l,costs,count);
//}

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
OLGA_INLINE void GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>:: setSmoothCost(SmoothCostFnT *sc)
{
	if (!m_smoothcostFn) { delete m_smoothcostFn; m_smoothcostFn = NULL;}
	m_smoothcostFn = sc;
}

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
OLGA_INLINE void GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>::setupDataCostsExpansion(SiteID size,LabelID alpha_label,EnergyT *e,VarID *variables)
{
	for ( SiteID i = 0; i < size; ++i )
		addterm1_checked(e,variables[i],m_datacostFn->compute(m_lookupSiteVar[i],alpha_label),m_datacostFn->compute(m_lookupSiteVar[i],m_labeling[m_lookupSiteVar[i]]));
}

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
OLGA_INLINE void GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>::setupSmoothCostsExpansion(SiteID size,LabelID alpha_label,EnergyT *e,VarID *variables)
{
	SiteID i,nSite,site,n,nNum,*nPointer;
	EnergyTermType *weights;

	for ( i = size - 1; i >= 0; i-- )
	{
		site = m_lookupSiteVar[i];
		m_lookupSiteVar[site] = i;
		giveNeighborInfo(site,&nNum,&nPointer,&weights);
		for ( n = 0; n < nNum; n++ )
		{
			nSite = nPointer[n];
			if ( m_labeling[nSite] != alpha_label )
			{	
				if ( site < nSite )
					addterm2_checked(e,variables[i],variables[m_lookupSiteVar[nSite]],
					                 m_smoothcostFn->compute(site,nSite,alpha_label,alpha_label),
					                 m_smoothcostFn->compute(site,nSite,alpha_label,m_labeling[nSite]),
					                 m_smoothcostFn->compute(site,nSite,m_labeling[site],alpha_label),
					                 m_smoothcostFn->compute(site,nSite,m_labeling[site],m_labeling[nSite]),weights[n]);
			}
			else
				addterm1_checked(e,variables[i],
								 m_smoothcostFn->compute(site,nSite,alpha_label,m_labeling[nSite]),
				                 m_smoothcostFn->compute(site,nSite,m_labeling[site],m_labeling[nSite]),weights[n]);
		}
	}
}

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
void GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>::alpha_expansion(LabelID alpha_label)
{
	finalizeNeighbors();
	SiteID i,size = 0; 
	EnergyT *e = new EnergyT(ErrMessage, truncate_flag);
	
	for ( i = 0; i < m_num_sites; i++ )
	{
		if ( m_labeling[i] != alpha_label )
		{
			m_lookupSiteVar[size] = i;
			size++;
		}
	}
	if ( size > 0 ) 
	{
		VarID *variables = new VarID[size];

		for ( i = 0; i < size; i++ )
		{
			variables[i] = e ->add_variable();
			e->setVarID(variables[i],m_lookupSiteVar[i]);
		}
		
		setupDataCostsExpansion(size, alpha_label, e, variables);
		
		setupSmoothCostsExpansion(size, alpha_label, e, variables);

//		printf("\n-------------------------------------------------------------------------\n");
//		printf("current labeling:\n");
//		printf("-------------------------------------------------------------------------\n");
//		for ( i = 0; i < m_num_sites; i++ ) {printf("%lld\t",m_labeling[i]+1); if((i+1)%10==0) printf("\n");}
//		printf("\n-------------------------------------------------------------------------\n");
//		printf("-------------------------------------------------------------------------\n");
//		printf("current expansion label: %lld\n", alpha_label+1);
//		printf("-------------------------------------------------------------------------\n");
//		printf("before minimization: datacost -> %1.0f\tsmoothcost -> %1.0f\n", giveDataEnergy(), giveSmoothEnergy());

		e -> Init();

//		readOffGraph(e,variables,size);

		e -> minimize();
	
		// lookuptable is a smart way to map all labels differing from alpha_label!
		for ( i = 0,size = 0; i < m_num_sites; i++ )
		{
			if ( m_labeling[i] != alpha_label )
			{
				if ( e->get_var(variables[size]) == 0 )
					m_labeling[i] = alpha_label;

				size++;
			}
		}

//		printf("\n-------------------------------------------------------------------------\n");
//		printf("current labeling:\n");
//		printf("-------------------------------------------------------------------------\n");
//		for ( i = 0; i < m_num_sites; i++ ) {printf("%lld\t",m_labeling[i]+1); if((i+1)%10==0) printf("\n");}
//		printf("\n-------------------------------------------------------------------------\n");
//		printf("after minimization: datacost -> %1.0f\tsmoothcost -> %1.0f\n",giveDataEnergy(), giveSmoothEnergy());
//		readOffGraph(e,variables,size);

		delete variables;
	}

	delete e;
}

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
EnergyType GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>::oneExpansionIteration()
{
	permuteLabelTable();
	for (LabelID next = 0;  next < m_num_labels;  next++ )
		alpha_expansion(m_labelTable[next]);
	
	return compute_energy();
}

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
EnergyType GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>::expansion(int max_num_iterations)
{
	if (max_num_iterations == -1) max_num_iterations = MAX_INTT;
	int curr_cycle = 1;
	EnergyType new_energy,old_energy;

	new_energy = compute_energy();

	old_energy = (new_energy+1)*10; // BAGON changed init value to exceed current energy by factor of 10 (thanks to A. Khan)

	while ( old_energy > new_energy  && curr_cycle <= max_num_iterations)
	{
		old_energy = new_energy;
		new_energy = oneExpansionIteration();	
		curr_cycle++;	
	}

	return new_energy;
}

//-------------------------------------------------------------------//
//                  METHODS for SWAP MOVES                           //  
//-------------------------------------------------------------------//

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
void GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>::setupDataCostsSwap(SiteID size, LabelID alpha_label, LabelID beta_label, EnergyT *e,VarID *variables, SiteID *sites)
{
	for ( SiteID i = 0; i < size; i++ )
		e->add_term1(variables[i],m_datacostFn->compute(sites[i],alpha_label),
		               m_datacostFn->compute(sites[i],beta_label) );
}

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
void GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>::setupSmoothCostsSwap(SiteID size, LabelID alpha_label, LabelID beta_label, EnergyT *e,VarID *variables, SiteID *sites)
{
	SiteID i,nSite,site,n,nNum,*nPointer;
	EnergyTermType *weights;

	for ( i = size - 1; i >= 0; i-- )
	{
		site = m_lookupSiteVar[i];
		m_lookupSiteVar[site] = i;
		giveNeighborInfo(site,&nNum,&nPointer,&weights);
		for ( n = 0; n < nNum; n++ )
		{
			nSite = nPointer[n];
			if ( m_labeling[nSite] == alpha_label || m_labeling[nSite] == beta_label )
			{
				if (site < nSite)
					addterm2_checked(e,variables[i],m_lookupSiteVar[nSite],
					                 m_smoothcostFn->compute(site,nSite,alpha_label,alpha_label),
					                 m_smoothcostFn->compute(site,nSite,alpha_label,beta_label),
				    	             m_smoothcostFn->compute(site,nSite,beta_label,alpha_label),
				        	         m_smoothcostFn->compute(site,nSite,beta_label,beta_label),weights[n]);
			}
			else
				addterm1_checked(e,variables[i],
								 m_smoothcostFn->compute(site,nSite,alpha_label,m_labeling[nSite]),
				                 m_smoothcostFn->compute(site,nSite,beta_label, m_labeling[nSite]),weights[n]);
		}
	}
}

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
EnergyType GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>::swap(int max_num_iterations)
{
	if (max_num_iterations == -1) max_num_iterations = MAX_INTT;

	int curr_cycle = 1;
	EnergyType new_energy,old_energy;

	new_energy = compute_energy();
	old_energy = (new_energy+1)*10; // BAGON changed init value to exceed current energy by factor of 10 (thanks to A. Khan)  

	while ( old_energy > new_energy  && curr_cycle <= max_num_iterations)
	{
		old_energy = new_energy;
		new_energy = oneSwapIteration();	
		curr_cycle++;	
	}

	return(new_energy);
}

//--------------------------------------------------------------------------------

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
EnergyType GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>:: GCoptimization::oneSwapIteration()
{
	permuteLabelTable();

	for (LabelID next = 0;  next < m_num_labels;  next++ )
		for (LabelID next1 = m_num_labels - 1;  next1 >= 0;  next1-- )
			if ( m_labelTable[next] < m_labelTable[next1] )
				alpha_beta_swap(m_labelTable[next],m_labelTable[next1]); 

	return(compute_energy());
}

//---------------------------------------------------------------------------------

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
void GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>::alpha_beta_swap(LabelID alpha_label, LabelID beta_label)
{
	SiteID i,size = 0;
	EnergyT *e = new EnergyT("", false); // BAGON - no truncation for ab-swap
	SiteID *pixels = new SiteID[m_num_sites];
	

	for ( i = 0; i < m_num_sites; i++ )
	{
		if ( m_labeling[i] == alpha_label || m_labeling[i] == beta_label)
		{
			pixels[size]    = i;
			m_lookupSiteVar[i] = size;
			size++;
		}
	}

	if ( size == 0 )
	{
		delete e;
		delete [] pixels;
		return;
	}

	// Create binary variables for each remaining site, add the data costs,
	// and compute the smooth costs between variables.
	VarID *variables = new VarID[size];
	for ( i = 0; i < size; i++ )
		variables[i] = e ->add_variable();

	setupDataCostsSwap(size, alpha_label, beta_label, e, variables);
	setupSmoothCostsSwap(size, alpha_label, beta_label, e, variables);
	e -> minimize();

	for ( i = 0; i < size; i++ )
		if ( e->get_var(variables[i]) == 0 )
			m_labeling[pixels[i]] = alpha_label;
		else m_labeling[pixels[i]] = beta_label;


	delete [] variables;
	delete [] pixels;
	delete e;
}

template <typename EnergyTermType, typename EnergyType, typename SiteID, typename LabelID>
void GCoptimization<EnergyTermType, EnergyType, SiteID, LabelID>::readOffGraph(EnergyT *e, VarID *varibles, int num_sites)
{
	for (int i = 0; i < num_sites; i ++)
		e->readVar(varibles[i]);
}

//-------------------------------------------------------------------
// DataCostFnSparse methods
//-------------------------------------------------------------------

template <typename EnergyTermType, typename SiteID, typename LabelID>
DataCostFnSparse<EnergyTermType, SiteID, LabelID>::DataCostFnSparse(const DataCostFnSparse& src)
: m_num_sites(src.m_num_sites)
, m_num_labels(src.m_num_labels)
, m_buckets_per_label(src.m_buckets_per_label)
, m_buckets(0)
{
	assert(!src.m_buckets); // not implemented
}

template <typename EnergyTermType, typename SiteID, typename LabelID>
DataCostFnSparse<EnergyTermType, SiteID, LabelID>::~DataCostFnSparse()
{
	if (m_buckets) {
		for (LabelID l = 0; l < m_num_labels; ++l)
			if (m_buckets[l*m_buckets_per_label].begin)
				delete [] m_buckets[l*m_buckets_per_label].begin;
		delete [] m_buckets;
	}
}

template <typename EnergyTermType, typename SiteID, typename LabelID>
void DataCostFnSparse<EnergyTermType, SiteID, LabelID>::set(LabelID l, const SparseDataCost* costs, SiteID count)
{
	// Create the bucket if necessary, and copy all the costs
	//
	if (!m_buckets) {
		m_buckets = new DataCostBucket[m_num_labels*m_buckets_per_label];
		memset(m_buckets, 0, m_num_labels*m_buckets_per_label*sizeof(DataCostBucket));
	}

	DataCostBucket* b = &m_buckets[l*m_buckets_per_label];
	if (b->begin)
		delete [] b->begin;
	SparseDataCost* next = new SparseDataCost[count];
	memcpy(next,costs,count*sizeof(SparseDataCost));

	//
	// Scan the list of costs and remember pointers to delimit the 'buckets', i.e. where 
	// ranges of SiteIDs lie along the array. Buckets can be empty (begin == end).
	//
	const SparseDataCost* end  = next+count;
	SiteID prev_site = -1;
	for (int i = 0; i < m_buckets_per_label; ++i) {
		b[i].begin = b[i].predict = next;
		SiteID end_site = (i+1)*cSitesPerBucket;
		while (next < end && next->site < end_site) {
			if (next->site < 0 || next->site >= m_num_sites)
				throw GCException("Invalid site id given for sparse data cost; must be within range.");
			if (next->site <= prev_site)
				throw GCException("Sparse data costs must be sorted in increasing order of SiteID");
			prev_site = next->site;
			++next;
		}
		b[i].end = next;
	}
}

template <typename EnergyTermType, typename SiteID, typename LabelID>
EnergyTermType DataCostFnSparse<EnergyTermType, SiteID, LabelID>::search(DataCostBucket& b, SiteID s)
{
	// Perform binary search for requested SiteID
	//
	const SparseDataCost* L = b.begin;
	const SparseDataCost* R = b.end-1;
	if ( R - L == m_num_sites )
		return b.begin[s].cost; // special case: this particular label is actually dense
	do {
		const SparseDataCost* mid = (const SparseDataCost*)((((size_t)L+(size_t)R) >> 1) & cDataCostPtrMask);
		if (s < mid->site)
			R = mid-1;         // eliminate upper range
		else if (mid->site < s)
			L = mid+1;         // eliminate lower range
		else {
			b.predict = mid+1;
			return mid->cost;  // found it!
		}
	} while (R - L > cLinearSearchSize);
	
	// Finish off with linear search over the remaining elements
	//
	do {
		if (L->site >= s) {
			if (L->site == s) {
				b.predict = L+1;
				return L->cost;
			}
			break;
		}
	} while (++L <= R);
	b.predict = L;

	return GCO_MAX_ENERGYTERM; // the site belongs to this bucket but with no cost specified
}

template <typename EnergyTermType, typename SiteID, typename LabelID>
EnergyTermType DataCostFnSparse<EnergyTermType, SiteID, LabelID>::compute(SiteID s, LabelID l)
{
	DataCostBucket& b = m_buckets[l*m_buckets_per_label + (s >> cLogSitesPerBucket)];
	if (b.begin == b.end)
		return GCO_MAX_ENERGYTERM;
	if (b.predict < b.end) {
		// Check for correct prediction
		if (b.predict->site == s)
			return (b.predict++)->cost; // predict++ for next time
		
		// If the requested site is missing from the site ids near 'predict'
		// then we know it doesn't exist in the bucket, so return INF
		if (b.predict->site > s && b.predict > b.begin && (b.predict-1)->site < s)
			return GCO_MAX_ENERGYTERM;
	}
	if ( (size_t)b.end - (size_t)b.begin == cSitesPerBucket*sizeof(SparseDataCost) )
		return b.begin[s-b.begin->site].cost; // special case: this particular bucket is actually dense!

	return search(b,s);
}

#endif

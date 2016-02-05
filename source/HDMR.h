/*
High Dimensional Model Representation
======================================
This class is intended to be used in conjunction with the AdaptiveSparseGrid library of "Xiang Ma, xm25@cornell.edu".
The interface of the two class are through SGwrite & SGread.
Status Note: Work in progress.
aryan.eftekhari@usi.ch/gmail.com
*/

#ifndef HDMR_h_is_included
#define HDMR_h_is_included

#include "SparseGrid.h"

class HDMR {
  public:

	// Properties
	unsigned int gridPointCount;

	void (*problemFunc)(double*, int, double*);

	SGread** sgread;
	SGwrite** sgwrite;

	typedef vector<double> DataContainerDouble;
	typedef map <string, DataContainerDouble>  ObjDouble;
	typedef map <string, int>  ArrInt;
	typedef map <string, vector<double>> var;

	struct ProbParam {
		int dim;
		int dof;
	} probParam;

	struct SgParam {
		int maxLevel;
		int gridType;
		double cutOff;
	} sgParam;

	struct HdmrParam {
		int maxOrder;
		int sampleSize;
		double cutOff;
	} hdmrParam;

	struct ComputePool {
		int rank = -1;
		int size = -1;

		int blockID = 0;
		int blockSize = 1;

		int threadID = 0;
		int tpb = 1;
		int opt_tpb = 1;

		vector<int> workload;

		MPI_Comm subComm = MPI_COMM_WORLD;
	} computePool;

	struct Job {
		int* task;
		int size;
	}* job;

	struct Data {
		double* xBar;
		double* fxBar;
	} data;

	struct RunTime {
		int verbose = 0;
		int writeState = 1;
		double start = 0;
		double end = 0;
	} runTime;

	struct Cache {
		ArrInt taskLookUpIndex;
		var interpolate;
		var integral;
		std::vector <string> blackList;
	} localCache;

	// Methods
	HDMR(int verbose = 0);
	~HDMR();



	/**
	* [HDMR::write Decompose the problem into a series of surplus files]
	*/
	int write(
	    void (*problemFunc_)(double*, int, double*), int problemDim_, int problemDoF_ ,
	    int maxOrder_ , int sampleSize_ , double hdmrCutOff_ ,
	    int maxLevel_ , int gridType_ , double sgCutOff_ ,
	    string folderName = "surplus/", vector<double>preDefinedXbar = {});

	void read(string folderName);

	void interpolate(double* xSet, double* fvalueSet, int pointCount = 1);


	/**
	* [HDMR::debug Print out current system information]
	*/
	void debug(int showComputePool = 1, int showProb = 1, int showSG = 1, int showHDMR = 1, int showJob = 1 , int showOther = 1);

  private:
	const char* hline = "\n=================================================================\n";

	/**
	 * Initialize process allocation
	 * @param computeArch            [Architecture of computePool]
	 * @param optimalThreadsPerBlock [Optimal Number of threads per block (Default =2)]
	 */
	void computeAlloc(string computeArch = "flat", int optimalThreadsPerBlock = 2);

	/**
	* [HDMR::xbar Estimate the value of argMin_{xbar} ||f(xbar) - 1/N sum_{N} f(x_n)||_2 for x_n in R^n [0,1] ]
	* @param numRuns       [Number of runs of Monte Carlo simulation]
	* @param dynamicDomain [Incrementally shrink the sample domain R^n [0,1], increase convergence rate (Experimental!!!!)]
	*/
	void xbar();
	/**
	* [HDMR::calc_xbar Generate the samples and for each xBar run]
	* @param lower_bound [Lower bound of search domain]
	* @param upper_bound [Upper bound of search domain]
	*/
	void calc_xbar(double lower_bound, double upper_bound);

	/**
	* [HDMR::jobAlloc Allocate the job to be processed]
	*/
	void writeJobAlloc();
	int readJobAlloc(string folderName = "surplus/");


	inline void traverseIterpolate(vector<int> alphabet, int k, double* fvalue, int op);

	/**
	 * Traverse the integral tree
	 * @param  alphabet [Alphabet set]
	 * @param  k        [Subset size]
	 * @param  fvalue   [Integral funciton value]
	 * @param  op       [Operation]
	 * @return          [1=ignore this alphabet]
	 */
	inline int  traverseIntegral  (vector<int> alphabet, int k, double* fvalue, int op);


	/**
	* [HDMR::combination (Wrap for calc_combination) Generate the unique combinations]
	* @param n      [Total set length]
	* @param k      [Selection length]
	* @param rtnBuffer [Buffer to write results to]
	*/
	void combination(vector<int> alphabet, int k, int* rtnBuffer);
	/**
	* [HDMR::calc_combination Generate the unique combinations]
	* @param alphabet   [Alphabet if combinations]
	* @param buffer     [Buffer to write results to]
	* @param order      [Selection length or order]
	* @param tempBuffer [Temporary buffer for min calculations]
	* @param offset     [Offset for recursive calls]
	* @param start      [Starting point of alphabet]
	* @param index      [Index of depth]
	*/
	void calc_combination(vector<int> alphabet, int* rtnBuffer, const int order, int* tempBuffer, int* offset, int start = 0, int index = 0);


	/**
	* [HDMR::nCk Calculate number of possible combinations]
	* @param  n [Total set length]
	* @param  k [Selection length]
	* @return   [Count of unique combinations Note! (1,2)==(2,1)]
	*/
	unsigned int nCk(int n, int k);

	double l2norm(double* x, int size);

	void cleanUp();
};

#endif
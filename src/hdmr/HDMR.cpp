/*
High Dimensional Model Representation
======================================
This class is intended to be used in conjunction with the AdaptiveSparseGrid library of "Xiang Ma, xm25@cornell.edu".
The interface of the two class are through SGwrite & SGread.
Status Note: Work in progress.
aryan.eftekhari@usi.ch/gmail.com
*/

#include "SGwrite.h"
#include "SGread.h"
#include "HDMR.h"
#include "util/combination.h"
#include "util/stringchop.h"
#include "util/linalg.h"
#include "util/print.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <random>
#include <string>
#include <map>
#include <algorithm>
#include <unistd.h> // Sleep
#include <unordered_map>


using std::cout;
using std::string;
using std::vector;
using std::map;
using std::fill;
using std::unordered_map;
using namespace std_plus;



/*
 * Constructur, Destructor and clear
*/
HDMR::HDMR(int verbose_) {
	runTime.verbose = verbose_;
}
HDMR::~HDMR() {
	clear();
	delete[] sgread;
}
void HDMR::clear() {

	MPI_Barrier(MPI_COMM_WORLD);

	probParam = {};
	sgParam = {};
	hdmrParam = {};
	computePool = {};
	job = {};
	runTime = {};
	cache.fval.clear();
	cache.comFun.clear();
}








/*
 * Write
 */
int HDMR::write( //  For HDMR
    void (*problemFunc_)(double*, int, double*),
    int problemDim_,
    int problemDoF_ ,
    int sgMaxLevel_ ,
    double sgCutOff_ ,
    int sgGridType_ ,
    int hmdrMaxOrder_ ,
    double hdmrCutOff_ ,
    vector<double>xBar_,
    int processPerNode_,
    string folderName_) {

	// Clear all variables
	clear();

	// Check for consistency between maxOrder and probelm Dim
	if (hmdrMaxOrder_ > problemDim_) {
		err("ERROR HDMR MaxOrder cannot be > Problem Dim");
	}

	// Set runTime mode and folder name
	folderName = folderName_;
	runTime.mode = "HDMR_WRITE";
	runTime.interpolationPoints = 0;

	// Start timer
	startTimer();

	// Initialize Problem Parameters
	problemFunc          = problemFunc_;
	probParam.dim        = problemDim_;
	probParam.dof        = problemDoF_;

	// Initialize Sparse grid parameters
	sgParam.maxLevel     = sgMaxLevel_;
	sgParam.gridType     = sgGridType_;
	sgParam.cutOff       = sgCutOff_;

	// Initialize HDMR parameters
	hdmrParam.maxOrder   = hmdrMaxOrder_;
	hdmrParam.cutOff     = hdmrCutOff_;
	hdmrParam.xBar.resize(probParam.dim);
	hdmrParam.fxBar.resize(probParam.dof);
	hdmrParam.xBar       = xBar_;
	problemFunc(&hdmrParam.xBar[0], probParam.dim, &hdmrParam.fxBar[0]);

	// Set procsss per pool, the default is set to a single node
	if (processPerNode_ == -1) {
		MPI_Comm_size(MPI_COMM_WORLD, &computePool.processPerNode);
	} else {
		computePool.processPerNode = processPerNode_;
	}

	// Allocate compute resources
	allocCompute();

	// Create working director
	if (computePool.grank == 0) {
		string command = "";

		command = "rm -rf " + folderName + " > nul 2>&1"; // Supress ouput
		system(command.c_str());

		command = "mkdir " + folderName + " > nul 2>&1"; // Supress ouput
		system(command.c_str());
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// Initlize sparse grids on all processes
	sgwrite = new SGwrite*[hdmrParam.maxOrder];
	for (int d = 1; d <= hdmrParam.maxOrder; ++d) {
		sgwrite[d - 1] = new SGwrite(problemFunc, d, probParam.dof , sgParam.maxLevel, sgParam.cutOff, sgParam.gridType , runTime.verbose);
		sgwrite[d - 1]->setAnchor(hdmrParam.xBar);
		sgwrite[d - 1]->resetMPI(computePool.subComm);
	}

	// Allocate full job on all processes
	allocJobs();

	// Generate the active jobs to be processed
	if (hdmrParam.cutOff == 0.0) {
		allocActiveJobs(0);
	} else {
		allocActiveJobs(1);
	}

	// Generate lower dimension surplus files for all active Jobs
	vector<int> activeDim;
	string fileName = "";
	for (int d = 1; d <= hdmrParam.maxOrder ; ++d) {

		// Resize active Dim to current hdmr order
		activeDim.resize(d);

		// Go thorugh job.active list
		for (int i = 0; i < job.active[d].size(); ++i) {

			// Round robin taks allocation
			if (i % computePool.nodeSize ==  computePool.nodeRank) {

				// Get active dimensions and set file path
				activeDim = job.active[d][i];
				fileName = folderName + vector_join(activeDim, ".") + ".data";

				//Build write and cleanup sparse grid
				runTime.interpolationPoints += sgwrite[d - 1]->build(activeDim);
				sgwrite[d - 1]->write(fileName);
				sgwrite[d - 1]->Cleanup();
			}
		}
	}

	// Sum up all interpolation points accross all process
	if (computePool.rank != 0) {
		runTime.interpolationPoints = 0;
	}
	MPI_Allreduce(MPI_IN_PLACE, &runTime.interpolationPoints, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

	// Write index and log file
	writeIndexFile();

	// Delete all instance of sgwrite
	MPI_Barrier(MPI_COMM_WORLD);
	delete[] sgwrite;

	// End timer
	endTimer();

	return runTime.interpolationPoints;
}
int HDMR::write( // For SG
    void (*problemFunc_)(double*, int, double*),
    int problemDim_,
    int problemDoF_ ,
    int sgMaxLevel_ ,
    double sgCutOff_ ,
    int sgGridType_ ,
    string folderName_) {

	// Clear all variables
	clear();

	// Set runTime mode and folder name
	folderName = folderName_;
	runTime.mode = "SG_WRITE";
	runTime.interpolationPoints = 0;

	// Start timer
	startTimer();

	//Initialize Problem Parameters
	problemFunc          = problemFunc_;
	probParam.dim        = problemDim_;
	probParam.dof        = problemDoF_;

	//Initialize Sparse grid parameters
	sgParam.maxLevel     = sgMaxLevel_;
	sgParam.gridType     = sgGridType_;
	sgParam.cutOff       = sgCutOff_;

	// Set procsss per pool - there is only 1 node with all process inside it!
	MPI_Comm_size(MPI_COMM_WORLD, &computePool.processPerNode);

	// Allocate compute scheme
	allocCompute();

	// Create working director
	if (computePool.grank == 0) {
		string command = "";

		command = "rm -rf " + folderName + " > nul 2>&1"; // Supress ouput
		system(command.c_str());

		command = "mkdir " + folderName + " > nul 2>&1"; // Supress ouput
		system(command.c_str());
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// Create single SG instance and compute
	sgwrite = new SGwrite*[1];
	sgwrite[0] = new SGwrite(problemFunc, probParam.dim, probParam.dof , sgParam.maxLevel, sgParam.cutOff, sgParam.gridType , runTime.verbose);
	sgwrite[0]->resetMPI(computePool.subComm);
	runTime.interpolationPoints = sgwrite[0]->build();
	sgwrite[0]->write(folderName + "surplus.data");

	// Write index and log file
	writeIndexFile();

	// Delete single instance of sgwrite
	MPI_Barrier(MPI_COMM_WORLD);
	delete sgwrite[0];

	// End timer
	endTimer();

	return runTime.interpolationPoints;
}








/*
* Read
*/
int HDMR::read(string folderName_) {

	// Clear Variables
	clear();

	// Initilize value
	folderName = folderName_;
	runTime.interpolationPoints = 0;

	// Start Timer
	startTimer();

	// Load index file
	readIndexFile();

	// Allocate parralelization scheme
	allocCompute();

	if (runTime.mode == "SG_READ") {

		// Initilize SGread instance on all processess
		sgread = new SGread*[1];
		sgread[0] = new SGread(runTime.verbose);
		sgread[0]->resetMPI(computePool.subComm);
		sgread[0]->read(folderName + "surplus.data");
	} else {

		// Initilize SGread instance on all processes for every active job
		sgread = new SGread*[job.active_size];
		for (int k = 0; k < job.active_size; ++k) {
			sgread[k] = new SGread(runTime.verbose);
			sgread[k]->resetMPI(computePool.subComm);
		}

		int k = 0;
		for (int d = 1; d < job.active.size(); ++d) {

			// Go thorugh job.active list
			for (int i = 0; i < job.active[d].size(); ++i) {

				// Round robin task allocation --- NOT DOING THIS
				//	if (i % computePool.nodeSize ==  computePool.nodeRank) {

				sgread[k]->read(folderName + vector_join(job.active[d][i], ".") + ".data");
				k++;
			}
		}
	}

	// End timer
	endTimer();

	MPI_Barrier(MPI_COMM_WORLD);
	return runTime.interpolationPoints;
}









/*
* Interpolation
*/
void HDMR::interpolate(double* xSet, double* valSet , int pointCount) {

	// Start Timer
	startTimer();

	// Force all valSet to equal zero
	fill(valSet, valSet + probParam.dof * pointCount, 0.0);
	MPI_Barrier(MPI_COMM_WORLD);

	// Call different interoplation function based on curret runtime mode
	if (runTime.mode == "SG_READ") {
		interpolate_SG(&xSet[0], &valSet[0], pointCount);
	} else if (runTime.mode == "HDMR_READ") {
		interpolate_HDMR(&xSet[0], &valSet[0], pointCount);
	} else {
		err("ERROR Invalid interpolation runTime.mode " + runTime.mode);
	}

	// End timer
	endTimer();
	MPI_Barrier(MPI_COMM_WORLD);
}
void HDMR::interpolate_SG(double* xSet, double* valSet , int pointCount) {

	// Cycle through each point
	int xShift, valShift;
	for (int i = 0; i < pointCount; ++i) {

		// round robin style task allocation
		if (i % computePool.nodeSize ==  computePool.nodeRank) {
			xShift = i * probParam.dim;
			valShift = i * probParam.dof;
			sgread[0]->interpolateValue(&xSet[xShift], &valSet[valShift]);
		}
	}

	// Reduce in place all tasks
	MPI_Allreduce(MPI_IN_PLACE, &valSet[0], probParam.dof * pointCount, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}
void HDMR::interpolate_HDMR(double* xSet, double* valSet , int pointCount) {

	int xShift, valShift;
	vector<int> alphabet;
	vector<int> sub_alphabet;

	// Allocate cache
	allocCache();

	// Cycle through each point
	for (int r = 0; r < pointCount; ++r) {

		if (r % computePool.nodeSize ==  computePool.nodeRank) {

			// Set shifting constants
			xShift = r * probParam.dim;
			valShift = r * probParam.dof;

			// Add the final -f_0 for at the start;
			linalg_add(&valSet[valShift], hdmrParam.fxBar);

			// Allocate the cache.fval on current process
			allocCacheFval(&xSet[xShift]);

			// *** Kernel of the code!
			// Calculate Component functions

			// Cycel through jobs
			for (int d = 1; d < job.active.size(); ++d) {

				alphabet.resize(d);
				for (int i = 0; i < job.active[d].size(); ++i) {

					// Initilize alphabit at the activeJob and set the component function intial value
					alphabet = job.active[d][i] ;
					cache.comFun [ alphabet ] = cache.fval[job.active[d][i]];

					// Generate all combinations of lower level combination of alphabit
					for (int k = 1; k < d; ++k) {

						// Intilize sub alphabet
						sub_alphabet.resize(d - k);
						sub_alphabet.assign(alphabet.begin(), alphabet.end() - k);

						do {
							// Subtract all lower combinations fval from
							linalg_less(cache.comFun [ alphabet ], cache.comFun [ sub_alphabet ]);
						} while (next_combination(alphabet.begin(), alphabet.end(), sub_alphabet.begin(), sub_alphabet.end()));
					}

					linalg_less(cache.comFun [ alphabet ], hdmrParam.fxBar);
				}
			}

			for (auto comFun : cache.comFun) {

				for (int i = 0; i < comFun.second.size(); ++i) {
					valSet[valShift + i] += comFun.second[i];
				}
			}
		}

		//Reset cache to zero
		cache.fval.clear();
		cache.comFun.clear();
	}

	MPI_Allreduce(MPI_IN_PLACE, valSet, pointCount * probParam.dof, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}








/*
* Index files
*/
void HDMR::writeIndexFile() {

	// Write index and log file
	if (computePool.grank == 0) {

		// Create index.data file
		fstream indexFile;
		indexFile.open(folderName + "index.data", ios::out | ios::trunc);

		indexFile << runTime.mode << "\n";
		indexFile << runTime.interpolationPoints << "\n";

		indexFile << probParam.dim << "\n";
		indexFile << probParam.dof << "\n";

		indexFile << sgParam.maxLevel << "\n";
		indexFile << sgParam.gridType << "\n";
		indexFile << sgParam.cutOff << "\n";

		if (runTime.mode == "HDMR_WRITE") {
			indexFile << hdmrParam.maxOrder << "\n";
			indexFile << hdmrParam.cutOff << "\n";

			indexFile << job.all_size << "\n";
			indexFile << job.active_size << "\n";

			for (int i = 0; i < hdmrParam.xBar.size(); ++i) {
				indexFile << hdmrParam.xBar[i] << "\t";
			}
			indexFile << "\n";

			for (int i = 0; i < hdmrParam.fxBar.size(); ++i) {
				indexFile << hdmrParam.fxBar[i] << "\t";
			}
			indexFile << "\n";

			// Write each job
			for (int i = 0; i < job.active.size(); ++i) {
				for (int j = 0; j < job.active[i].size(); ++j) {

					// Write order of dimension
					indexFile << i << "\t";
					for (int r = 0; r < job.active[i][j].size(); ++r) {
						indexFile << job.active[i][j][r] << "\t";
					}
					indexFile << "\n";
				}
			}
		}
		indexFile.close();

		// Create human readable log file
		fstream logFile;
		logFile.open(folderName + "log", ios::out | ios::trunc);
		logFile << "Run time mode : " << runTime.mode << "\n";
		logFile << "Interpolation points : " << runTime.interpolationPoints << "\n";
		logFile << "Average Write Time : " << runTime.avgTime << "\n";

		logFile << "DOF : " << probParam.dof << "\n";
		logFile << "Dim : " << probParam.dim << "\n";

		logFile << "SG Max Level : " << sgParam.maxLevel << "\n";
		logFile << "SG Grid Type : " << sgParam.gridType << "\n";
		logFile << "SG Cutoff : " << sgParam.cutOff << "\n";

		if (runTime.mode == "HDMR_WRITE") {
			logFile << "HDMR maxOrder : " << hdmrParam.maxOrder << "\n";
			logFile << "HDMR cutOff : " << hdmrParam.cutOff << "\n";

			logFile << "Total Job Size : " << job.all_size << "\n";
			logFile << "Total Active Job Size : " << job.active_size << "\n";

			logFile << "xBar : [ ";
			for (int i = 0; i < hdmrParam.xBar.size(); ++i) {
				logFile << hdmrParam.xBar[i] << " ";
			}
			logFile << "]\n";

			logFile << "fxBar : [ ";
			for (int i = 0; i < hdmrParam.fxBar.size(); ++i) {
				logFile << hdmrParam.fxBar[i] << " ";
			}
			logFile << "]\n";

			for (int d = 0; d < job.active.size(); ++d) {
				logFile << "HDMR Active Component Function [order , count] : " << d << " , " << job.active[d].size()  << "\n";
			}
		}
		logFile.close();
	}
}
void HDMR::readIndexFile() {

	// Load index file
	string lastMode = "";
	fstream indexFile;

	// Create index.data file
	indexFile.open(folderName + "index.data", ios::in);
	if (indexFile.is_open()) {
		indexFile >> lastMode;
		indexFile >> runTime.interpolationPoints;

		indexFile >> probParam.dim;
		indexFile >> probParam.dof;

		indexFile >> sgParam.maxLevel;
		indexFile >> sgParam.gridType;
		indexFile >> sgParam.cutOff;

		if (lastMode == "HDMR_WRITE") {
			indexFile >> hdmrParam.maxOrder ;
			indexFile >> hdmrParam.cutOff ;

			indexFile >> job.all_size ;
			indexFile >> job.active_size;


			// Load xbar and fxBar
			double valueD;
			hdmrParam.xBar.resize(probParam.dim);
			hdmrParam.fxBar.resize(probParam.dof);

			for (int i = 0; i < probParam.dim; ++i) {
				indexFile >> valueD;
				hdmrParam.xBar[i] = valueD;
			}

			for (int i = 0; i < probParam.dof; ++i) {
				indexFile >> valueD;
				hdmrParam.fxBar[i] = valueD;
			}

			// Load job list
			int value, d, j;
			job.active.resize(hdmrParam.maxOrder + 1);

			indexFile >> d;
			indexFile >> value;

			job.active[d].push_back(vector<int>());
			j = job.active[d].size() - 1;
			job.active[d][j].push_back(value);

			while ( indexFile ) {
				indexFile >> d;

				// Check character
				if (indexFile.peek() == EOF) {
					break;
				}

				job.active[d].push_back(vector<int>());
				j = job.active[d].size() - 1;

				for (int i = 0; i < d; ++i) {
					indexFile >> value;
					job.active[d][j].push_back(value);
				}
			}

			job.all = job.active;
		}
	}
	indexFile.close();

	// Set current runTime mode base on last mode
	if (lastMode == "SG_WRITE") {
		runTime.mode = "SG_READ";
	} else if (lastMode == "HDMR_WRITE")  {
		runTime.mode = "HDMR_READ";
	} else {
		err("ERROR invalid runtime mode in " + folderName);
	}
}









/*
* Job Allocation
*/
void HDMR::allocJobs() {

	// Initilize the job list
	job.all.resize(hdmrParam.maxOrder + 1);
	job.all[0].push_back(vector<int>());
	job.all[0][0].push_back(-1);
	job.all_size++;

	// Create base alphabet
	vector<int> alphabet(probParam.dim);
	for (int i = 0; i < alphabet.size(); ++i) {
		alphabet[i] = i;
	}
	vector<int> sub_alphabet;

	// Generate all combinations of activeDim and save them to job
	for (int i = 1; i <= hdmrParam.maxOrder; ++i) {

		// Initilize step counter and create base of sub alphabet
		int n = 0;
		sub_alphabet.resize(i);
		sub_alphabet.assign(alphabet.begin(), alphabet.begin() + i);
		do {
			// Generate a matrix of vectors storing the jobs
			job.all[i].push_back(vector<int> (i));
			job.all[i][n].assign(sub_alphabet.begin(), sub_alphabet.end());
			n++;
			job.all_size++;
		} while (next_combination(alphabet.begin(), alphabet.end(), sub_alphabet.begin(), sub_alphabet.end()));
	}
}
void HDMR::allocActiveJobs(int method) {

	// Copy all jobs to active jobs
	if (method == 0) {
		job.active = job.all;
		job.active_size = job.all_size;
	} else if (method == 1) {

	}
}








/*
* Allocaet Cache
*/
void HDMR::allocCache() {
	cache.fval.clear();
	cache.comFun.clear();
}

void HDMR::allocCacheFval(double* x) {

	int activeJobIndex;
	vector<double> xPartial;
	vector<double> fvalPartial(probParam.dof, 0.0);

	// Initilize indexing and offset variables
	activeJobIndex = 0;

	// Cycle through each order of the active job (zero order has already beed set to fval)
	for (int d = 1; d <= hdmrParam.maxOrder ; ++d) {

		// Resize lower order input x
		xPartial.resize(d);

		// Cycle through job.active
		for (int i = 0; i < job.active[d].size(); ++i) {

			// Allocate lower order x based on the active index from active job list
			for (int j = 0; j < d; ++j) {
				xPartial[j] = x[job.active[d][i][j]];
			}
			sgread[activeJobIndex]->interpolateValue(&xPartial[0], &fvalPartial[0]);
			cache.fval[job.active[d][i]] = fvalPartial;
			activeJobIndex ++;
		}
	}
}
void HDMR::allocCacheFval_prl(double* x) {

	/*
		int activeJobIndex, cacheShift;
		vector<double> xPartial;

		//Copy fbar into cache only for rank zeros (so we dont double count)
		if (computePool.grank == 0) {
			copy(hdmrParam.fxBar.begin(), hdmrParam.fxBar.end(), cache.fval);
		}

		// Initilize indexing and offset variables
		activeJobIndex = 0;
		cacheShift = probParam.dof;

		// Cycle through each order of the active job (zero order has already beed set to fval)
		for (int d = 1; d <= hdmrParam.maxOrder ; ++d) {

			// Resize lower order input x
			xPartial.resize(d);

			// Cycle through job.active
			for (int i = 0; i < job.active[d].size(); ++i) {

				// round robin style task allocation
				if (i % computePool.nodeSize ==  computePool.nodeRank) {

					// Allocate lower order x based on the active index from active job list
					for (int j = 0; j < d; ++j) {
						xPartial[j] = x[job.active[d][i][j]];
					}

					// Interpolate the value and save it to cache
					sgread[activeJobIndex]->interpolateValue(&xPartial[0], &cache.fval[cacheShift]);

					// Only keep the local zero rank so we dont double count int MPI reduce
					if (computePool.rank != 0) {
						fill(&cache.fval[cacheShift], &cache.fval[cacheShift] + probParam.dof, 0.0);
					}
				}

				// Increment indexing and offset variable
				cacheShift += probParam.dof;
				activeJobIndex ++;
			}
		}

		// Reduce the parital results of fval cache place all, thus all nodes have the same thing
		MPI_Allreduce(MPI_IN_PLACE, &cache.fval[0], job.active_size * probParam.dof, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


			cout << endl;

			for (int i = 0; i < probParam.dim; ++i) {
				cout << x[i] << " ";
			}
			cout << endl;
			cout << "~~~~~~~~~~~~~~~~~~~" << endl;
			for (int i = 0; i < job.active_size; ++i) {
				for (int j = 0; j < probParam.dof; ++j) {
					cout << cache.fval[i * probParam.dof + j] << " ";
				}
				cout << endl;
			}
			cout << computePool.grank << " -> " << cache.fval[13] << endl;
			cout << hline;
		*/


}







/*
* Compute Allocation
*/
void HDMR::allocCompute() {

	if (runTime.mode == "SG_WRITE") {
		// 1 node with all process in that node

		computePool.subComm = MPI_COMM_WORLD;

		MPI_Comm_size(MPI_COMM_WORLD, &computePool.gsize);
		MPI_Comm_rank(MPI_COMM_WORLD, &computePool.grank);

		computePool.size = computePool.gsize;
		computePool.rank = computePool.grank;

		computePool.nodeRank = 0;
		computePool.nodeSize = 1;
		computePool.processPerNode = computePool.gsize;

	} else if (runTime.mode == "SG_READ") {
		// Mutiple nodes with processPerNode=1

		MPI_Comm_size(MPI_COMM_WORLD, &computePool.gsize);
		MPI_Comm_rank(MPI_COMM_WORLD, &computePool.grank);

		computePool.processPerNode = 1;
		computePool.nodeSize = computePool.gsize / computePool.processPerNode ;

		computePool.nodeRank = computePool.grank % computePool.nodeSize;

		MPI_Comm_split(MPI_COMM_WORLD , computePool.nodeRank , computePool.grank , &computePool.subComm );

		MPI_Comm_size(computePool.subComm , &computePool.size);
		MPI_Comm_rank(computePool.subComm , &computePool.rank);

	} else if (runTime.mode == "HDMR_WRITE" ) {
		// Multiple nodes with processPerNode defined by intput

		MPI_Comm_size(MPI_COMM_WORLD, &computePool.gsize);
		MPI_Comm_rank(MPI_COMM_WORLD, &computePool.grank);

		computePool.nodeSize = computePool.gsize / computePool.processPerNode ;

		if (computePool.nodeSize < 1) {err("ERROR Invalid nodeSize, processPerNode cannot be > global process size");}

		computePool.nodeRank = computePool.grank % computePool.nodeSize;

		MPI_Comm_split(MPI_COMM_WORLD , computePool.nodeRank , computePool.grank , &computePool.subComm );

		MPI_Comm_size(computePool.subComm , &computePool.size);
		MPI_Comm_rank(computePool.subComm , &computePool.rank);

	} else if (runTime.mode == "HDMR_READ") {

		// Mutiple nodes with processPerNode=1

		MPI_Comm_size(MPI_COMM_WORLD, &computePool.gsize);
		MPI_Comm_rank(MPI_COMM_WORLD, &computePool.grank);

		computePool.processPerNode = 1;
		computePool.nodeSize = computePool.gsize / computePool.processPerNode ;

		computePool.nodeRank = computePool.grank % computePool.nodeSize;

		MPI_Comm_split(MPI_COMM_WORLD , computePool.nodeRank , computePool.grank , &computePool.subComm );

		MPI_Comm_size(computePool.subComm , &computePool.size);
		MPI_Comm_rank(computePool.subComm , &computePool.rank);

	} else {
		err("ERROR In allocCompute(), runTime.mode not defined");
	}

	MPI_Barrier(MPI_COMM_WORLD);
}







/*
* Timer
*/
void HDMR::startTimer() {
	runTime.time = MPI_Wtime() * -1.0;
}
void HDMR::endTimer() {

	// Calculate time & average run time
	runTime.time += MPI_Wtime();
	runTime.avgTime = runTime.time ;
	runTime.maxTime = runTime.time ;
	runTime.minTime = runTime.time ;

	MPI_Allreduce(MPI_IN_PLACE, &runTime.avgTime, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	runTime.avgTime = runTime.avgTime / computePool.gsize;

	MPI_Allreduce(MPI_IN_PLACE, &runTime.maxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &runTime.minTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
}







/*
* Debug
*/
void HDMR::debug( string heading, int showRunTime, int showComputePool, int showProb, int showSG_HDMRparam, int showJob) {

	if (computePool.grank == 0) {cout << hline;};

	MPI_Barrier(MPI_COMM_WORLD);
	sleep( (computePool.nodeRank * computePool.nodeSize + computePool.rank) / computePool.gsize);

	cout << "\n " + heading + " - State Dump [ Global Rank: " << computePool.grank << " of " << computePool.gsize << " ]\n";
	cout << "================================================\n";

	if (showRunTime) {
		if (showRunTime == 1) {
			if (computePool.grank == 0) {
				cout << "\n> RunTime Param\n";
				cout << "  System Mode: " << runTime.mode << endl;
				cout << "  Verbose: " << runTime.verbose << endl;
				cout << "  Total Process Avg Time: " << runTime.avgTime << endl;
				cout << "  Total Process Max Time: " << runTime.maxTime << endl;
				cout << "  Total Process Min Time: " << runTime.minTime << endl;
				cout << "  Total Interpolation Points: " << runTime.interpolationPoints << endl;
			}
		} else {
			cout << "\n> RunTime Param\n";
			cout << "  System Mode: " << runTime.mode << endl;
			cout << "  Verbose: " << runTime.verbose << endl;
			cout << "  Process Time: " << runTime.time << endl;
			cout << "  Total Process Avg Time: " << runTime.avgTime << endl;
			cout << "  Total Process Max Time: " << runTime.maxTime << endl;
			cout << "  Total Process Min Time: " << runTime.minTime << endl;
			cout << "  Total Interpolation Points: " << runTime.interpolationPoints << endl;
		}
	}


	if (showComputePool) {
		cout << "\n> Compute Pool\n";
		cout << "  [group][local]: "
		     << "[" << computePool.nodeRank << " of " << computePool.nodeSize << "] "
		     << "[" << computePool.rank      << " of " << computePool.size      << "]" << endl;
		cout << "  processPerNode: " << computePool.processPerNode << endl;
	}

	if (showProb && computePool.grank == 0) {
		cout << "\n> Problem Param\n";
		cout << "  Dim: " << probParam.dim << endl;
		cout << "  DOF: " << probParam.dof << endl;
	}

	if (showSG_HDMRparam && computePool.grank == 0) {
		cout << "\n> SG param\n";
		cout << "  Lmax: " << sgParam.maxLevel << endl;
		cout << "  Type: " << sgParam.gridType << endl;
		cout << "  CutOff/Epsilon: " << sgParam.cutOff << endl;


		if (runTime.mode == "HDMR_WRITE" || runTime.mode == "HDMR_READ") {
			cout << "\n> HDMR param\n";
			cout << "  maxOrder: " << hdmrParam.maxOrder << endl;
			cout << "  CutOff/Epsilon: " << hdmrParam.cutOff << endl;

			cout << "\n> Xbar & F(Xbar)\n";
			cout << "  Xbar : [ ";

			for (int i = 0; i < hdmrParam.xBar.size(); ++i) {
				cout << hdmrParam.xBar[i] << " ";
			}
			cout << "]" << endl;

			cout << "  f(Xbar) : [ ";
			for (int i = 0; i < hdmrParam.fxBar.size(); ++i) {
				cout << hdmrParam.fxBar[i] << " ";
			}
			cout << "]" << endl;
		}
	}

	if (showJob && computePool.grank == 0) {

		cout << "\n> Job List\n";
		if (showJob > 1) {
			cout << "  Full List:\n";
			for (int i = 0; i < job.all.size(); ++i) {
				cout << "  [" << i << "][" << job.all[i].size() << "] ";
				for (int j = 0; j < job.all[i].size(); ++j) {
					for (int r = 0; r < job.all[i][j].size(); ++r) {
						cout << job.all[i][j][r];
					}
					cout << " ";
				}
				cout << endl;
			}
			cout << endl;
		}

		cout << "  Active List:\n";
		for (int i = 0; i < job.active.size(); ++i) {
			cout << "  [" << i << "][" << job.active[i].size() << "] ";
			for (int j = 0; j < job.active[i].size(); ++j) {
				for (int r = 0; r < job.active[i][j].size(); ++r) {
					cout << job.active[i][j][r];
				}
				cout << " ";
			}
			cout << endl;
		}
	}


	MPI_Barrier(MPI_COMM_WORLD);
	if (computePool.grank == 0) {sleep(1); cout << hline;};
}


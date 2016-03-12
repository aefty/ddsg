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
#include "cpp_util/include.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <map>
#include <algorithm>
#include <unistd.h> // Sleep

using std::cout;
using std::string;
using std::vector;
using std::map;
using std::find;
using std::fill;
using namespace cpp_util;

/*
 * Constructur, Destructor and clear
*/
HDMR::HDMR(int verbose_) {
	runTime.verbose = verbose_;
}
HDMR::~HDMR() {

	for (int i = 0; i < job.active_size; ++i) {
		delete sgread[i];
	}
	delete[] sgread;

	clear();
}
void HDMR::clear() {

	MPI_Barrier(MPI_COMM_WORLD);

	folderName = {};
	l2convergence = {};

	probParam = {};
	sgParam = {};
	hdmrParam = {};
	computePool = {};
	job = {};
	runTime = {};

	resetCache();
}
void HDMR::resetCache() {
	cache.fval.clear();
	cache.fcomFun.clear();
	cache.ival.clear();
	cache.icomFun.clear();
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

	// Check for consistency between maxOrder and probelm Dim
	if (hmdrMaxOrder_ > problemDim_) {
		err("ERROR HDMR MaxOrder cannot be > Problem Dim");
	}

	// Clear all variables
	clear();

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
	setComputePool();

	// Create working directory
	if (computePool.grank == 0) {
		string command = "";

		command = "rm -rf " + folderName + " > nul 2>&1"; // Supress ouput
		system(command.c_str());

		command = "mkdir " + folderName + " > nul 2>&1"; // Supress ouput
		system(command.c_str());
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// Allocate full job on all processes
	setAllJobs();


	sgwrite = new SGwrite*[hdmrParam.maxOrder];
	for (int d = 1; d <= hdmrParam.maxOrder; ++d) {
		sgwrite[d - 1] = new SGwrite(problemFunc, d, probParam.dof , sgParam.maxLevel, sgParam.cutOff, sgParam.gridType , runTime.verbose);
		sgwrite[d - 1]->setAnchor(hdmrParam.xBar);
	}

	// Generate the active jobs to be processed
	if (hdmrParam.cutOff == 0.0) {
		runTime.interpolationPoints = setActiveJobs_noAdaptivity();
	} else {
		runTime.interpolationPoints = setActiveJobs_integralAdaptivty();
	}

	// Write index and log file
	writeIndexFile();


	// Delete all instance of sgwrite
	MPI_Barrier(MPI_COMM_WORLD);

	// Delete instance of sgwrite
	for (int d = 1; d <= hdmrParam.maxOrder; ++d) {
		delete sgwrite[d - 1];
	}
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
	setComputePool();

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
	sgwrite[0]->build();
	runTime.interpolationPoints = sgwrite[0]->write(folderName + "surplus.data");

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
	setComputePool();

	if (runTime.mode == "SG_READ") {

		// Initilize SGread instance on all processess
		sgread = new SGread*[1];
		sgread[0] = new SGread(runTime.verbose);
		sgread[0]->resetMPI(computePool.subComm);
		runTime.interpolationPoints = sgread[0]->read(folderName + "surplus.data");
	} else {


		// NOTES: ALL active files need to be loaded on all process.
		// This is required as each procss will later work indpenendalty
		// to interpolate each given point. Each interpolation call will
		// access all surplus files. As such we need each file and associated
		// SGread object to be avaiable on all process.


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

				// NO prallelization here ! No Round robin task allocation
				runTime.interpolationPoints += sgread[k]->read(folderName + vector_join(job.active[d][i], ".") + ".data");
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

	resetCache();

	// Cycle through each point
	for (int r = 0; r < pointCount; ++r) {

		if (r % computePool.nodeSize ==  computePool.nodeRank) {

			// Set shifting constants
			xShift = r * probParam.dim;
			valShift = r * probParam.dof;

			// Add the final -f_0 for at the start;
			linalg_add(&valSet[valShift], hdmrParam.fxBar);

			// Allocate the cache.fval on current process
			setCacheFval(&xSet[xShift]);

			// *** Kernel of the code!
			// Calculate Component functions

			// Cycel through jobs
			for (int d = 1; d < job.active.size(); ++d) {

				alphabet.resize(d);
				for (int i = 0; i < job.active[d].size(); ++i) {

					// Initilize alphabit at the activeJob and set the component function intial value
					alphabet = job.active[d][i] ;
					cache.fcomFun [ alphabet ] = cache.fval[job.active[d][i]];

					// Generate all combinations of lower level combination of alphabit
					for (int k = 1; k < d; ++k) {

						// Intilize sub alphabet
						sub_alphabet.resize(d - k);
						sub_alphabet.assign(alphabet.begin(), alphabet.end() - k);

						do {
							// Subtract all lower combinations fval from
							linalg_less(cache.fcomFun [ alphabet ], cache.fcomFun [ sub_alphabet ]);
						} while (next_combination(alphabet.begin(), alphabet.end(), sub_alphabet.begin(), sub_alphabet.end()));
					}

					linalg_less(cache.fcomFun [ alphabet ], hdmrParam.fxBar);
				}
			}

			for (auto fcomFun : cache.fcomFun) {
				for (int i = 0; i < fcomFun.second.size(); ++i) {
					valSet[valShift + i] += fcomFun.second[i];
				}
			}
		}

		//Reset cache to zero
		cache.fval.clear();
		cache.fcomFun.clear();
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
void HDMR::setAllJobs() {

	// Initilize the job list
	job.all.resize(hdmrParam.maxOrder + 1);
	job.all[0].push_back(vector<int>());
	job.all[0][0].push_back(-1);
	job.all_size++;

	//job.full.resize(hdmrParam.maxOrder + 1);
	//job.full[0].push_back(list<vector<int>>);

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
int HDMR::setActiveJobs_noAdaptivity() {

	// Generate lower dimension surplus files for all active Jobs
	vector<int> activeDim;
	string fileName = "";
	int interpolationPoints = 0;
	int jobIndex = 1; // Ignore job 0, that just the constant

	// Set active job to be all the jobs
	job.active = job.all;
	job.active_size = job.all_size;

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

				// Set computPool
				sgwrite[d - 1]->resetMPI(computePool.subComm);

				//Build and write sparse grid
				sgwrite[d - 1]->build(activeDim);
				interpolationPoints += sgwrite[d - 1]->write(fileName);
				sgwrite[d - 1]->clear();
			}
		}
	}

	// Sum up all interpolation points accross all process
	if (computePool.rank != 0) {
		interpolationPoints = 0;
	}
	MPI_Allreduce(MPI_IN_PLACE, &interpolationPoints, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

	l2convergence = -1;

	return interpolationPoints;
}
int HDMR::setActiveJobs_integralAdaptivty() {

	// Generate lower dimension surplus files for all active Jobs
	vector<int> activeDim;
	vector<int> sub_activeDim;

	vector<double> sum_icomfun(probParam.dof, 0.0);

	vector<double> integral(probParam.dof, 0.0);
	vector<double> integral_last(probParam.dof, 0.0);

	vector <vector<int>> cadidateDim;
	vector <int> rejectDimIndex;

	int interpolationPoints = 0;
	double nu = 100;
	double rho = 100;
	double rho_last = 100;

	vector<double> ival;
	vector<double> mpiContainer_icomFun;

	string fileName = "";

	// Set active job to all job
	job.active = job.all;

	integral_last = hdmrParam.fxBar;

	for (int d = 1; d <= hdmrParam.maxOrder; ++d) {

		// Check if we have any jobs left
		if (job.active[d].size() < 1) {
			break;
		}

		// Check if we have converged
		if (rho < hdmrParam.cutOff) {

			for (int i = d; i <= hdmrParam.maxOrder; ++i) {
				job.active[i].clear();
			}
			break;
		}

		// Set candindate dimenension set and container for active dimension
		cadidateDim = job.active[d];

		activeDim.clear();
		activeDim.resize(d);

		rejectDimIndex.clear();
		rejectDimIndex.resize(cadidateDim.size(), 0);

		// Clear and resize MPI container
		ival.clear();
		ival.resize(cadidateDim.size()*probParam.dof, 0.0);

		mpiContainer_icomFun.clear();
		mpiContainer_icomFun.resize(cadidateDim.size()*probParam.dof, 0.0);

		// Parralel processing of job list
		for (int i = 0; i < cadidateDim.size(); ++i) {
			if (i % computePool.nodeSize ==  computePool.nodeRank) {

				// Get active dimensions
				activeDim = cadidateDim[i];
				fileName = folderName + vector_join(activeDim, ".") + ".data";

				// Instatiate the sparse grid
				sgwrite[d - 1]->resetMPI(computePool.subComm);
				sgwrite[d - 1]->build(activeDim);
				sgwrite[d - 1]->integrateDomain(&ival[i * probParam.dof]);

				// Set the relvant component of mpiContainer_icomFun to ival
				// icomFun = ival -  sum_{ all permutations of =u }(icomFun_u)
				// Initilize icomFun = ival ... for relevant components in computepool
				copy(&ival[i * probParam.dof], &ival[i * probParam.dof] + probParam.dof  , &mpiContainer_icomFun[i * probParam.dof]);

				//Intilize the sum_icomfun vector to be zero
				fill(sum_icomfun.begin(), sum_icomfun.end(), 0.0);

				// Calculate : sum_{ all permutations of =u }(icomFun_u)
				if (d > 1) {

					for (int k = 1; k < d; ++k) {
						// Genearate sub_activeDim
						sub_activeDim.resize(d - k);
						sub_activeDim.assign(activeDim.begin(), activeDim.end() - k);

						// Sum up all lower combinations of the icomFun
						// Calculate : sum_{ all permutations of =u }(icomFun_u)
						do {
							linalg_add(sum_icomfun, cache.icomFun[sub_activeDim]);
						} while (next_combination(activeDim.begin(), activeDim.end(), sub_activeDim.begin(), sub_activeDim.end()));
					}
				}

				// subtract the zeroth order icomFun (which is integral of fxBar =  fxbar)
				linalg_add(sum_icomfun, hdmrParam.fxBar);

				// Calculate  = ival - sum_{ all permutations of =u }(icomFun_u)
				linalg_less(&mpiContainer_icomFun[i * probParam.dof], sum_icomfun);

				nu = linalg_l2(&mpiContainer_icomFun[i * probParam.dof], probParam.dof) / linalg_l2(sum_icomfun);

				printd(activeDim, "activeDim");
				printd(&ival[i * probParam.dof], probParam.dof, "ival");
				printd(sum_icomfun, "sum_icomfun");
				printd(&mpiContainer_icomFun[i * probParam.dof], probParam.dof, "icomFun");

				printd(nu, "Nu");

				// Add dimension index to reject list
				if (nu < hdmrParam.cutOff) {
					rejectDimIndex[i] = 1;
					printd("Reject Dim !!!!!!!!!!");
				}

				// Write out surplus file only if its not a rejected dimension or if its first order
				if (!rejectDimIndex[i] || d == 1 ) {
					interpolationPoints += sgwrite[d - 1]->write(fileName);
					printd("Accept Dim");
				}

				// clear object
				sgwrite[d - 1]->clear();
			}
		}

		// Clear the contents of the local rank 0, so we dont double count in reduceAll and Share the the integral globally
		if (computePool.rank != 0) {
			fill(mpiContainer_icomFun.begin(), mpiContainer_icomFun.end() , 0.0);
		}

		MPI_Allreduce(MPI_IN_PLACE, &mpiContainer_icomFun[0], mpiContainer_icomFun.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &rejectDimIndex[0], rejectDimIndex.size(), MPI_INT, MPI_MAX, MPI_COMM_WORLD);

		// Allocate the icomFun to all local cache
		for (int i = 0; i < cadidateDim.size(); ++i) {
			cache.icomFun[cadidateDim[i]].assign(&mpiContainer_icomFun[i * probParam.dof], &mpiContainer_icomFun[i * probParam.dof] + probParam.dof);
		}


		printd(rejectDimIndex, "rejectDimIndex");
		printd();

		// Remove hdmr jobs that are not required
		for (int i = 0; i < cadidateDim.size(); ++i) {

			// Set alphabet for cobination
			activeDim = cadidateDim[i];

			if (rejectDimIndex[i]) {
				if (d == 1) {
					removeDimFromJobs(activeDim, 2);
				} else {
					removeDimFromJobs(activeDim, d);
				}
			}
		}

		// Initilize temp as zero and sum up all components
		integral = hdmrParam.fxBar;
		for (auto icomFun : cache.icomFun) {
			linalg_add(integral, icomFun.second);
		}

		if (computePool.grank == 0) {
			printd(integral_last, "integral_last");
			printd(integral, "integral");
		}

		rho = 1.0 / linalg_l2(integral_last);
		linalg_less(integral_last, integral);

		printd(integral_last, "-Diff");

		rho *= linalg_l2(integral_last);

		if (computePool.grank == 0) {
			printd(rho, "RHO");
			printd();
		}

		integral_last = integral;
	}

	l2convergence = rho;

	// Update job.active_size and new hdmrParam.maxOrder
	for (int d = 0; d < job.active.size(); ++d) {
		job.active_size += job.active[d].size();
		if (job.active[d].size() != 0) {
			hdmrParam.maxOrder = d;
		}
	}

	// Sum up all interpolation points accross all process
	if (computePool.rank != 0) {
		interpolationPoints = 0;
	}
	MPI_Allreduce(MPI_IN_PLACE, &interpolationPoints, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

	return interpolationPoints;
}



void HDMR::removeDimFromJobs(vector<int>& activeDim, int d_start) {

	int match;
	vector<vector<int>> temp;

	for (int d = d_start; d < job.active.size(); ++d) {

		temp.clear();

		for (int i = 0; i < job.active[d].size(); ++i) {
			match = 0;
			for (int k = 0; k < activeDim.size(); ++k) {
				match += contains(job.active[d][i], activeDim[k]);
			}

			if (match != activeDim.size() ) {
				printd(job.active[d][i], "ACCEPT THIS");
				temp.push_back(job.active[d][i]);
			}
		}
		job.active[d].clear();
		job.active[d] = temp;

		printd(job.active[d], "job.active[d]");
	}
}




void HDMR::setCacheFval(double* x) {

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



/*
* Compute Allocation
*/
void HDMR::setComputePool() {

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
		err("ERROR In setComputePool(), runTime.mode not defined");
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
void HDMR::debug( string heading, int rankToShow, int showRunTime, int showComputePool, int showProb, int showSG_HDMRparam, int showJob) {

	if (computePool.grank == 0) {cout << hline;};

	MPI_Barrier(MPI_COMM_WORLD);

	if (rankToShow == -1) {
		sleep( (computePool.nodeRank * computePool.nodeSize + computePool.rank) / computePool.gsize);
	} else {
		if (computePool.grank != rankToShow) {
			return;
		}
	}

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
				cout << "  HDMR L2 Convergence: " << l2convergence << endl;
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
			cout << "  HDMR L2 Convergence: " << l2convergence << endl;
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
		cout << "  Grid Type: " << sgParam.gridType << endl;
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

	if (computePool.grank == 0) {sleep(1); cout << hline;};
}


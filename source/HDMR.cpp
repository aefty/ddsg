/*
High Dimensional Model Representation
======================================
This class is intended to be used in conjunction with the AdaptiveSparseGrid library of "Xiang Ma, xm25@cornell.edu".
The interface of the two class are through SGwrite & SGread.
Status Note: Work in progress.
aryan.eftekhari@usi.ch/gmail.com
*/

#include "HDMR.h"
#include "SGwrite.h"
#include "SGread.h"

#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <random>
#include <string>
#include <map>
#include <algorithm>
#include <unistd.h> // Sleep


using namespace std;

// General

HDMR::HDMR(int verbose) {
	runTime.verbose = verbose;
}

HDMR::~HDMR() {
	cleanUp();
}

void HDMR::cleanUp() {

	MPI_Barrier(MPI_COMM_WORLD);
	job = {};
	data = {};

	computePool = {};
	probParam = {};
	hdmrParam = {};
	sgParam = {};

	problemFunc = NULL;

	runTime = {};

	localCache.taskLookUpIndex.clear();
	localCache.interpolate.clear();

	localCache.blackList.clear();
	localCache.integral.clear();

	if (runTime.writeState == 1) {
		for (int k = 0; k < hdmrParam.maxOrder; ++k) {
			delete sgwrite[k];
		}
	} else {
		int current_k = computePool.blockID;
		int workload = computePool.workload[current_k];

		for (int i = 0; i < workload; ++i) {
			delete sgread[i];
		}
	}
}

void HDMR::computeAlloc(string computeArch, int optimalThreadsPerBlock ) {

	//Initialize MPI parameters if not already done
	if (computePool.size == -1) {
		MPI_Comm_rank(MPI_COMM_WORLD, &computePool.rank);
		MPI_Comm_size(MPI_COMM_WORLD, &computePool.size);
	}

	//Setup computePool architecture
	if (computeArch == "flat") { //Compute resource mapped to optimalThreadsPerBlock

		// If mpi Size is greater or "equal" 2 then use different setting than default
		if (computePool.size >= optimalThreadsPerBlock) {
			computePool.opt_tpb = optimalThreadsPerBlock; // Threads per Chip

			//Split Communicator into N groups.
			computePool.blockID = computePool.rank / computePool.opt_tpb;
			computePool.blockSize = int (ceil( double (computePool.size) / double(computePool.opt_tpb)));

			MPI_Comm_split(MPI_COMM_WORLD , computePool.blockID , computePool.rank , &computePool.subComm );
			MPI_Comm_rank(computePool.subComm , &computePool.threadID );
			MPI_Comm_size(computePool.subComm , &computePool.tpb );
		}
	} else if (computeArch == "bigPool") {

		// Allocate one single block with all ranks equalt to a thread. d
		computePool.blockID = 0;
		computePool.blockSize = 1;

		computePool.threadID = computePool.rank;
		computePool.tpb = computePool.size;
		computePool.opt_tpb = computePool.size;

		computePool.workload.clear();
		computePool.workload.push_back(0);
		for (int k = 0; k < hdmrParam.maxOrder; ++k) {
			computePool.workload[0] += job[k].size;
		}

		computePool.subComm = MPI_COMM_WORLD;
	}
}

//Write Routines
inline int HDMR::traverseIntegral(vector<int> alphabet, int k, double* fvalue, int op) {

	int comboCount = nCk(alphabet.size(), k);
	int comboBuffer[comboCount * k];
	int ignore = 0;
	string currentIndex = "";
	vector<int> subAlphabet;

	// Degrade alphabet into total of n=comboCount and of length k
	combination(alphabet, k, comboBuffer);

	// Cycle through each combination set
	for (int i = 0; i < comboCount; ++i) {

		// Make surplus file name & construct the correct x value
		currentIndex = "";
		for (int r = 0; r < k - 1; ++r) {
			currentIndex += to_string(comboBuffer[i * k + r]) + ".";
		}
		currentIndex += to_string(comboBuffer[i * k + k - 1]);

		// If black listed ignore by return 1
		if (std::find(localCache.blackList.begin(), localCache.blackList.end(), currentIndex) != localCache.blackList.end()) {
			return 1;
		}

		// Apply the operation op on the new interpolation
		for (int j = 0; j < probParam.dof; ++j) {
			fvalue[j] +=  localCache.integral[currentIndex][j] * op ;
		}

		// cout << "	" << currentIndex << " " << localCache.integral[currentIndex][0] << "(" << op << ")" << endl;

		// Construct a new vector holding the sub alphabet of the above combination
		subAlphabet.clear();
		for (int j = 0; j < k; ++j) {
			subAlphabet.push_back(comboBuffer[i * k + j]);
		}

		// Recrusive call to traverse tree for each K level. If ignore=1, return 1 ... ending the recursion
		for (int r = k - 1; r > 0; --r) {
			ignore = traverseIntegral(subAlphabet, r, fvalue, (op * -1));
			if (ignore == 1) { return 1; }
		}

		// Deduct fxBar from fValeue
		for (int j = 0; j < probParam.dof ; ++j) {
			fvalue[j]  += (op * -1.0)  * localCache.integral["0"][j];
		}

		// cout << "	" << "f0" << " " <<  localCache.integral["0"][0]  << "(" <<  -1 * op << ")" << endl;
	}

	return 0;
}

void HDMR::writeJobAlloc() {

	double f1value[probParam.dof];
	double f0value[probParam.dof];
	double deltaf[probParam.dof];
	double rho = 0;
	double sumrho = 0;
	int ignore = 0;
	string currentIndex;
	vector<double> fullPointMap(probParam.dim);
	vector<int> subAlphabet;
	vector<int> alphabet(probParam.dim);

	if (computePool.rank == 0 && runTime.verbose ) {cout << endl;}

	// Initlize stuff
	for (int i = 0; i < probParam.dof; ++i) {
		localCache.integral["0"].push_back(data.fxBar[i]);
	}

	for (int i = 0; i < alphabet.size(); ++i) {
		alphabet[i] = i + 1;
	}

	//Create and array of SGwrite classes for each HDMR order (on every rank)
	sgwrite = new SGwrite*[hdmrParam.maxOrder];

	// Create array of jobs for each HDMR order
	job = new Job[hdmrParam.maxOrder];
	for (int k = 0; k < hdmrParam.maxOrder; ++k) {
		sgwrite[k] = new SGwrite(problemFunc, k + 1 , probParam.dof , sgParam.maxLevel, sgParam.cutOff, sgParam.gridType , 0);
		job[k].size = 0;
	}

	for (int k = 1; k < hdmrParam.maxOrder + 1; ++k) {

		int comboCount =  nCk(alphabet.size(), k) ;
		int comboBuffer[ comboCount * k ];

		// Allocate the job properties
		job[k - 1].task = new int[k * comboCount]();

		// Generate combination of each different tasks
		combination(alphabet, k, comboBuffer);

		// Cycle through each parital funciton
		for (int i = 0; i < comboCount; ++i) {

			// Define current index string
			currentIndex = "";
			for (int j = 0; j < k - 1; ++j) {
				currentIndex += to_string(comboBuffer[i * k + j]) + ".";
			}
			currentIndex += to_string(comboBuffer[i * k + k - 1]);

			//Copy xBar into it full map
			//Select the combination index from the task list and change fullPointMap and also make fileName
			for (int j = 0; j < probParam.dim; ++j) {
				fullPointMap[j] = data.xBar[j];
			}
			for (int j = 0; j < k ; ++j) {
				fullPointMap[ comboBuffer[i * k + j] - 1] = -1;
			}
			for (int j = 0; j < probParam.dof; ++j) {
				f1value[j] = 0;
				f0value[j] = data.fxBar[j];
			}

			//Build sparse grid & integrate
			sgwrite[k - 1]->build(fullPointMap);
			sgwrite[k - 1]->integrateDomain(f1value);
			sgwrite[k - 1]->Cleanup();

			// add it to the cache
			for (int j = 0; j < probParam.dof; ++j) {
				localCache.integral[currentIndex].push_back(f1value[j]);
			}

			// cout <<  currentIndex << "->" << f1value[0] << endl;

			// Write sub alphabet
			subAlphabet.clear();
			for (int j = 0; j < k ; ++j) {
				subAlphabet.push_back( comboBuffer[i * k + j]);
			}

			for (int subk = k - 1; subk > 0; --subk) {
				ignore = traverseIntegral(subAlphabet, subk , f0value, 1);
				if (ignore == 1) {break;}
			}

			if (ignore != 1) {
				// calculate difference
				for (int j = 0; j < probParam.dof; ++j) {
					deltaf[j] = f1value[j] - f0value[j];
				}

				rho = l2norm(deltaf, probParam.dof) / l2norm(f0value, probParam.dof);
				if (rho < hdmrParam.cutOff) { localCache.blackList.push_back(currentIndex); }

				sumrho += rho;

				bool blackListCheck = std::find(localCache.blackList.begin(), localCache.blackList.end(), currentIndex) == localCache.blackList.end();

				if (blackListCheck || k == 1) {

					for (int j = 0; j < k; ++j) {
						job[k - 1].task[ (job[k - 1].size) * k + j ] =  comboBuffer[i * k + j];
					}

					job[k - 1].size++;
				}
			}
		}

		// VERBOSE OUTPUT
		if (computePool.rank == 0 && runTime.verbose ) {
			if (sumrho == 0 ||  job[k - 1].size == 0) {
				cout << "[" << k << "] Order Terminated: [" << sumrho << "/" << job[k - 1].size << "] {" << job[k - 1].size << "/" << comboCount << "}" << endl;
			} else {
				cout << "[" << k << "] Order Rho: " << sumrho / job[k - 1].size << " {" << job[k - 1].size << "/" << comboCount << "}" << endl;
			}
		}

		// Break if the the last elvel
		if (job[k - 1].size == 0) { break; }
		sumrho = 0;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	return;
}

int HDMR::write(
    void (*problemFunc_)(double*, int, double*), int problemDim_, int problemDoF_ ,
    int maxOrder_ , int sampleSize_ , double hdmrCutOff_ ,
    int maxLevel_ , int gridType_ , double sgCutOff_ ,
    string folderName , vector<double>preDefinedXbar) {


	// Initialize Global Timer - Start
	runTime.writeState = 1;
	runTime.start = MPI_Wtime();

	//Initialize Problem Parameters
	problemFunc          = problemFunc_;
	probParam.dim        = problemDim_;
	probParam.dof        = problemDoF_;

	//Initialize HDMR parameters
	hdmrParam.maxOrder   = maxOrder_;
	hdmrParam.sampleSize = sampleSize_;
	hdmrParam.cutOff     = hdmrCutOff_;

	//Initialize Sparse grid parameters
	sgParam.maxLevel     = maxLevel_;
	sgParam.gridType     = gridType_;
	sgParam.cutOff       = sgCutOff_;

	//Initialize xBar and f(xBar) array
	data.xBar = new double[probParam.dim]();
	data.fxBar = new double[probParam.dof]();


	//Initialize process allocation default with optimalThreadsPerBlock_
	if (computePool.size == -1) {
		computeAlloc("flat", 2);
	}

	// HDMR order cant be bigger or equal to problem dim ...
	if (computePool.rank == 0 && (hdmrParam.maxOrder >= probParam.dim)) {
		cout << hline;
		cout << "\n Error: HDMR maxOrder (" << hdmrParam.maxOrder << ") must be < problem Dim (" << probParam.dim << ")\n";
		cout << hline;
		exit(0);
	}

	// set predfined xbar or calculate new xbar.
	if (hdmrParam.sampleSize < 1) {

		if (preDefinedXbar.size() == probParam.dim ) {
			for (int d = 0; d < probParam.dim; ++d) {
				data.xBar[d] += preDefinedXbar[d];
			}
			problemFunc(data.xBar, probParam.dim, data.fxBar);
		} else {
			cout << hline;
			cout << "\n Error: For 0 sampling Iteration xBar must be provided \n";
			cout << hline;
			exit(0);
		}
	} else {
		xbar();
	}

	writeJobAlloc();

	//General variables
	int index;
	unsigned pointsCount = 0;

	string fileName = "";
	string localIndexFileString = "";

	//Initialize fullPointMap Vector and copy xBar into it
	vector<double> fullPointMap(probParam.dim);
	memcpy(&fullPointMap[0], data.xBar, probParam.dim * sizeof(double));

	// Loop over every "Job Order" and go over the ith task in the "Job order"
	for (int k = 0; k < hdmrParam.maxOrder; ++k) {

		for (int i = 0; i < job[k].size; ++i) {

			// Allocate each task to specific block (The entire bock works on create the surplus files)
			if ( i % computePool.blockSize == computePool.blockID) {

				// Write order of gird (k is zero based so we add 1)
				if (computePool.threadID == 0) {
					localIndexFileString += to_string(k + 1) + "\t";
				}

				//Select the combination index from the task list and change fullPointMap also make fileName
				for (int j = 0; j < k + 1 ; ++j) {
					index = job[k].task[i * (k + 1) + j];
					fileName +=  to_string(index) + ".";
					fullPointMap[index - 1] = -1;

					if (computePool.threadID == 0) {
						localIndexFileString += to_string(index) + "\t";
					}
				}

				if (computePool.threadID == 0) {
					localIndexFileString += "\n";
				}

				fileName += "data";

				//Set the MPI environment for current operations
				sgwrite[k]->mpiCOMM = computePool.subComm;
				sgwrite[k]->rank = computePool.threadID;
				sgwrite[k]->size = computePool.tpb;

				//Build sparse grid & Adaptive for HDMR
				pointsCount += sgwrite[k]->build(fullPointMap);

				// Remove double counting
				if (computePool.threadID != 0) {
					pointsCount = 0;
				}

				// Write out surplus file & Get number of points, for sums we only keep threadID 0 - Everyone else is set to 0
				sgwrite[k]->write(folderName + fileName);

				// Reset and cleanup Variables
				fileName = "";
				sgwrite[k]->Cleanup();
				memcpy(&fullPointMap[0], data.xBar, probParam.dim * sizeof(double));
			}
		}
	}

	//	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &pointsCount, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

	//Write to file add delay to make synchronous (Better way is to write with rank=0 using MIP_gather)
	sleep(computePool.rank * 1);
	fstream output;

	// Make rank 0  write the header of the file
	if (computePool.rank == 0) {

		// Clear out file if it exists
		output.open(folderName + "index.data", ios::out | ios::trunc);

		output << probParam.dim << "\t";
		output <<  probParam.dof << "\t";

		output << hdmrParam.maxOrder << "\t";
		output <<  hdmrParam.sampleSize << "\t";
		output <<  hdmrParam.cutOff << "\t";

		output << sgParam.maxLevel << "\t";
		output << sgParam.gridType << "\t";
		output << sgParam.cutOff << "\t";

		output << pointsCount << "\n";

		// Output xBar
		for (int i = 0; i < probParam.dim; ++i) {output << setprecision(15) << data.xBar[i] << "\t";}

		output << "\n";

		// Output f(xBar)
		for (int i = 0; i < probParam.dof; ++i) {output << setprecision(15) << data.fxBar[i] << "\t";}

		output << "\n";

	} else {
		output.open(folderName + "index.data", ios::out | ios::app);
	}

	output << localIndexFileString;
	output.close();

	// Initialize Global Timer - End
	runTime.end = MPI_Wtime();

	//We set global barrier to make sure all all files have been written
	MPI_Barrier(MPI_COMM_WORLD);
	return pointsCount;
}

// Read Routines

inline void HDMR::traverseIterpolate(vector<int> alphabet, int k, double* fvalue , int op) {

	int comboCount = nCk(alphabet.size(), k);
	int comboBuffer[comboCount * k];
	vector<int> subAlphabet;
	string surplusFile = "";

	// Degrade alphabet into total of n=comboCount and of length k
	combination(alphabet, k, comboBuffer);

	// Cycle through each combination set
	for (int i = 0; i < comboCount; ++i) {

		// Make surplus file name & construct the correct x value
		surplusFile = "";
		for (int j = 0; j < (k); ++j) {
			surplusFile += to_string(comboBuffer[i * (k) + j]) + ".";
		}
		surplusFile += "data";

		// Apply the operation op on the new interpolation
		for (int j = 0; j < probParam.dof; ++j) {
			fvalue[j] +=  localCache.interpolate[surplusFile][j] * op ;
		}

		// Construct a new vector holding the sub alphabet of the above combination
		subAlphabet.clear();
		for (int j = 0; j < k; ++j) {
			subAlphabet.push_back(comboBuffer[i * k + j]);
		}

		for (int r = k - 1; r > 0; --r) {
			traverseIterpolate(subAlphabet, r, fvalue, (op * -1) );
		}

		// Deduct fxBar from fValeue
		for (int j = 0; j < probParam.dof ; ++j) {
			fvalue[j]  += (op * -1.0)  * data.fxBar[j];
		}
	}

	return;
}

void HDMR::interpolate(double* xSet, double* fvalueSet , int pointCount) {

	runTime.writeState = 0;
	runTime.start = MPI_Wtime();
	string surplusFile = "";
	vector<int> alphabet;
	int offset;
	vector<string> lookup;

	int workload = 0;
	for (int i = 0; i < computePool.workload.size(); ++i) {
		workload += computePool.workload[i];
	}

	// Set fvalueSet to zero
	for (int i = 0; i < pointCount * probParam.dof; ++i) { fvalueSet[i] = 0.0; }


	// Cycle through each interpolation x
	for (int s = 0; s < pointCount; ++s) {

		// Select the current x and fvalue
		double* x = &xSet[s * probParam.dim];
		double* fvalue = &fvalueSet[s * probParam.dof];
		double* temp = new double[workload * probParam.dof]();

		offset = 0;
		localCache.interpolate.clear();
		lookup.clear();

		for (int k = 1; k  <= hdmrParam.maxOrder ; ++k) {

			//sleep(computePool.rank * 2);
			double* xBuffer = new double[k];

			// Interpolate on root level (l=order)
			for (int i = 0; i < job[k - 1].size ; ++i) {

				// Construct surplusFile name and Select the correct index of X
				surplusFile = "";
				for (int j = 0; j < k  ; ++j) {
					xBuffer[j] = x[job[k - 1].task[i * k  + j] - 1];
					surplusFile += to_string(job[k - 1].task[i * k + j]) + ".";
				}
				surplusFile += "data";
				lookup.push_back(surplusFile);

				// Split block workload equally among threads (round robin)
				if ( i % computePool.tpb == computePool.threadID) {
					sgread[localCache.taskLookUpIndex[surplusFile]]->interpolateValue(xBuffer, &temp[probParam.dof * (offset + i)]);
				}
			}

			offset += job[k - 1].size;
			delete[] xBuffer;
		}

		MPI_Allreduce(
		    MPI_IN_PLACE,
		    temp,
		    workload * probParam.dof,
		    MPI_DOUBLE,
		    MPI_SUM,
		    MPI_COMM_WORLD);

		// Copy to local cache
		for (int i = 0; i < lookup.size(); ++i) {
			for (int j = 0; j < probParam.dof ; ++j) {
				localCache.interpolate[lookup[i]].push_back(temp[i * probParam.dof  + j]);
			}
		}

		//Cycle component of each root level (If sub level exist)
		for (int k = hdmrParam.maxOrder; k > 0; --k) {

			// Cycle through each job at the current level
			for (int i = 0; i < job[k - 1].size; ++i) {

				// Split block workload equally among threads
				if ( i % computePool.tpb == computePool.threadID) {

					// Construct alphabet for combinatoric
					alphabet.clear();
					for (int j = 0; j < k; ++j) {
						alphabet.push_back(job[k - 1].task[i * k  + j]);
					}

					traverseIterpolate(alphabet, k, fvalue, 1.0);
				}
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
		delete[] temp;
	}

	MPI_Allreduce(
	    MPI_IN_PLACE,
	    fvalueSet,
	    probParam.dof * pointCount,
	    MPI_DOUBLE,
	    MPI_SUM,
	    MPI_COMM_WORLD);

	// Add final fxBar
	for (int s = 0; s < pointCount; ++s) {
		for (int j = 0; j < probParam.dof ; ++j) {
			fvalueSet[s * probParam.dof + j] += data.fxBar[j];
		}
	}

	runTime.end = MPI_Wtime();
	return;
}

int HDMR::readJobAlloc(string folderName) {

	ifstream indexFile(folderName + "index.data");
	int k = 0;
	int taskCount = 0;

	if (indexFile.is_open()) {

		//Load problem parameters
		indexFile >> probParam.dim;
		indexFile >> probParam.dof;

		//Load HDMR parameters
		indexFile >>  hdmrParam.maxOrder;
		indexFile >>  hdmrParam.sampleSize;
		indexFile >>  hdmrParam.cutOff;

		indexFile >>  sgParam.maxLevel;
		indexFile >>  sgParam.gridType;
		indexFile >>  sgParam.cutOff;

		indexFile >>  gridPointCount;

		//Initialize xBar and f(xBar) array
		data.xBar = new double[probParam.dim]();
		data.fxBar = new double[probParam.dof]();

		//Load xBar
		for (int i = 0; i < probParam.dim; ++i) {
			indexFile >> data.xBar[i];
		}

		//Load fxBar
		for (int i = 0; i < probParam.dof; ++i) {
			indexFile >> data.fxBar[i];
		}
	}

	// Allocate job list
	job = new Job[hdmrParam.maxOrder];

	for (int k = 0; k < hdmrParam.maxOrder; ++k) {
		job[k].size = 0;
		job[k].task = new int[ nCk(probParam.dim, k + 1) * (k + 1) ]();
	}

	// Allocate tasks to correct job
	while (indexFile) {
		// Read the order of job and allocate task
		indexFile >> k;

		for (int i = 0; i < k; ++i) {
			indexFile >> job[k - 1].task[ job[k - 1].size * k + i];
		}

		//Increase size on each new line (negate one for the last newline at the EOF)
		job[k - 1].size ++;

		if (indexFile.peek() == EOF) {
			job[k - 1].size --;
		}
	}

	indexFile.close();

	// Get Total taskCount
	for (int k = 0; k < hdmrParam.maxOrder; ++k) {
		taskCount += job[k].size;
	}

	return taskCount;
}

void HDMR::read(string folderName) {

	// Initialize class state and timer
	runTime.writeState = 0;
	runTime.start = MPI_Wtime();

	// Allocate job form index file and initialize computePool
	readJobAlloc(folderName);

	// Allocate compute resources
	computeAlloc("bigPool");

	// General Variables
	string surplusFile = "";
	int offset = 0;

	// Instantiate SGread class for every for every group size
	sgread = new SGread* [computePool.workload[0]];

	for (int k = hdmrParam.maxOrder; k > 0; --k) {

		// Cycle through all job on each level
		for (int i = 0; i < job[k - 1].size; ++i) {

			// Construct surplus file name
			surplusFile = "";
			for (int j = 0; j < k ; ++j) {
				surplusFile += to_string(job[k - 1].task[i * k + j]) + ".";
			}
			surplusFile += "data";

			// Initialize localCache
			localCache.taskLookUpIndex[surplusFile] = offset + i;

			// Allocated every task to every thread (Even though we might not use that task on the thread)
			sgread[localCache.taskLookUpIndex[surplusFile]] = new SGread(0);

			//Set the MPI environment for current operations
			sgread[localCache.taskLookUpIndex[surplusFile]]->mpiCOMM = computePool.subComm;
			sgread[localCache.taskLookUpIndex[surplusFile]]->rank = 0;
			sgread[localCache.taskLookUpIndex[surplusFile]]->size = 1;
			sgread[localCache.taskLookUpIndex[surplusFile]]->read(folderName + surplusFile);
		}

		// Increment offset of sgread index by the last job size
		offset += job[k - 1].size;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	runTime.end = MPI_Wtime();
	return;
}






// Support Routines

void HDMR::xbar() {

	//Domain size
	double lower_bound = 0;
	double upper_bound = 1;

	//Dynamic domain size shift
	double lower_bound_shift = 0;
	double upper_bound_shift = 0;


	calc_xbar(lower_bound, upper_bound);

	MPI_Barrier(MPI_COMM_WORLD);

	//xBar right now is a sum over "numSums" and all "ranks"
	MPI_Allreduce(MPI_IN_PLACE, &data.xBar[0], probParam.dim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	// Now we normalize by dividing by "numSums" and all "ranks"
	for (int i = 0; i < probParam.dim; ++i) {
		data.xBar[i] = data.xBar[i] / double(computePool.size);
	}

	// Evaluate & Store f(xBar)
	problemFunc(&data.xBar[0], probParam.dim, &data.fxBar[0]);
}

void HDMR::calc_xbar( double lower_bound, double upper_bound) {

	//Array allocation for sampling
	double* sampleValue = new double [probParam.dof * hdmrParam.sampleSize]();
	double* sampleSpace = new double [probParam.dim * hdmrParam.sampleSize]();
	double* sampleAvg   = new double [probParam.dof]();

	double errNorm  = 0.0;
	double minNorm  = 0.0;
	int centroidInx = 0;

	//Generate flat/uniform distribution
	random_device rd;
	mt19937 generator(rd());
	uniform_real_distribution<double> distribution(lower_bound, upper_bound);

	//Find function mean for n="sampleSize" random points
	for (int i = 0; i < hdmrParam.sampleSize; ++i) {

		//Assemble Random vector
		for (int d = 0; d < probParam.dim; ++d) {
			sampleSpace[i * probParam.dim + d] = distribution(generator);
		}

		//Evaluate function at random vector
		problemFunc(&sampleSpace[i * probParam.dim], probParam.dim, &sampleValue[i * probParam.dof] );

		//Calculate mean
		for (int dof = 0; dof < probParam.dof ; ++dof) {
			sampleAvg[dof] += sampleValue[i * probParam.dof + dof] / double (hdmrParam.sampleSize);
		}
	}

	//Find coordinate of mean
	for (int i = 0; i < hdmrParam.sampleSize; ++i) {

		errNorm = 0.0;

		//Calculate L2 Norm
		for (int dof = 0; dof < probParam.dof ; ++dof) {
			errNorm += pow(sampleAvg[dof] - sampleValue[i * probParam.dof + dof], 2) ;
		}


		//If norm is less than before select the current index
		if (errNorm < minNorm || i == 0 ) {
			centroidInx = i;
			minNorm = errNorm;
		}
	}

	for (int d = 0; d < probParam.dim; ++d) {
		data.xBar[d] += sampleSpace[centroidInx * probParam.dim + d];
	}

	delete[] sampleAvg;
	delete[] sampleSpace;
	delete[] sampleValue;
}

double HDMR::l2norm(double* x, int size) {
	double norm = 0;

	for (int i = 0; i < size; ++i) {
		norm += x[i] * x[i];
	}

	return sqrt(norm);
}


// Combinatorial Routines
unsigned int HDMR::nCk(int n , int k) {
	if (k > n) {
		return 0;
	}

	if (k * 2 > n) {
		k = n - k;
	}

	if (k == 0) {
		return 1;
	}

	int result = n;

	for ( int i = 2; i <= k; ++i ) {
		result *= (n - i + 1);
		result /= i;
	}

	return result;
}

void HDMR::combination(vector<int> alphabet, int k, int* rtnBuffer) {

	//Setup vars & buffers
	int s = nCk(alphabet.size(), k);
	int* tempBuffer = new int[k]();
	int* offset = new int[1]();

	calc_combination(alphabet, &rtnBuffer[0], k, &tempBuffer[0], &offset[0]);

	delete tempBuffer;
	delete offset;
}

void HDMR::calc_combination(vector<int> alphabet, int* rtnBuffer, const int order , int* tempBuffer, int* offset,  int start, int index ) {

	if (index == order) {
		for (int j = 0; j < order; j++) {
			rtnBuffer[offset[0] + j] = tempBuffer[j];
		}

		offset[0] += order;
		return;
	}

	// replace index with all possible elements. The condition
	// "end-i+1 >= r-index" makes sure that including one element
	// at index will make a combination with remaining elements
	// at remaining positions
	for (int i = start; i <= alphabet.size() - 1 && alphabet.size() - i >= order - index; i++) {
		tempBuffer[index] = alphabet[i];
		calc_combination(alphabet, &rtnBuffer[0], order, &tempBuffer[0], &offset[0], i + 1, index + 1);
	}
}





// Debut rutines

void HDMR::debug(int showComputePool, int showProb, int showSG, int showHDMR, int showJob , int showOther) {

	MPI_Barrier(MPI_COMM_WORLD);
	sleep(computePool.rank * 2);

	cout << hline;
	cout << "HDMR Class State Dump [ Rank: " << computePool.rank << " ]\n";

	if (showComputePool) {
		cout << "\n> Compute Pool\n";
		cout << "  Rank: " << computePool.rank << endl;
		cout << "  Size: " << computePool.size << endl;
		cout << "  blockID: " << computePool.blockID << endl;
		cout << "  blockSize: " << computePool.blockSize << endl;
		cout << "  opt_tpb: " << computePool.opt_tpb << endl;
		cout << "  threadID: " << computePool.threadID << endl;
		cout << "  tbp: " << computePool.tpb << endl;

		if (!runTime.writeState) {
			cout << "  workload: [ ";

			for (int i = 0; i < computePool.workload.size(); ++i) {
				cout << computePool.workload[i] << " ";
			}

			cout << "]" << endl;
		}
	}

	if (showProb) {
		cout << "\n> Problem Param\n";
		cout << "  Dim: " << probParam.dim << endl;
		cout << "  DOF: " << probParam.dof << endl;
	}

	if (showSG) {
		cout << "\n> Spares Grid\n";
		cout << "  Lmax: " << sgParam.maxLevel << endl;
		cout << "  Type: " << sgParam.gridType << endl;
		cout << "  CutOff/Epsilon: " << sgParam.cutOff << endl;
	}

	if (showHDMR) {
		cout << "\n> HDMR Reduction\n";
		cout << "  maxOrder: " << hdmrParam.maxOrder << endl;
		cout << "  sampleSize: " << hdmrParam.sampleSize << endl;
		cout << "  CutOff/Epsilon: " << hdmrParam.cutOff << endl;
	}


	if (showJob) {
		cout << "\n> Job Allocation\n";

		for (int k = 0; k < hdmrParam.maxOrder; ++k) {
			cout << " Job Order: " << k + 1  << endl;
			cout << "  Size: " << job[k].size  << "/" << nCk(probParam.dim, k + 1) << endl;
			cout << "  Task: ";

			for (int i = 0; i < job[k].size; ++i) {
				cout << "{ ";

				for (int j = 0; j < k + 1; ++j) {
					cout << job[k].task[i * (k + 1) + j] << " ";
				}

				cout << "} ";
			}

			cout <<  endl;
		}


		cout << "\n> HDMR integral\n";
		for (auto i : localCache.integral) {
			cout << "  w(" << i.first << ")=" << i.second[0] << endl;
		}


		cout << "\n> HDMR BlackList\n";
		for (int i = 0; i < localCache.blackList.size(); ++i) {
			cout << "  " << localCache.blackList[i] << endl;
		}
	}

	if (showOther) {
		cout << "\n> Xbar & F(Xbar)\n";
		cout << "  Xbar : [ ";

		for (int i = 0; i < probParam.dim; ++i) {
			cout << data.xBar[i] << " ";
		}

		cout << "]" << endl;

		cout << "  f(Xbar) : [ ";

		for (int i = 0; i < probParam.dof; ++i) {
			cout << data.fxBar[i] << " ";
		}

		cout << "]" << endl;

	}

	cout << "\n> Run Time\n";
	cout << "  Wall-Time: " << runTime.end - runTime.start <<  endl;

	if (computePool.rank == computePool.size - 1) {
		cout << hline;
	}
}

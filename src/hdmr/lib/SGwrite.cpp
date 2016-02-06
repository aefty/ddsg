#include "SGwrite.h"
#include <stdio.h>
#include <array>
#include <string>
#include <algorithm>

using namespace std;

SGwrite::SGwrite( void (*problem_)(double*, int, double*), int dim_, int dof_, int Lmax_, double epsilon_, int gridType_ , int verbose_): AdaptiveSparseGrid(dim_, Lmax_, epsilon_, gridType_ ) {

	// Instantiate Values of AdaptiveSparseGrid
	dim = dim_;
	Lmax = Lmax_;
	epsilon = epsilon_;
	TotalDof = dof_;

	problem = problem_;
	verbose = verbose_;
	print = verbose;

	gridType = gridType_;

	//Allocate memory of variable  in AdaptiveSparseGrid library (required!)
	surplus = new double[TotalDof];
}

SGwrite::~SGwrite() {
	if (fullPointMap.size()) {
		delete[] fullPoint;
	}
}

int SGwrite::build(vector<double> fullPointMap_) {


	// Allocate fullPointMap, if fullPointMap={}, we make all points free
	fullPointMap = fullPointMap_;
	if (fullPointMap.size() == 0) {
		for (int i = 0; i < dim; ++i) {
			fullPointMap.push_back(-1.0);
		}
	}

	// Allocate fullPoint array
	fullPoint = new double[fullPointMap.size()]();

	double tm1, tm2;

	// Construct sparse Grid
	tm1 = MPI_Wtime();
	if (verbose > 1) { cout << hline; }
	BuildAdaptiveSparseGrid(print);
	if (verbose > 1) { cout << hline; }
	tm2 = MPI_Wtime();


	// Verbose output
	int NoPoint = NumberOfPoints();
	if (rank == 0 && verbose > 1) {
		cout << hline;
		cout << "fullPointMap: [ ";

		for (int i = 0; i < fullPointMap.size(); ++i) {
			cout << fullPointMap[i] << " ";
		}

		cout << "]" << endl;
		cout << "Basis Function Type : " << type << endl;
		cout << "Construct the adaptive sparse grid use (sec) :" << (tm2 - tm1) << endl;
		cout << "The number of the grid is : " << NoPoint << endl;
		cout << "The depth  of the grid is : " << InterpolationLevel() << endl;
		cout << hline;
	}

	return NoPoint;
}

void SGwrite::write(string fileName) {

	char* fileName_char =  new char[fileName.length() + 1];
	strcpy(fileName_char, fileName.c_str());

	// Write Surplus
	AdaptiveARRAY<int> index;
	index.redim(TotalDof);
	for ( int i = 1; i <= TotalDof; i++) {
		index(i) = i - 1;
	}
	StoreSurplus(index, fileName_char);

	delete[] fileName_char;
}

void SGwrite::integrateDomain(double* fvalue, double op) {

	if (op == 0.0) {// Write over value
		SpIntegrate( &fvalue[0]);
	} else { // Add or Subtract value base don op
		double tempFval[TotalDof];
		SpIntegrate(&tempFval[0]);
		for (int i = 0; i < TotalDof; ++i) {
			fvalue[i] += tempFval[i] * op;
		}
	}
}

void SGwrite::EvaluateFunctionAtThisPoint( AdaptiveARRAY<double>* x) {

	//Combine AdaptiveSparseGrid grid point in R^Dim [0,1] with fullPointMap
	int index = 0;
	for (int i = 0; i < fullPointMap.size(); ++i) {

		if (fullPointMap[i] < 0) {
			fullPoint[i] = x->pData[index];
			index++;
		} else {
			fullPoint[i] = fullPointMap[i];
		}
	}

	problem(&fullPoint[0], fullPointMap.size(), surplus);
}
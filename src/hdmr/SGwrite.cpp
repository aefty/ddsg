#include "SGwrite.h"
#include <stdio.h>
#include <array>
#include <string>
#include <cstring>
#include <algorithm>

using std::cout;
using std::string;
using std::vector;
using std::map;


SGwrite::SGwrite(
    void (*problem_)(double*, int, double*),
    int dim_,
    int dof_,
    int Lmax_,
    double epsilon_,
    int gridType_ ,
    int verbose_): AdaptiveSparseGrid(dim_, Lmax_, epsilon_, gridType_ ) {

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

	// Initilize masked x value ... this will be updated in the case we use build(activeDim)
	maskedX.resize(dim);

	// Initilize xbar value ... this will be updated in the case we use build(activeDim)
	xBar.resize(dim);
}

SGwrite::~SGwrite() {

}


// #AE
void SGwrite::resetMPI(MPI_Comm mpiCOMM_) {

	mpiCOMM = mpiCOMM_;
	MPI_Comm_rank(mpiCOMM, &rank);
	MPI_Comm_size(mpiCOMM, &size);
}


void SGwrite::setAnchor(vector<double>& xBar_) {
	// Set anchor point
	xBar.resize(xBar_.size());
	xBar = xBar_;
}

void SGwrite::build(vector<int>& activeDim_) {

	// copy active Dimentions and xBar
	activeDim.resize(activeDim_.size());
	activeDim = activeDim_;

	// Initilize masked x value
	maskedX.resize(xBar.size());

	build();
}



void SGwrite::build() {

	double tm1, tm2;

	// If activeDim is not set, than set every element as active
	if (activeDim.empty()) {
		activeDim.resize(dim);
		for (int i = 0; i < activeDim.size(); ++i) {
			activeDim[i] = i;
		}
	}

	// Construct sparse Grid
	tm1 = MPI_Wtime();
	if (rank == 0 && verbose)  { cout << hline; }
	BuildAdaptiveSparseGrid(print);
	if (rank == 0 && verbose) { cout << hline; }
	tm2 = MPI_Wtime();

	// Verbose output
	int NoPoint = NumberOfPoints();

	if (rank == 0 && verbose) {
		cout << hline;
		cout << "activeDim: [ ";

		for (int i = 0; i < activeDim.size(); ++i) {
			cout << activeDim[i] << " ";
		}

		cout << "]" << endl;
		cout << "Basis Function Type : " << gridType << endl;
		cout << "Construct the adaptive sparse grid use (sec) :" << (tm2 - tm1) << endl;
		cout << "The number of the grid is : " << NoPoint << endl;
		cout << "The depth  of the grid is : " << InterpolationLevel() << endl;
		cout << hline;
	}
}




int SGwrite::write(string surplusFileName) {

	char surplusFileName_char [surplusFileName.length() + 1];
	strcpy(surplusFileName_char, surplusFileName.c_str());

	AdaptiveARRAY<int> index;
	index.redim(TotalDof);
	for ( int i = 1; i <= TotalDof; i++) {
		index(i) = i - 1;
	}
	StoreSurplus(index, surplusFileName_char);

	return NumberOfPoints();
}
void SGwrite::integrateDomain(double* fvalue) {

	SpIntegrate( &fvalue[0]);
}

void SGwrite::EvaluateFunctionAtThisPoint( AdaptiveARRAY<double>* x ) {

	// Copy points to masked value
	maskedX = xBar;

	/*
		cout << "xBar : ";
		for (int i = 0; i < maskedX.size(); ++i) {
			cout << maskedX[i] << " ";
		}
		cout << endl;




		cout << "activeDim : ";
		for (int i = 0; i < activeDim.size(); ++i) {
			cout << activeDim[i] << " ";
		}
		cout << endl;




		cout << "Input : ";
		for (int i = 0; i < dim; ++i) {
			cout << x->pData[i] << " ";
		}
		cout << endl;
	*/

	// apply the mask
	for (int i = 0; i < activeDim.size(); ++i) {
		maskedX[activeDim[i]] = x->pData[i];
	}

	/*
		cout << "maskedX : ";
		for (int i = 0; i < maskedX.size(); ++i) {
			cout << maskedX[i] << " ";
		}
		cout << endl << endl;
	*/

	problem(&maskedX[0], xBar.size(), surplus);
}

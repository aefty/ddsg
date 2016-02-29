#include "SGread.h"
#include <stdio.h>
#include <array>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;

SGread::SGread (int verbose_) : Post() {
	verbose = verbose_;
}

SGread::~SGread() {
	surplusFileName = "";
}




void SGread::resetMPI(MPI_Comm mpiCOMM_) {
	mpiCOMM = mpiCOMM_;
	MPI_Comm_rank(mpiCOMM, &rank);
	MPI_Comm_size(mpiCOMM, &size);
}



int SGread::read(string surplusFileName) {

	char surplusFileName_char [surplusFileName.length() + 1];
	strcpy(surplusFileName_char, surplusFileName.c_str());
	LoadData(surplusFileName_char);

	return nno;
}


void SGread::interpolateValue(double* x, double* fvalue, double op) {

	if (op == 0.0) {// Write over value
		Interpolate(x, fvalue);
	} else { // Add or Subtract value base don op
		double tempFval[TotalDof];
		Interpolate(x, tempFval);

		for (int i = 0; i < TotalDof; ++i) {
			fvalue[i] += tempFval[i] * op;
		}
	}
	return;
}

void SGread::integrateDomain(double* fvalue, double op) {

	if (op == 0.0) {// Write over value
		Integrate( &fvalue[0]);
	} else { // Add or Subtract value base don op
		double tempFval[TotalDof];
		Integrate(&tempFval[0]);
		for (int i = 0; i < TotalDof; ++i) {
			fvalue[i] += tempFval[i] * op;
		}
	}

	return;
}
#include "SGread.h"
#include <stdio.h>
#include <array>
#include <string>
#include <cstring>
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


void SGread::interpolateValue(double* x, double* fvalue) {
	Interpolate(x, fvalue);
}

void SGread::integrateDomain(double* fvalue) {
	Integrate(fvalue);
}
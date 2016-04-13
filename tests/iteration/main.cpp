#include <iostream>
#include <ostream>
#include <vector>
#include <math.h>
#include <random>
#include <iomanip> // setprecision

#include "../../TestFunctions/TestFunctions.cpp"
#include "../../src/include.h"

int main(int argc, char* argv[]) {


	MPI_Init (&argc, &argv);

	// Input form command line
	string folderName = argv[1];


	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	HDMR* hdmr = new HDMR();

	hdmr->read(folderName);
	hdmr->debug("debug");

	int numberOfpoints = 1;

	vector<double> x(hdmr->probParam.dim, 0.5 );

	double* fval = new double[hdmr->probParam.dof * numberOfpoints];
	//generateRandom(x, hdmr->probParam.dim, 1);

	for (int i = 0; i < 10; ++i) {
		hdmr->interpolate( &x[0], fval, numberOfpoints);


		cout << "f(";
		for (int r = 0; r < hdmr->probParam.dim; ++r) {
			cout << x[r] <<  " ";
		}
		cout << ")=";


		cout << "=[ ";
		for (int r = 0; r < hdmr->probParam.dof; ++r) {
			cout << fval[r] <<  " ";
		}
		cout << "]" << endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	//delete[] x;
	delete[] fval;
	delete hdmr;

	MPI_Finalize ();
	return 0;
}





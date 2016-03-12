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

	// Problem Parameters
	int dim = 100;
	int dof = 1;

	// HMDR Parameters
	int HDMRmaxOrder = 2;
	double HDMRcutOff = 0.00001;

	// SG Parameters
	int SGmaxLevel = 7;
	double SGcutOff = 0.0;
	int SGgridType = 1;

	int processPerGroup = 0;
	string method = "";

	// Input form command line
	method = argv[1];
	processPerGroup = atoi(argv[2]);
	dim = atoi(argv[3]);
	SGmaxLevel = atoi(argv[4]);
	SGgridType = atoi(argv[5]);
	HDMRmaxOrder = atoi(argv[6]);
	HDMRcutOff = atof(argv[7]);

	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Genearte test points
	int numberOfpoints = 10000;
	double* x = new double[dim * numberOfpoints];
	double* fval = new double[dof * numberOfpoints];
	generateRandom(x, dim, numberOfpoints);

	HDMR* hdmr = new HDMR();
	vector<double> xBar(dim, 0.5);

	if (method == "sg") {

		hdmr->write(f4, dim, dof, SGmaxLevel, SGcutOff, SGgridType);
		hdmr->debug("SG.WRITE", 0, 1, 1, 1, 1, 0);

		hdmr->read("surplusSG/");
		//	hdmr->debug("SG.READ");

		hdmr->interpolate(x, fval, numberOfpoints);
		hdmr->debug("SG.INTERPOLATE");
	} else if (method == "hdmr") {
		hdmr->write(f4, dim, dof, SGmaxLevel, SGcutOff, SGgridType, HDMRmaxOrder, HDMRcutOff, xBar, processPerGroup);
		hdmr->debug("HDMR.WRITE", 0, 1, 0, 0, 1, 1);

		hdmr->read("surplusHDMR/");
		//	hdmr->debug("HDMR.READ");

		hdmr->interpolate(x, fval, numberOfpoints);
		hdmr->debug("HDMR.INTERPOLATE",  0, 1, 0, 0, 0, 1);
	} else {
		cout << "Invalid input method ... sg or hdmr" << endl;
		MPI_Finalize ();
		return 0;
	}


	if (rank == 0) {
		double error =  checkError(f4, x, dim, fval, dof, numberOfpoints, 0);
		cout << "Total Percentage Error Avg : " << error << endl;
	}


	MPI_Barrier(MPI_COMM_WORLD);

	delete[] x;
	delete[] fval;
	delete hdmr;

	MPI_Finalize ();
	return 0;
}





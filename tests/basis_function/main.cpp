#include <iostream>
#include <ostream>
#include <vector>
#include <math.h>
#include <random>
#include <iomanip> // setprecision

#include "../../TestFunctions/TestFunctions.cpp"
#include "../../src/include.h"


void funNull(double* x ,  int xDim, double* val ) {
	val[0] = x[0];
}

int main(int argc, char* argv[]) {


	MPI_Init (&argc, &argv);

	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Not of this matters, its just o intilize the class
	int dim = 10;
	int dof = 1;
	double SGcutOff = 0.0;
	int SGmaxLevel = 3;
	int SGgridType = 1;

	// Genearl Variables
	int level = 1;
	double val;

	// Step size for plot density
	double x, xp;
	double n = 5;
	double h;

	// Set variables
	if (argc == 3) {
		SGgridType = atoi(argv[1]);
		SGmaxLevel = atoi(argv[2]);
	} else if (argc == 4) {
		SGgridType = atoi(argv[1]);
		SGmaxLevel = atoi(argv[2]);
		n = atoi(argv[3]);

	} else {
		cout << "INPUT SHOULD BE : <SGgridType><SGmaxLevel><OPTION points_in_X>" << endl;
		return 0;
	}

	h = pow(2, -n);



	// Intilize SG
	SGwrite* sgwrite = new SGwrite(funNull, dim , dof , SGmaxLevel, SGcutOff, SGgridType , 0);
	sgwrite->mpiCOMM = MPI_COMM_WORLD;
	sgwrite->rank = rank;
	sgwrite->size = size;


	std::cout.precision(4);
	std::cout.setf( std::ios::fixed, std:: ios::floatfield );


	cout << "GRID TYPE :" << SGgridType << endl;
	cout << "Max Level :" << SGmaxLevel << endl;

	x = 0.0;
	cout << "Level: " << level << " [index|NodalPoint|x]" << endl;

	while (x <= 1.0) {
		val = sgwrite->BasisFunction(x, 1, 0);
		xp = sgwrite->IndextoCoordinate(1, 0);
		if (val > 0) {
			cout <<  1 << "|" << xp << "|" << x << " ";
			for (double i = 0; i < val; i += 0.05) {
				cout << ".";
			}
			cout << val << endl;
		}

		x += h;
	};


	level++;

	if (SGgridType == 1 || SGgridType == 4 ) {
		int l = level;
		cout << "=================================================" << endl;
		cout << "Level: " << 2 << " [index|NodalPoint|x]" << endl;
		cout << "=================================================" << endl;
		for (int  j = 1; j < pow(2, l) ; j += 2) {

			x = 0.0;

			while (x <= 1.0) {
				val = sgwrite->BasisFunction(x, l, j);
				xp = sgwrite->IndextoCoordinate(l, j);
				if (val > 0 ) {
					cout <<  j << "|" << xp << "|" << x << " ";
					for (double i = 0; i < val; i += 0.05) {
						cout << ".";
					}
					cout << val << endl;
				}

				x += h;
			}
		}
		cout << endl;
		level++;
	}


	for (int l = level; l <= SGmaxLevel; ++l) {
		cout << "================================================+" << endl;
		cout << "Level:" << l << " [index|NodalPoint|x]" << endl;
		cout << "=================================================" << endl;
		for (int  j = 2; j <= pow(2, l); j += 2) {
			x = 0.0;
			xp = sgwrite->IndextoCoordinate(l, j);

			while (x <= 1.0 ) {
				val = sgwrite->BasisFunction(x, l, j);

				if (val > 0 ) {
					cout <<  j << "|" << xp << "|" << x << " ";
					for (double i = 0; i < val; i += 0.05) {
						cout << ".";
					}
					cout << val << endl;
				}
				x += h;
			}
		}
		cout << endl;
	}



	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize ();

	delete sgwrite;

	return 0;
}





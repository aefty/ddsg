/**
 * Collection of test function
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <random>
#include <string>
#include <map>
#include <iomanip>
#include <algorithm>
#include <unistd.h>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::fill;

double PI =  4 * atan(1.0);

// Sum of squares
void f0(double* x ,  int dim, double* val) {
	val[0] = 0.0;

	for (int i = 0; i < dim; ++i) {
		val[0] += x[i] * x[i];
	};


	val[1] = 0.0;

	for (int i = 0; i < dim; ++i) {
		val[1] += x[i] * x[i] * x[i];
	};

}

// OSCILLATORY
void f1(double* x ,  int dim, double* val ) {
	val[0] = 0.0;

	for (int i = 0; i < dim; ++i) {
		val[0] += 5.0 * x[i];
	}

	val[0] = cos(2.0 * PI * 0.5 + val[0]) ;
}

// PRODUCT PEAK
void f2(double* x ,  int dim, double* val ) {
	val[0] = 1.0;

	for (int i = 0; i < dim; ++i) {
		val[0] *= 1.0 / (1.0 / (5.0 * 5.0) + (x[i] - 0.5) * (x[i] - 0.5)) ;
	}
}

// CORNER PEAK
void f3(double* x ,  int dim, double* val ) {
	val[0] = 0.0;

	for (int i = 0; i < dim; ++i) {
		val[0] += 5.0 * x[i] ;
	}

	val[0] = pow((1.0 + val[0]), -1.0 * (dim + 1));
}

// GAUSSIAN
void f4(double* x ,  int dim, double* val ) {
	val[0] = 0.0;

	for (int i = 0; i < dim; ++i) {
		val[0] += 5.0 * 5.0 * (x[i] - 0.5) * (x[i] -  0.5) ;
	}

	val[0] = exp(-1.0 * val[0]);
}

// LOG GAUSSIAN
void f4_1(double* x ,  int dim, double* val ) {
	val[0] = 0.0;

	for (int i = 0; i < dim; ++i) {
		val[0] += 5.0 * 5.0 * (x[i] - 0.5) * (x[i] -  0.5) ;
	}

	val[0] = -1.0 * val[0];
}

// GAUSSIAN NON-SYMETRICAL
void f4_2(double* x ,  int dim, double* val ) {
	val[0] = 0.0;

	for (int i = 0; i < dim; ++i) {
		val[0] += ( 5.0  * (double(i +  1.0) / dim)) * ( 5.0 * (double(i +  1.0) / dim))  * (x[i] - 0.5) * (x[i] -  0.5) ;
	}

	val[0] = exp(-1.0 * val[0]);
}

// CONTINUOUS
void f5(double* x ,  int dim, double* val ) {
	val[0] = 0.0;

	for (int i = 0; i < dim; ++i) {

		val[0] += 5 * fabs(x[i] - 0.5) ;
	}

	val[0] = exp(-1.0 * val[0]);
}

// MAX
void f6(double* x ,  int dim, double* val ) {
	val[0] = 0.0;
	double offset = 0.0;

	for (int i = 0; i < dim; ++i) {
		val[0] += fmax(0.31 * (x[i] + offset),  (x[i] + offset) *  (x[i] + offset)) ;
	}
}


void generateRandom(double* x, int dim, int numberOfpoints, int SEED = 100) {

	std::default_random_engine rd (SEED);
	std::uniform_real_distribution<double> dist(0.0, 1.0);

	for (int i = 0; i < numberOfpoints; ++i) {

		for (int j = 0; j < dim; ++j) {
			x[i * dim + j] = dist(rd);
		}

	}
}


double checkError(
    void (*problemFunc)(double*, int, double*),
    double* x,
    int dim,
    double* fval,
    int dof,
    int numberOfpoints,
    int verbose_percision = 0) {

	vector<double> errorVec(dof, 0.0);;
	double error = 0.0;
	double exactVal[dof];

	if (verbose_percision > 0) {
		cout << "Average Error Measure" << endl;
		cout << "=======================" << endl;
		std::cout.precision(verbose_percision);
		std::cout.setf( std::ios::fixed, std:: ios::floatfield );
	}

	for (int i = 0; i < numberOfpoints; ++i) {

		if (verbose_percision > 0) {
			cout << "f(" << x[i * dim ];
			for (int j = 1; j < dim; ++j) {
				cout << "," << x[i * dim + j];
			}
			cout << ") = ";
		}

		problemFunc(x + i * dim, dim, exactVal);

		for (int j = 0; j < dof; ++j) {
			errorVec[j] += fabs( fval[i * dof + j] - exactVal[j] ) / fabs(exactVal[j])  ;
		}

		if (verbose_percision > 0) {
			cout << "[" << fval[i * dof + 0];
			for (int j = 1; j < dof; ++j) {
				cout << " "  << fval[i * dof + j];
			}
			cout << "]";

			cout << " [" << exactVal[0];
			for (int j = 1; j < dof; ++j) {
				cout << " "  << exactVal[j];
			}
			cout << "] ";

			cout << "[" << fabs( fval[i * dof + 0] - exactVal[0] ) / fabs(exactVal[0]);
			for (int j = 1; j < dof; ++j) {
				cout << " "  << fabs( fval[i * dof + j] - exactVal[j] ) / fabs(exactVal[j]);
			}
			cout << "] " << endl;
		}

	}

	for (int j = 0; j < dof; ++j) {
		error += errorVec[j] ;
	}

	return error / numberOfpoints;
}
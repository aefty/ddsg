#include <iostream>
#include <ostream>
#include <vector>
#include <math.h>
#include <random>
#include <iomanip> // setprecision

#include "../../core/HDMR/lib/SGwrite.cpp"
#include "../../core/HDMR/lib/SGread.cpp"

#include "../../core/HDMR/lib/HDMR.cpp"

using namespace std;

// GLOBAL VARIABLES
int mpirank;
int mpisize;
string outputFolder;
double PI =  4 * atan(1.0);
int SEED = 100;
double** c;
double** w;
double* x;

int numFunc = 6;

// OSCILLATORY
void f1(double* x ,  int xDim, double* val ) {
	val[0] = 0.0;

	for (int i = 0; i < xDim; ++i) {
		val[0] += c[0][i] * x[i] ;
	}

	val[0] = cos(2.0 * PI * w[0][0] + val[0]) ;
}

// PRODUCT PEAK
void f2(double* x ,  int xDim, double* val ) {
	val[0] = 1.0;

	for (int i = 0; i < xDim; ++i) {
		val[0] *= 1.0 / (1.0 / (c[1][i] * c[1][i]) + (x[i] - w[1][i]) * (x[i] - w[1][i])) ;
	}
}

// CORNER PEAK
void f3(double* x ,  int xDim, double* val ) {
	val[0] = 0.0;

	for (int i = 0; i < xDim; ++i) {
		val[0] += c[2][i] * x[i] ;
	}

	val[0] = pow((1.0 + val[0]), -1.0 * (xDim + 1));
}

// LOG GAUSSIAN
void f4(double* x ,  int xDim, double* val ) {
	val[0] = 0.0;

	for (int i = 0; i < xDim; ++i) {
		val[0] += c[3][i] * c[3][i] * (x[i] - w[3][i]) * (x[i] - w[3][i]) ;
	}

	val[0] = -1.0 * val[0];
}

// CONTINUOUS
void f5(double* x ,  int xDim, double* val ) {
	val[0] = 0.0;

	for (int i = 0; i < xDim; ++i) {
		val[0] += c[4][i] * fabs(x[i] - w[4][i]) ;
	}

	val[0] = exp(-1.0 * val[0]);
}

// MAX
void f6(double* x ,  int xDim, double* val ) {
	val[0] = 0.0;
	double offset = 0.0;

	for (int i = 0; i < xDim; ++i) {
		val[0] += max(0.31 * (x[i] + offset),  (x[i] + offset) *  (x[i] + offset)) ;
	}
}

void testWriteHDMR(int dim, int dof, int HDMRmaxOrder, int HDMRsampleSize, double HDMRcutOff, int SGmaxLevel, int SGgridType, double SGcutOff) {

	if (mpirank == 0) {
		cout << "\n+HDMR Write: Start" << endl;
		cout << "==================" << endl;
	}

	int* HDMRinterpolationPoints = new int[numFunc]();
	vector<double> v(dim, 0.5);

	if (mpirank == 0) {
		cout << "+HDMR Write: ";
	}

	// HDMR Code
	HDMR* writeModelf1 =  new HDMR(1);
	HDMR* writeModelf2 =  new HDMR(1);
	HDMR* writeModelf3 =  new HDMR(1);
	HDMR* writeModelf4 =  new HDMR(1);
	HDMR* writeModelf5 =  new HDMR(1);
	HDMR* writeModelf6 =  new HDMR(1);

	HDMRinterpolationPoints[0] = writeModelf1->write(f1 , dim , dof , HDMRmaxOrder ,  HDMRsampleSize ,  HDMRcutOff , SGmaxLevel ,  SGgridType ,  SGcutOff , "surplus_hdmr/f1/", v);
	HDMRinterpolationPoints[1] = writeModelf2->write(f2 , dim , dof , HDMRmaxOrder ,  HDMRsampleSize ,  HDMRcutOff , SGmaxLevel ,  SGgridType ,  SGcutOff , "surplus_hdmr/f2/", v);
	HDMRinterpolationPoints[2] = writeModelf3->write(f3 , dim , dof , HDMRmaxOrder ,  HDMRsampleSize ,  HDMRcutOff , SGmaxLevel ,  SGgridType ,  SGcutOff , "surplus_hdmr/f3/", v);
	HDMRinterpolationPoints[3] = writeModelf4->write(f4 , dim , dof , HDMRmaxOrder ,  HDMRsampleSize ,  HDMRcutOff , SGmaxLevel ,  SGgridType ,  SGcutOff , "surplus_hdmr/f4/", v);
	HDMRinterpolationPoints[4] = writeModelf5->write(f5 , dim , dof , HDMRmaxOrder ,  HDMRsampleSize ,  HDMRcutOff , SGmaxLevel ,  SGgridType ,  SGcutOff , "surplus_hdmr/f5/", v);
	HDMRinterpolationPoints[5] = writeModelf6->write(f6 , dim , dof , HDMRmaxOrder ,  HDMRsampleSize ,  HDMRcutOff , SGmaxLevel ,  SGgridType ,  SGcutOff , "surplus_hdmr/f6/", v);


	MPI_Barrier(MPI_COMM_WORLD);

	if (0) {
		writeModelf1->debug();
		writeModelf2->debug();
		writeModelf3->debug();
		writeModelf4->debug();
		writeModelf5->debug();
		writeModelf6->debug();
	}

	if (mpirank == 0) {
		cout << "Successful" << endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (mpirank == 0) {
		// gridpoints Write
		fstream gridpointsWrite_output;
		gridpointsWrite_output.open(outputFolder + "/gridpoints_hdmr.csv", ios::out | ios::trunc);

		for (int i = 0; i < numFunc; ++i) {
			gridpointsWrite_output << HDMRinterpolationPoints[i] << "\t";
		}
		gridpointsWrite_output.close();

		cout << "Output gridpoints_hdmr.csv : Successful" << endl;
	}

	delete[] HDMRinterpolationPoints;

	delete writeModelf1;
	delete writeModelf2;
	delete writeModelf3;
	delete writeModelf4;
	delete writeModelf5;
	delete writeModelf6;
}

void testWriteSG(int dim, int dof, int SGmaxLevel, int SGgridType, double SGcutOff) {

	if (mpirank == 0) {
		cout << "\n+SG Write: Start" << endl;
		cout << "==================" << endl;
	}

	int* SGinterpolationPoints   = new int[numFunc]();
	vector<double> v(dim, 0.5);


	if (mpirank == 0) {
		cout << "+SG Write: ";
	}

	// SGWrite Code
	SGwrite* sgwritef1 = new SGwrite(f1, dim , dof , SGmaxLevel, SGcutOff, SGgridType , 0);
	SGinterpolationPoints[0] = sgwritef1->build();
	sgwritef1->write("surplus_sg/f1/surplus.data");

	SGwrite* sgwritef2 = new SGwrite(f2, dim , dof , SGmaxLevel, SGcutOff, SGgridType , 0);
	SGinterpolationPoints[1] = sgwritef2->build();
	sgwritef2->write("surplus_sg/f2/surplus.data");

	SGwrite* sgwritef3 =  new SGwrite(f3, dim , dof , SGmaxLevel, SGcutOff, SGgridType , 0);
	SGinterpolationPoints[2] = sgwritef3->build();
	sgwritef3->write("surplus_sg/f3/surplus.data");

	SGwrite* sgwritef4 =  new SGwrite(f4, dim , dof , SGmaxLevel, SGcutOff, SGgridType , 0);
	SGinterpolationPoints[3] = sgwritef4->build();
	sgwritef4->write("surplus_sg/f4/surplus.data");

	SGwrite* sgwritef5 =  new SGwrite(f5, dim , dof , SGmaxLevel, SGcutOff, SGgridType , 0);
	SGinterpolationPoints[4] = sgwritef5->build();
	sgwritef5->write("surplus_sg/f5/surplus.data");

	SGwrite* sgwritef6 =  new SGwrite(f6, dim , dof , SGmaxLevel, SGcutOff, SGgridType , 0);
	SGinterpolationPoints[5] = sgwritef6->build();
	sgwritef6->write("surplus_sg/f6/surplus.data");

	if (mpirank == 0) {
		cout << "Successful" << endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (mpirank == 0) {
		// gridpoints Write
		fstream gridpointsWrite_output;
		gridpointsWrite_output.open(outputFolder + "/gridpoints_sg.csv", ios::out | ios::trunc);

		for (int i = 0; i < numFunc; ++i) {
			gridpointsWrite_output << SGinterpolationPoints[i] << "\t" ;
		}
		gridpointsWrite_output.close();

		cout << "Output gridpoints_sg.csv : Successful" << endl;
	}

	delete[] SGinterpolationPoints;

	delete sgwritef1;
	delete sgwritef2;
	delete sgwritef3;
	delete sgwritef4;
	delete sgwritef5;
	delete sgwritef6;
}


void testReadHDMR(int dim, int dof,  int pointCout, int HDMRmaxOrder, int HDMRsampleSize, double HDMRcutOff, int SGmaxLevel, int SGgridType, double SGcutOff) {

	if (mpirank == 0) {
		cout << "\n+HDMR Read: Start" << endl;
		cout << "==================" << endl;
	}

	int* HDMRinterpolationPoints = new int[numFunc]();
	double* interValHDMR  = new double[dof * pointCout * numFunc]();

	if (mpirank == 0) {
		cout << "+HDMR Read: ";
	}

	HDMR* readModelf1 = new HDMR(1);
	HDMR* readModelf2 = new HDMR(1);
	HDMR* readModelf3 = new HDMR(1);
	HDMR* readModelf4 = new HDMR(1);
	HDMR* readModelf5 = new HDMR(1);
	HDMR* readModelf6 = new HDMR(1);

	readModelf1->read("surplus_hdmr/f1/");
	readModelf2->read("surplus_hdmr/f2/");
	readModelf3->read("surplus_hdmr/f3/");
	readModelf4->read("surplus_hdmr/f4/");
	readModelf5->read("surplus_hdmr/f5/");
	readModelf6->read("surplus_hdmr/f6/");

	if (mpirank == 0) {
		cout << "Successful" << endl;
	}


	MPI_Barrier(MPI_COMM_WORLD);

	if (0) {
		readModelf1->debug();
		readModelf2->debug();
		readModelf3->debug();
		readModelf4->debug();
		readModelf5->debug();
		readModelf6->debug();
	}

	if (mpirank == 0) {
		cout << "+HDMR interpolate: ";
	}

	readModelf1->interpolate(x, &interValHDMR[0 * pointCout * dof], pointCout);
	readModelf2->interpolate(x, &interValHDMR[1 * pointCout * dof], pointCout);
	readModelf3->interpolate(x, &interValHDMR[2 * pointCout * dof], pointCout);
	readModelf4->interpolate(x, &interValHDMR[3 * pointCout * dof], pointCout);
	readModelf5->interpolate(x, &interValHDMR[4 * pointCout * dof], pointCout);
	readModelf6->interpolate(x, &interValHDMR[5 * pointCout * dof], pointCout);


	if (mpirank == 0) {
		cout << "Successful" << endl;
	}

	if (mpirank == 0) {

		// hdmrInterpolate
		fstream hdmrInterpolate_output;
		hdmrInterpolate_output.open(outputFolder + "/hdmrInterpolate.csv", ios::out | ios::trunc);
		for (int i = 0; i < pointCout; ++i) {
			for (int f = 0; f < numFunc; ++f) {
				for (int d = 0; d < dof; ++d) {
					hdmrInterpolate_output << setprecision(15) << interValHDMR[f * pointCout * dof + i * dof + d] << "\t";
				}
			}
			hdmrInterpolate_output << "\n";
		}

		cout << "Output hdmrInterpolate.csv : Successful" << endl;
		hdmrInterpolate_output.close();
	}

	delete[] HDMRinterpolationPoints;
	delete[] interValHDMR;

	delete readModelf1;
	delete readModelf2;
	delete readModelf3;
	delete readModelf4;
	delete readModelf5;
	delete readModelf6;

}


void testReadSG(int dim, int dof,  int pointCout, int SGmaxLevel, int SGgridType, double SGcutOff) {

	if (mpirank == 0) {
		cout << "\n+SG Read: Start" << endl;
		cout << "==================" << endl;
	}

	int* SGinterpolationPoints   = new int[numFunc]();

	double* exactVal      = new double[dof * pointCout * numFunc]();
	double* interValSG    = new double[dof * pointCout * numFunc]();

	if (mpirank == 0) {
		cout << "+Exact Value Calc:";
	}
	for (int i = 0; i < pointCout; ++i) {
		f1(&x[i * dim], dim, &exactVal[0 * pointCout * dof + i * dof]);
		f2(&x[i * dim], dim, &exactVal[1 * pointCout * dof + i * dof]);
		f3(&x[i * dim], dim, &exactVal[2 * pointCout * dof + i * dof]);
		f4(&x[i * dim], dim, &exactVal[3 * pointCout * dof + i * dof]);
		f5(&x[i * dim], dim, &exactVal[4 * pointCout * dof + i * dof]);
		f6(&x[i * dim], dim, &exactVal[5 * pointCout * dof + i * dof]);
	}

	if (mpirank == 0) {
		cout << "Successful" << endl;
	}

	if (mpirank == 0) {
		cout << "+SG Read: ";
	}

	SGread* sgreadModelf1 = new SGread(0);
	sgreadModelf1->mpiCOMM = MPI_COMM_WORLD;
	sgreadModelf1->rank = mpirank;
	sgreadModelf1->size = mpisize;

	SGread* sgreadModelf2 = new SGread(0);
	sgreadModelf2->mpiCOMM = MPI_COMM_WORLD;
	sgreadModelf2->rank = mpirank;
	sgreadModelf2->size = mpisize;

	SGread* sgreadModelf3 = new SGread(0);
	sgreadModelf3->mpiCOMM = MPI_COMM_WORLD;
	sgreadModelf3->rank = mpirank;
	sgreadModelf3->size = mpisize;

	SGread* sgreadModelf4 = new SGread(0);
	sgreadModelf4->mpiCOMM = MPI_COMM_WORLD;
	sgreadModelf4->rank = mpirank;
	sgreadModelf4->size = mpisize;

	SGread* sgreadModelf5 = new SGread(0);
	sgreadModelf5->mpiCOMM = MPI_COMM_WORLD;
	sgreadModelf5->rank = mpirank;
	sgreadModelf5->size = mpisize;

	SGread* sgreadModelf6 = new SGread(0);
	sgreadModelf6->mpiCOMM = MPI_COMM_WORLD;
	sgreadModelf6->rank = mpirank;
	sgreadModelf6->size = mpisize;

	SGinterpolationPoints[0] = sgreadModelf1->read("surplus_sg/f1/surplus.data");
	SGinterpolationPoints[1] = sgreadModelf2->read("surplus_sg/f2/surplus.data");
	SGinterpolationPoints[2] = sgreadModelf3->read("surplus_sg/f3/surplus.data");
	SGinterpolationPoints[3] = sgreadModelf4->read("surplus_sg/f4/surplus.data");
	SGinterpolationPoints[4] = sgreadModelf5->read("surplus_sg/f5/surplus.data");
	SGinterpolationPoints[5] = sgreadModelf6->read("surplus_sg/f6/surplus.data");

	if (mpirank == 0) {
		cout << "Successful" << endl;
	}


	if (mpirank == 0) {
		cout << "+SG Interpolate: ";
	}

	for (int i = 0; i < pointCout; ++i) {
		sgreadModelf1->interpolateValue(&x[i * dim], &interValSG[0 * pointCout * dof + i * dof]);
		sgreadModelf2->interpolateValue(&x[i * dim], &interValSG[1 * pointCout * dof + i * dof]);
		sgreadModelf3->interpolateValue(&x[i * dim], &interValSG[2 * pointCout * dof + i * dof]);
		sgreadModelf4->interpolateValue(&x[i * dim], &interValSG[3 * pointCout * dof + i * dof]);
		sgreadModelf5->interpolateValue(&x[i * dim], &interValSG[4 * pointCout * dof + i * dof]);
		sgreadModelf6->interpolateValue(&x[i * dim], &interValSG[5 * pointCout * dof + i * dof]);
	}

	if (mpirank == 0) {
		cout << "Successful" << endl;
	}

	if (mpirank == 0) {

		// exact
		fstream exact_output;
		exact_output.open(outputFolder + "/exact.csv", ios::out | ios::trunc);
		for (int i = 0; i < pointCout; ++i) {
			for (int f = 0; f < numFunc; ++f) {
				for (int d = 0; d < dof; ++d) {
					exact_output << setprecision(15) << exactVal[f * pointCout * dof + i * dof + d] << "\t";
				}
			}
			exact_output << "\n";
		}
		cout << "Output exact.csv : Successful" << endl;
		exact_output.close();

		// sgInterpolate
		fstream sgInterpolate_output;
		sgInterpolate_output.open(outputFolder + "/sgInterpolate.csv", ios::out | ios::trunc);
		for (int i = 0; i < pointCout; ++i) {
			for (int f = 0; f < numFunc; ++f) {
				for (int d = 0; d < dof; ++d) {
					sgInterpolate_output << setprecision(15) << interValSG[f * pointCout * dof + i * dof + d] << "\t";
				}
			}
			sgInterpolate_output << "\n";
		}
		cout << "Output sgInterpolate.csv : Successful" << endl;
		sgInterpolate_output.close();
	}

	delete[] SGinterpolationPoints;

	delete sgreadModelf1;
	delete sgreadModelf2;
	delete sgreadModelf3;
	delete sgreadModelf4;
	delete sgreadModelf5;
	delete sgreadModelf6;

	delete[] interValSG;
	delete[] exactVal;
}


int main(int argc, char* argv[]) {

	// Problem Parameters
	int dim;
	int dof = 1;
	int pointCout;

	// HMDR Parameters
	int HDMRmaxOrder;
	int HDMRsampleSize;
	double HDMRcutOff;

	// SG Parameters
	int SGmaxLevel;
	double SGcutOff;
	int SGgridType;

	int action;

	MPI_Init (&argc, &argv);
	if (argc < 9) {
		cout << "./main <dim> <pointCout> <HDMRmaxOrder> <HDMRsampleSize> <HDMRcutOff> <SGmaxLevel> <SGcutOff> <SGgridType> <outputFolder> <action>" << endl;
		MPI_Finalize ();
		return 0;
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpisize);

	dim             = atoi(argv[1]);
	pointCout       = atoi(argv[2]);
	HDMRmaxOrder    = atoi(argv[3]);
	HDMRsampleSize  = atoi(argv[4]);
	HDMRcutOff      = atof(argv[5]);
	SGmaxLevel      = atoi(argv[6]);
	SGcutOff        = atof(argv[7]);
	SGgridType      = atoi(argv[8]);
	outputFolder    = argv[9];
	action          = atoi(argv[10]);

	// Allocate constants for test functions
	double sum;
	c = new double*[numFunc];
	w = new double*[numFunc];
	for (int i = 0; i < numFunc; ++i) {
		c[i] = new double[dim]();
		w[i] = new double[dim]();
	}

	x = new double [dim * pointCout]();

	//double normalize;
	double normalize [numFunc];
	normalize[0] = 1.5;
	normalize[1] = dim;
	normalize[2] = 1.85;
	normalize[3] = 7.03;
	normalize[4] = 20.4;

	default_random_engine rd (SEED);
	uniform_real_distribution<double> dist(0.0, 1.0);

	// Write random c[i] and normalize to the normalizing factor
	for (int k = 0; k < numFunc; ++k) {
		sum = 0.0;
		for (int i = 0; i < dim; i++) {
			c[k][i] = dist(rd) ;
			sum += c[k][i];
		}

		for (int i = 0; i < dim; i++) {
			c[k][i] = (c[k][i] / sum) * normalize[k];
		}

		// Write random w[i] and normalize to 1
		sum = 0.0;
		for (int i = 0; i < dim; i++) {
			w[k][i] = dist(rd) ;
			sum += w[k][i];
		}

		for (int i = 0; i < dim; i++) {
			w[k][i] = (w[k][i] / sum) * 1.0;
		}
	}

	// Write random x[i]
	for (int i = 0; i < dim * pointCout; i++) {
		x[i] = dist(rd);
	}


	MPI_Barrier(MPI_COMM_WORLD);


	if (mpirank == 0) {

		fstream index_output;
		index_output.open(outputFolder + "/index.txt", ios::out | ios::trunc);

		index_output << "==================================" << endl;
		index_output << "dim: "            << dim            << endl;
		index_output << "HDMRmaxOrder: "   << HDMRmaxOrder   << endl;
		index_output << "HDMRsampleSize: " << HDMRsampleSize << endl;
		index_output << "HDMRcutOff: "     << HDMRcutOff     << endl;
		index_output << "SGmaxLevel: "     << SGmaxLevel     << endl;
		index_output << "SGgridType: "     << SGgridType     << endl;
		index_output << "SGcutOff: "       << SGcutOff       << endl;
		index_output << "PointCout: "      << pointCout      << endl;
		index_output << "==================================" << endl;
		cout << "Output index.txt : Successful" << endl;
		index_output.close();

		fstream x_output;
		x_output.open(outputFolder + "/x.csv", ios::out | ios::trunc);
		for (int i = 0; i < pointCout; ++i) {
			for (int d = 0; d < dim; ++d) {
				x_output << setprecision(15) << x[i * dim + d] << "\t";
			}
			x_output << "\n";
		}
		cout << "Output x.csv : Successful" << endl;
		x_output.close();

		fstream c_output;
		c_output.open(outputFolder + "/c.csv", ios::out | ios::trunc);
		for (int f = 0; f < numFunc; ++f) {
			for (int d = 0; d < dim; ++d) {
				c_output << setprecision(15) << c[f][d] << "\t";
			}
			c_output << "\n";
		}
		cout << "Output c.csv : Successful" << endl;
		c_output.close();

		fstream w_output;
		w_output.open(outputFolder + "/w.csv", ios::out | ios::trunc);
		for (int f = 0; f < numFunc; ++f) {
			for (int d = 0; d < dim; ++d) {
				w_output << setprecision(15) << w[f][d] << "\t";
			}
			w_output << "\n";
		}
		cout << "Output w.csv : Successful" << endl;
		c_output.close();
	}


	if (action == 1) {
		testWriteSG(dim, dof, SGmaxLevel, SGgridType, SGcutOff);
		testReadSG(dim, dof, pointCout, SGmaxLevel, SGgridType, SGcutOff);
	} else if (action == 2) {
		testWriteHDMR(dim, dof, HDMRmaxOrder, HDMRsampleSize, HDMRcutOff, SGmaxLevel, SGgridType, SGcutOff);
		testReadHDMR(dim, dof, pointCout, HDMRmaxOrder, HDMRsampleSize, HDMRcutOff, SGmaxLevel, SGgridType, SGcutOff);
	} else if (action == 3) {
		testWriteSG(dim, dof, SGmaxLevel, SGgridType, SGcutOff);
		testReadSG(dim, dof, pointCout, SGmaxLevel, SGgridType, SGcutOff);
		testWriteHDMR(dim, dof, HDMRmaxOrder, HDMRsampleSize, HDMRcutOff, SGmaxLevel, SGgridType, SGcutOff);
		testReadHDMR(dim, dof, pointCout, HDMRmaxOrder, HDMRsampleSize, HDMRcutOff, SGmaxLevel, SGgridType, SGcutOff);
	} else {
		cout << "\nAction " <<  action << " Not valid" << endl;
		exit(0);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Clean-up
	for (int i = 0; i < numFunc; ++i) {
		delete[] c[i];
		delete[] w[i];
	}

	delete[] c;
	delete[] w;
	delete[] x;

	MPI_Finalize ();
	return 0;
}





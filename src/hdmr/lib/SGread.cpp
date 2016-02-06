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
	delete[] fin;
	delete[] fout;
	surplusFileName = "";
}


int SGread::read(string fileName) {

	surplusFileName = fileName;

	char fileName_char [fileName.length() + 1];
	strcpy(fileName_char, fileName.c_str());
	LoadData(fileName_char);


	fin = new double[dim];
	fout = new double[TotalDof];

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



void SGread::interpolateDomain(int h,  string fileName) {

	//Mesh density
	double step = 1.0 / (1 << h);

	// Alphabet for permutation (which is mesh steps)
	vector<double> alphabet;
	for (double i = 0.0; i <= 1.0; i += step) {
		alphabet.push_back(i);
	}

	permutation(alphabet, dim);

	return;
}



int SGread::permutation(vector<double> alphabet, int k) {
	vector<int> index(alphabet.size());

	int count = 0;
	calc_permutation(alphabet, k, index, count);

	return count;
}


void SGread::calc_permutation(const vector<double>& alphabet, int k,  vector<int>& index, int& count, int depth) {
	if (depth == k) {
		++count;

		//Construct x vector
		for (auto i = 0; i < k; ++i) {
			fin[i] = alphabet[index[i]];
		}

		Interpolate(fin, fout);

		if (rank == 0) {
			cout <<  fout[0] << "\t";
		}


		std::cout << "\n";
		return;
	}

	for (auto i = 0; i < alphabet.size(); ++i) {
		index[depth] = i;
		calc_permutation(alphabet,  k, index, count, depth + 1);
	}

	return;
}



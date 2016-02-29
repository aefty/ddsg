#ifndef __PRINT_H__
#define __PRINT_H__

namespace std_plus {

  void print_vec(vector<int>& A, string varName = "") {
    cout << varName << ":[ ";
    for (int i = 0; i < A.size() - 1; ++i) {
      cout <<  A[i] << ",";
    }
    cout << A[A.size() - 1] << " ]" << endl;
  }

  void print_vec(vector<double>& A, string varName = "") {
    cout << varName << ":[ ";
    for (int i = 0; i < A.size() - 1; ++i) {
      cout <<  A[i] << ",";
    }
    cout << A[A.size() - 1] << " ]" << endl;
  }


  void print_vec(double* A, int size, string varName = "") {
    cout << varName << ":[ ";
    for (int i = 0; i < size - 1; ++i) {
      cout <<  A[i] << ",";
    }
    cout << A[size - 1] << " ]" << endl;
  }

}
#endif
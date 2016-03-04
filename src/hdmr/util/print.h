#ifndef __PRINT_H__
#define __PRINT_H__

namespace std_plus {




  void print_vec(vector<vector<int>>& A, string varName = "") {
    cout << varName << ":" << endl;
    for (int i = 0; i < A.size(); ++i) {
      for (int j = 0; j < A[i].size(); ++j) {
        cout <<  A[i][j] << " ";
      }
      cout << endl;
    }
    cout << endl;
  }


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


  void print_vec(const double* A, int size, string varName = "") {
    cout << varName << ":[ ";
    for (int i = 0; i < size - 1; ++i) {
      cout <<  A[i] << ",";
    }
    cout << A[size - 1] << " ]" << endl;
  }


  void print_vec(const int* A, int size, string varName = "") {
    cout << varName << ":[ ";
    for (int i = 0; i < size - 1; ++i) {
      cout <<  A[i] << ",";
    }
    cout << A[size - 1] << " ]" << endl;
  }


}
#endif
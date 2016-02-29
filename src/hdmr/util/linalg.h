#ifndef __LINALG_H__
#define __LINALG_H__

namespace std_plus {

  void linalg_dot(vector<double>& A , vector<double>& B, double* rtrn) {
    rtrn[0] = 0.0;

    for (int i = 0; i < A.size(); ++i) {
      rtrn[0] += A[i] * B[i];
    }
  }

  void linalg_sdot(double a, vector<double>& A ,  vector<double>& rtrn) {
    for (int i = 0; i < A.size(); ++i) {
      rtrn[i] =  a * A[i];
    }
  }

  void linalg_less(vector<double>& A,  vector<double>& B) {
    for (int i = 0; i < A.size(); ++i) {
      A[i] -= B[i];
    }
  }


  void linalg_add(double* A, vector<double>& B) {
    for (int i = 0; i < B.size(); ++i) {
      A[i] += B[i];
    }
  }


  void linalg_add(vector<double>& A,  vector<double>& B) {
    for (int i = 0; i < A.size(); ++i) {
      A[i] += B[i];
    }
  }


  void linalg_add(vector<double>& A,  vector<double>& B, vector<double>& rtrn) {
    for (int i = 0; i < A.size(); ++i) {
      rtrn[i] = A[i]  + B[i];
    }
  }

  void linalg_add(double a, vector<double>& A,  double b, vector<double>& B, vector<double>& rtrn) {
    for (int i = 0; i < A.size(); ++i) {
      rtrn[i] = a * A[i]  + b * B[i];
    }
  }
}

#endif
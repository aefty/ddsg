#ifndef __LINALG_H__
#define __LINALG_H__

namespace std_plus {




  double linalg_dot(vector<double>& A , vector<double>& B) {
    double rtrn = 0;
    for (int i = 0; i < A.size(); ++i) {
      rtrn += A[i] * B[i];
    }

    return rtrn;
  }


  double linalg_l2(vector<double>& A) {
    double rtrn = 0;
    for (int i = 0; i < A.size(); ++i) {
      rtrn += A[i] * A[i];
    }
    return pow(rtrn, 0.5);
  }

  double linalg_l2(double* A, int d) {
    double rtrn = 0;
    for (int i = 0; i < d; ++i) {
      rtrn += A[i] * A[i];
    }
    return pow(rtrn, 0.5);
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


  void linalg_less(double* A,  vector<double>& B) {
    for (int i = 0; i < B.size(); ++i) {
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
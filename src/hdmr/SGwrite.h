#ifndef SGwrite_h_is_included
#define SGwrite_h_is_included

#include "SparseGrid.h"

class SGwrite: public AdaptiveSparseGrid {
 public:
  // Properties
  void (*problem)(double*, int, double*);
  int verbose;

  vector<double> maskedX;
  vector<double> xBar;
  vector<int> activeDim;


  // Methods
  SGwrite(void (*problem_)(double*, int, double*), int dim_ , int dof_ , int Lmax_, double epsilon_, int gridType_ , int verbose_ = 0);
  ~SGwrite();


  void setAnchor(vector<double>& xBar);

  int build(vector<int>& activeDim_ );
  int build();

  void write(string surplusFileName);

  void integrateDomain(double* fvalue, double op = 0.0);

  void resetMPI(MPI_Comm mpiCOMM_);

 private:
  // Properties
  const char* hline = "\n========== == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == \n";

  // Methods
  virtual void EvaluateFunctionAtThisPoint( AdaptiveARRAY<double>* x);
};
#endif
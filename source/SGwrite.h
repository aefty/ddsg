#ifndef SGwrite_h_is_included
#define SGwrite_h_is_included

#include "SparseGrid.h"

class SGwrite: public AdaptiveSparseGrid {
 public:
  // Properties
  void (*problem)(double*, int, double*);
  int verbose;

  vector<double> fullPointMap;
  double* fullPoint;

  // Methods
  SGwrite(void (*problem_)(double*, int, double*), int dim_ , int dof_ , int Lmax_, double epsilon_, int gridType_ , int verbose_ = 0);
  ~SGwrite();

  int build( vector<double> fullPointMap_ = {});
  void write( string fileName = "surplus.data");

  void integrateDomain(double* fvalue, double op = 0.0);

 private:
  // Properties
  const char* hline = "\n================================================================================================\n";

  // Methods

  virtual void EvaluateFunctionAtThisPoint( AdaptiveARRAY<double>* x);
};
#endif
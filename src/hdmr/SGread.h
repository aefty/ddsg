#ifndef SGread_h_is_included
#define SGread_h_is_included

#include "SparseGrid.h"

class SGread: public Post {
 public:
  // Properties
  int verbose;
  string surplusFileName = "";


  // Methods
  SGread(int verbose_ = 0);
  ~SGread();

  int read(string surplusFileName);

  void interpolateValue(double* x, double* fvalue);
  void integrateDomain(double* value);

  void resetMPI(MPI_Comm mpiCOMM_);

 private:

};
#endif
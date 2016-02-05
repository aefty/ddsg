#ifndef SGread_h_is_included
#define SGread_h_is_included

#include "SparseGrid.h"

class SGread: public Post {
 public:
  // Properties
  int verbose;
  double* fin;
  double* fout;

  string surplusFileName = "";


  // Methods
  SGread(int verbose_ = 0);
  ~SGread();

  int read(string fileName = "surplus.data");

  void interpolateDomain(int h,  string fileName = "interpolate.data");
  void interpolateValue(double* x, double* fvalue, double op = 0.0);

  void integrateDomain(double* value, double op = 0.0);


 private:
  // Properties
  const char* hline = "\n================================================================================================\n";

  // Methods

  int permutation(vector<double> alphabet, int k);
  void calc_permutation(const vector<double>& alphabet, int k,  vector<int>& index, int& count, int depth = 0);



};
#endif
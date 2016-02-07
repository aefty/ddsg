/**
 * File has been edited by Aryan E. in sections noted by the tag  #AE:
 */


#ifndef Post_h_IS_INCLUDED
#define Post_h_IS_INCLUDED

#include "../DataStructure/DataStructure.h"
#include <fstream>

class Post {
 public:

  Post();
  ~Post() {};

  virtual void LoadData(char* filename);

  virtual void Interpolate(Array<double>& x, double* value);

  virtual void Interpolate(double* x, double* value);

  virtual void Integrate(double* value);

  virtual double IndextoCoordinate(int i, int j);


  //#AE: function that will call other basis functions based on gridType
  virtual double BasisFunction(double x, int i  , int j);

  virtual double LinearBasis(double x, int i, int j);
  virtual double LagrangePoly(double x, int i, int j);
  virtual void   GetChebyshevPoints(int i, Array<double>& points);

  //#AE: Flipup basis function based on Pfluger
  virtual double FlipUpBasis(double x, int i, int j);
  //#AE: QuadPoly
  virtual double PolyBasis(double x, int i, int j);


  virtual double LinearBasisVolumeIntegral(int level);
  virtual void   PlotSparseGrid(char *filename = "grid.plt");

  //! \param Number of dimensions
  int dim ;
  //! \param Total number of points
  int nno;
  //! \param Total degree of freedom
  int TotalDof;
  //! \param interpolation levle
  int Level;

  //#AE: param grid Type, defualt to 1
  int gridType = 1;

  //! Define MPI rank
  int rank;

  // ! Define MPI size
  int size;

  //! Define MPI Comm Scope
  MPI_Comm mpiCOMM = MPI_COMM_WORLD;

  //! Matrix to store the multi-index
  Matrix1<int> index;
  //! Matrix to store teh surplus
  Matrix1<double> surplus;

};
#endif
/*Class:Post

NAME:  Post - This is the class for offline serial interpolation.

DESCRIPTION:

  This class is used for offline interpolation.

  You need to provide a surplus file in the following format

  //! dimension,  number of point,  total degree of freedom,
  //! i_1..... i_dim  j_1 .... j_dim surplus_1 .... surplus_TotalDof


CONSTRUCTORS AND INITIALIZATION:

    The constructor takes no arguments.

    You need to load the surplus file to initialize the code.


MEMBER FUNCTIONS:


  Most member functions are self-explanatory.


End:
*/
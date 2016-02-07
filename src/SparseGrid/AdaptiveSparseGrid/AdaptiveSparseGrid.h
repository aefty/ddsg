/*
Filename : AdaptiveSparseGrid.h
Prepared by : Xiang Ma and Nicholas Zabaras
Version 1.0:  March 4,  2008
Version 2.0:  March 19, 2008
Version 5.0:  July  30, 2009
*/

#ifndef AdaptiveSparseGrid_h_is_included
#define AdaptiveSparseGrid_h_is_included

#include "../AdaptiveDataStructure/AdaptiveDataStructure.h"
#include "../BlockAllocator/BlockAllocator.h"
#include <set>
#include <deque>
#include <fstream>



class AdaptiveSparseGrid {

  public:

	AdaptiveSparseGrid(int dim, int Lmax, double epsilon = 0.0, int type = 0);

	//! A virtual deconstructor
	virtual ~AdaptiveSparseGrid();

	//! deque is a container from the STL to store the information of all of the points
	//! in the sparse grid
	deque<AdaptiveData*> SparseGrid;

	//! defines three BlockAllocators to allocate the memory for AdaptiveData and Multi-index
	BlockAllocator <AdaptiveData> AdaptiveDataAllocator;
	BlockAllocator <AdaptiveNodeData> AdaptiveNodeDataAllocator;
	BlockAllocator <AdaptiveARRAY<int> > AdaptiveCoordAllocator;

	//! set is a associative container to store all of the active points during adaptive refinment.
	//! The use of set is to ensure a unique interpolation point is genrated.
	//! The comparison of the coordinate is defined through the struct object AdaptiveIndexCompare
	//! from class AdaptiveDataStructure
	typedef set<AdaptiveData*, AdaptiveDataCompare> AdaptiveGrid;
	AdaptiveGrid ActiveIndex;

	virtual void   BuildAdaptiveSparseGrid(int print_ = 1, int restart_ = 0, char* filename = "surplus.plt");

	virtual void   ConstructAdaptiveSparseGrid();

	virtual void   Restart(char* filename);

	virtual void   BuildFirstLevel();

	virtual void   SetGridType(int type) {gridType = type;};

	//Subroutines for interpolation
	virtual void   SpInterpolate(AdaptiveARRAY<double> *x, double* value);
	virtual void   SpInterpolate(double* x, double* value);
	virtual void   SpInterpolateLevel(Matrix1<double>& px, double *temp);

	//!Subroutines for linear and polynomial basis functions
	virtual double BasisFunction(double x, int i, int j );
	virtual double LinearBasis(double x, int i, int j);
	virtual double FlipUpBasis(double x, int i, int j);
	virtual double LagrangePoly(double x, int i, int j);
	virtual void   GetChebyshevPoints(int i, Array<double>& points);

	//#AE add poly basis function
	virtual double PolyBasis(double x, int i, int j);


	virtual void   Refine( AdaptiveARRAY<int> *, AdaptiveARRAY<int> *);

	virtual void   FindOrInsertActiveIndex(AdaptiveARRAY<int> *i, AdaptiveARRAY<int> *j1, AdaptiveARRAY<int> *j2 );

	//!Subroutines for integration
	virtual void   SpIntegrate(double* value);
	virtual double LinearBasisVolumeIntegral(int level);
	virtual void   ComputeMeanAndVar(char *filename = "result.plt");
	virtual void   StoreMeanAndVar(double* value, char *filename) {};


	virtual void   PlotSparseGrid(char *filename = "grid.plt");

	virtual void   StoreSurplus(AdaptiveARRAY<int>& index, char *filename = "surplus.plt" );

	virtual void   Cleanup();

	//! Function to define when to refine the points; You may overload this function.
	virtual int    IsRefine(double* value, AdaptiveARRAY<int> *, AdaptiveARRAY<int> *);

	//! Functions you must overload
	virtual void   EvaluateFunctionAtThisPoint( AdaptiveARRAY<double> *x) = 0;


	//! Set the actions before and after storing the hierarchical surplused
	virtual void   BeforeStoreSurplus() {};
	virtual void   AfterStoreSurplus() {};


	//!Warning: For polynomial integration, its purpose is only to test the code, since it is very slow.
	//!Computing the intergral of the polynomial basis function using Gauss-Legendre quadrature.
	virtual void   calculate_wts_and_quad_pts(const double low, const double up);
	virtual double PolynomialVolumeIntegral(int i, int j);
	double *wts; ///< Weights for numerical integration  (1D polynomial)
	double *abs; ///< Abscissa for numerical integration (1D polynomial)
	int no_quad_pts;

	virtual int NumberOfPoints(); ///< Number of total collocation points

	virtual int NumberOfPointsInThisProc(); ///< Number of total collocation points in this process

	virtual int NumberOfActivePoints();     ///< Number of active collocation points at each level

	virtual int InterpolationLevel() {return L - 1;};

	virtual double IndextoCoordinate(int i, int j);

	virtual int    CoordinateToIndex(double x);


	Array<double> currentvalue;   ///< integration value

	virtual void UpdateError();   ///< Update the error at this level

	//! \param Current processor ID
	int rank;

	//! \param Number of precessors
	int size;

	//! \param Maximum interpolation level
	int Lmax;

	//! \param Number of dimensions
	int dim;

	//! \param Error threshold
	double epsilon;

	//! \param type of error indicator
	int type;

	//! \param current level of sparse grid
	int L;

	//! \param number of degree to store
	int NumberToStore;

	//! \param degree of freedom per node
	int DofPerNode;

	//! \param total number of degree  = NumberToStore * DofPerNode
	int TotalDof;

	//! \param  grid type 1:Clenshaw-Curtis 2: Chebyshev-Gauss-Lobatto 3:other
	int gridType;

	//! Array to store surplus
	double *surplus;

	//! \param restart
	int restart;

	//! \param whether to print intermediate results
	int print;

	//! Define MPI Comm Scope
	MPI_Comm mpiCOMM = MPI_COMM_WORLD;

};
#endif
/*Class:AdaptiveSparseGrid

NAME:  AdaptiveSparseGrid - Main class for constructing the adaptive sparse grid interpolation

DESCRIPTION:

  For details of algorithm, definition of the basis funciton,
  the reader is referred to the following papers:

  <b> Xiang Ma and N. Zabaras </b>
  An adaptive hierarchical sparse grid collocation algorithm for the solution
  of the stochastic differential equations, Journal of Computational Physics, Vol. 228, pp. 3084-3113, 2009.

  The interpolation space is the hypercube [0,1]^D. For other space, it must be transformed
  to this hypercube.

  This toolbox is fully parallel. You will need  the mpi library to compile
  the source code.

  Version 1.0: The construction of sparse grid in this code is through the multi-index (i,j) of the interpolation point.
  The storage of the points is distributed among the processors. The performance is not very optimal due
  to the communication between processors. However, if we want to solve large-scale problem, this code is useful
  to enlarge the memory storage.

  In this version, you can switch between two grids: Newton-Cotes and Clenshaw-Curtis. However, the Clenshaw-Curtis
  is only for test purpose whose speed is very slow. In addition, there is no adaptivity for Clenshaw-Curtis grid.
  Therefore, it is always recommending Newton-Cotes grid which is the default grid.

  Version 2.0: Resolve the memory problem when storing the surplus. Add Restart function to enable the capability to
  restart the construction for surplus file. For function "StoreSurplus", we store one more parameter: current interpolation
  level for the compatibility with function "Restart".

  Version 3.0: Modify the part which compute the hierarchical surplus. Remove the MPI_Bcast to send adaptivity information.
  Now, we use an arrary to store the information, which only involves one MPI_Allreduce operation. Now, it is nearly two times
  faster than last version.

  Version 4.0: Modify the datastructure part which is now compatible with the new g++ compiler. Add more functionality to
  the class Post.

  Version 5.0: Add a new error indicator which depends on the integration value. Rewrite some example codes.

  For stochastic modeling, we need to construct not only the surplus for the function itself but also the surplus
  for the square of the function in order to calculate the mean and variance.

   <b> WARNING </b>
  This the base class for the tool box. You must overload the function "EvaluateFunctionAtThisPoint".


CONSTRUCTORS AND INITIALIZATION:

    \param dim     The dimension of the function
	\param Lmax    The allowed maximum interpolation level
	\param epsilon The error criteria to control the adaptivity. By using
	               epsilon == 0, it goes back to the conventional sparse grid

    \param type    The different types of error adaptivity indicators. 0 = Maximum surplus; 1 = Integration value
                   The default is the maximum suplus. The type 0 error is suitable if you need accurate interpolation value. On the
				   other hand, if you only cares the integration value, then type 1 error is more appropriate.
    optional:
	\param gridtype The default is Newton-Cotes

	<b> WARNING </b>
	Remember to initilize the following variables in your class:
	//! \param number of degree to store
	int NumberToStore;
	//! \param degree of freedom per node
	int DofPerNode;
	//! \param total number of degree  = NumberToStore * DofPerNode
	int TotalDof;
	//! Array to store surplus
    double *surplus;

MEMBER FUNCTIONS:


  "ConstructAdaptiveSparseGrid" -  Main rountine to construct the adaptive sparse grid.

  "BuildFirstLevel"             -  We always refine the first level,which is just a single point.

  "SpInterpolate"               -  For the input interpolation point AdativeARRAY<double> *x,
                                   return the interpolating value using all of the current points in the sparse grid.

  "BuildAdaptiveSparseGrid"     -  Driver to construct the sparse grid

  "Restart"                     -  Load surplus from the file and initilize the sparse grid.

  "Refine"                      -  Generate the new interpolation points from next level
                                   if the adaptive criteria is satisfied.

  "IsRefine"                    -  Define your own adaptivity criteria.

  "FindOrInsertActiveIndex"     -  Check if the new point has alreay existed.
                                   If not, insert the new point into the active index.

  "ComputeMeanAndVar"           -  Compute mean and variance for the solution. When using this function,
                                   remember to provide not only the surplus for the function itself but
								   also the surplus for the square of the function.

  "StoreMeanAndVar "            -  This function is not implemented. You need to overload it.

  "VolumeIntegral"              -  Integrate the basis function at each interpolation point.

  "PlotSparseGrid"              -  Plot the sparse grid in Tecplot format.

  "StoreSurplus"                -  Store the surplus for offline processing.

  "NumberOfPoints"              -  Return the number of points in the current adaptive sparse grid.

  "InterpolationLevel"          -  Return the current interpolation level.

  "CoordinateToIndex"           -  Transform the coordinate to the corresponding multi-index i.

  "IndextoCoordinate"           -  Transform the multi-index (i,j) to coordinate of the interpolation point.

  "Cleanup"                     -  Clean up the memory space.


COPYRIGHT:

        Professor Nicholas Zabaras
        Materials Process Design & Control Laboratory
        Sibley School of Mechanical and Aerospace Engineering
        101 Frank H. T. Rhodes Hall
        Cornell University
        Ithaca, NY 14853-3801

        Email: zabaras@cornell.edu
        Phone: 607 255 9104
        Fax:   607 255 1222


AUTHOR:

  Xiang Ma, xm25@cornell.edu

End:
*/

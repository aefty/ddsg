/**
 * File has been edited by Aryan E. in sections noted by the tag  #AE:
 */

#include "AdaptiveSparseGrid.h"
#define M_PI  3.14159265358979323846

AdaptiveSparseGrid::AdaptiveSparseGrid(int d, int L, double eps, int gridType_ ) {
	dim = d;
	Lmax = L;
	epsilon = eps;
	//type = type_;

	surplus = NULL;

	//#AE: type_
	gridType = gridType_;

	SetGridType(gridType);
}

AdaptiveSparseGrid::~AdaptiveSparseGrid() {

	if (surplus) {
		delete[] surplus;
		surplus = NULL;
	}

	Cleanup();
}

//! Clean up the memory
void AdaptiveSparseGrid::Cleanup() {

	SparseGrid.clear();
	ActiveIndex.clear();
	AdaptiveCoordAllocator.ClearAll();
	AdaptiveNodeDataAllocator.ClearAll();
	AdaptiveDataAllocator.ClearAll();
}

void AdaptiveSparseGrid::BuildAdaptiveSparseGrid(int print_, int restart_, char* filename) {
	print = print_;
	restart = restart_;

	if ( ! restart) {
		//First Construct the first level
		BuildFirstLevel();
		//! Set the current interpolation level to 1
		L = 1;
	}

	// else load surplus from the file
	else {
		Restart(filename);
	}

	ConstructAdaptiveSparseGrid();
}

void AdaptiveSparseGrid::ConstructAdaptiveSparseGrid() {
	//! Define a buffer deque to store all of the intermediate information
	deque<AdaptiveData*> buffer;

	//! iterator for transverse the active index
	set<AdaptiveData*, AdaptiveDataCompare>::iterator it1;
	set<AdaptiveNodeData*, AdaptiveNodeDataCompare>::iterator it2;
	AdaptiveGrid OldIndex;

	AdaptiveData* pData ;
	AdaptiveNodeData* pNodeData;

	//! Here we still keep a copy of active index for each processor for the purpose
	//! of the easy implementation of the adaptivity

	while ( !ActiveIndex.empty() && L <= Lmax ) {
		int gnumber = NumberOfActivePoints();


		if (rank == 0 && print) {
			cout << "Now, it is in Level: " << L << endl;
			cout << "The active number of points are: " << gnumber << endl;
		}

		//Copy the ActiveIndex to an oldIndex
		OldIndex = ActiveIndex;

		ActiveIndex.clear();

		//! Define a temporary matrix to store all of the active point for parallel computation
		Matrix1<double> gpoint;
		gpoint.redim(gnumber, dim);
		double time1, time;
		time1 = MPI_Wtime();
		int row = 1;

		//! Extract the coordiante information
		for (it1 = OldIndex.begin(); it1 != OldIndex.end(); ++it1) {
			for (it2 = (*it1)->NodeData.begin(); it2 != (*it1)->NodeData.end(); ++it2) {
				for ( int i = 1; i <= dim ; i++) {
					gpoint(row, i)  = IndextoCoordinate((*(*it1)->index)(i), (*(*it2)->index)(i));
				}
				row++;
			}
		}

		// Parallel Implementation

		int node;

		AdaptiveARRAY<double> px;
		px.redim(dim);

		//! Distribute the points among the processors
		int numNodesPerProcessor = 0;

		for (int i = 0; i < gnumber; i++)
			if ( (i % size) == rank)
			{ numNodesPerProcessor++; }

		int colnumber = TotalDof * numNodesPerProcessor;
		//!Local value to store function value
		double*  local_value = new double[colnumber];
		//! An array to store the nodes which belong to this processor
		int*     mynodes = new int[numNodesPerProcessor];

		for ( int no = 1 ; no <= numNodesPerProcessor; no++) {
			node = rank + size * (no - 1) + 1;
			mynodes[no - 1] = node;

			for ( int i = 1 ; i <= dim ; i++)
			{ px(i) = gpoint(node, i); }

			EvaluateFunctionAtThisPoint(&px);

			for ( int i = 0; i < TotalDof ; i++) {
				local_value[(no - 1)*TotalDof + i]  = surplus[i] ;
			}
		}

		/**
		 * EDIT BY ARYAN
		 */
		MPI_Barrier(mpiCOMM);

		time = MPI_Wtime();

		if (rank == 0 && print ) { cout << "Parallel Calculation using " << time - time1 << endl; }

		gpoint.cleanup();

		//! Set the action before storing the surplus
		BeforeStoreSurplus();

		int index = 1;
		int cnt = 0, mynode;

		//! Extract the coordiante information
		for (it1 = OldIndex.begin(); it1 != OldIndex.end(); ++it1) {

			int number = (*it1)->NodeData.size();

			//! Define a temporary matrix to store all of the active point for parallel computation
			Matrix1<double> point;
			point.redim(number, dim);

			int row = 1;

			for (it2 = (*it1)->NodeData.begin(); it2 != (*it1)->NodeData.end(); ++it2) {
				for ( int i = 1; i <= dim ; i++) {
					point(row, i)  = IndextoCoordinate((*(*it1)->index)(i), (*(*it2)->index)(i));
				}

				row++;
			}


			//! If there are points in this processor, generate a new data to store information
			if ( numNodesPerProcessor != 0) {
				pData = AdaptiveDataAllocator.NewItem();
				pData->index = (*it1)->index;
				buffer.push_front(pData);
			}


			double* temp = new double[TotalDof * number];
			//! Calculate the hierarchical surplus
			SpInterpolateLevel(point, temp);

			//! Calculate the hierarchial surplus and generate the adaptivity information
			Array<int> local_refine; local_refine.redim(number); local_refine.fill(0);
			Array<int> global_refine; global_refine.redim(number);

			int temp_index = index;
			int temp_cnt = cnt;

			it2 = (*it1)->NodeData.begin();

			for ( int n = 1;  n <= number; n++) {

				if ( (temp_cnt != numNodesPerProcessor) && (temp_index == mynodes[temp_cnt]) ) {
					for ( int i = 0; i < TotalDof; i++) {
						local_value[temp_cnt * TotalDof + i] = surplus[i] = local_value[temp_cnt * TotalDof + i] - temp[(n - 1) * TotalDof + i];
					}

					//! Check the adaptivity criteria
					local_refine(n) = IsRefine(surplus, (*it1)->index, (*it2)->index);
					temp_cnt ++;
				}

				temp_index++;
				++it2;
			}

			//! Free space
			delete[] temp;


			int error = MPI_Allreduce(local_refine.pData, global_refine.pData, number, MPI_INT, MPI_SUM, mpiCOMM);
			//if (error) { cout << "ERORR 221 : " << error <<  "|" << rank << endl; }

			//! Free space
			local_refine.cleanup();

			it2 = (*it1)->NodeData.begin();

			// calculate herarchial surplus for each coordinate
			for ( int n = 1 ; n <= number; n++) {

				//! If this point belongs to this processor
				if ( (cnt != numNodesPerProcessor) && (index == mynodes[cnt]) ) {
					//Copy data to sparse grid
					pNodeData = (*it2);
					pNodeData->surplus = new double[TotalDof];

					for ( int i = 0; i < TotalDof; i++) {
						pNodeData->surplus[i] = local_value[cnt * TotalDof + i];
					}

					//! Refinement
					if ( global_refine(n))
					{ Refine((*it1)->index, (*it2)->index); }

					cnt ++;
					//! Insert the this nodadata to the buffer
					pData->NodeData.insert(pNodeData);
				}

				//! the point doesn't belong to this processor
				else {
					if ( global_refine(n))
					{ Refine((*it1)->index, (*it2)->index); }

					//! If this point doesn't belong to this processor, delete the point.
					AdaptiveCoordAllocator.DeleteItem((*it2)->index);
					AdaptiveNodeDataAllocator.DeleteItem(*it2);
				}

				++it2;
				index++;

			}

			//free space
			global_refine.cleanup();

			//! If there are no point for this processor, delete the data.
			if ( numNodesPerProcessor == 0) {
				AdaptiveCoordAllocator.DeleteItem((*it1)->index);
				AdaptiveDataAllocator.DeleteItem(*it1);
			}
		}// End for

		//free space
		delete[] local_value;
		delete[] mynodes;

		/**
		* EDIT BY ARYAN
		*/
		MPI_Barrier(mpiCOMM);

		time1 = MPI_Wtime();

		if (rank == 0 && print) { cout << "Surplus  Calculation using " << time1 - time << endl; }

		//! Insert all of the active points to the adaptive spare grid.
		SparseGrid.insert(SparseGrid.begin(), buffer.begin(), buffer.end());

		if ( type == 1)
		{ UpdateError(); }

		//free space
		buffer.clear();
		OldIndex.clear();
		L += 1;

		//! Set the action before storing the surplus
		{ AfterStoreSurplus(); }
	}//End for
}

void AdaptiveSparseGrid::UpdateError() {

	currentvalue.redim(TotalDof);
	SpIntegrate(currentvalue.pData);
}

int AdaptiveSparseGrid::IsRefine(double* value, AdaptiveARRAY<int>* ix, AdaptiveARRAY<int>* iy) {
	if ( type == 1) {

		double mul = 1.0;

		for ( int j = 1; j <= dim; j++) {

			mul *= LinearBasisVolumeIntegral((*ix)(j));

		}

		double sum1 = 0.0, sum2 = 0.0;

		for ( int j = 0; j < TotalDof; j++) {
			if ( currentvalue(j + 1) != 0) {
				sum1 += pow(mul * value[j], 2.0);
				sum2 += pow(currentvalue(j + 1), 2.0);
			}
		}

		double error;

		if ( sum2 <= 1e-15 )
		{ return 1; }

		else {
			error = sqrt(sum1) / sqrt(sum2);

			if ( error >= epsilon ) { return 1; }

			else { return 0; }
		}

	} else {
		double maximum = 0.0;

		for ( int i = 0 ; i < TotalDof; i++) {
			if ( fabs(value[i]) > maximum )  { maximum = fabs(value[i]); }
		}

		if ( maximum >= epsilon ) { return 1; }

		else { return 0; }
	}
}

void AdaptiveSparseGrid::BuildFirstLevel() {
	AdaptiveARRAY<int>* px1 = AdaptiveCoordAllocator.NewItem(dim);
	AdaptiveARRAY<int>* px2 = AdaptiveCoordAllocator.NewItem(dim);


	px1 -> fill(1);
	px2 -> fill(0);

	//! The first level is always a single point.
	AdaptiveARRAY<double> x;
	x.redim(dim);
	x.fill(0.5);

	L = 0;

	//! For the fisrt point, we always assign it to rank 0
	if (rank == 0) {
		EvaluateFunctionAtThisPoint(&x);
		AdaptiveData* pData = AdaptiveDataAllocator.NewItem();
		AdaptiveNodeData* pNodeData = AdaptiveNodeDataAllocator.NewItem();
		pData->index = px1;
		pNodeData->index = px2;
		pNodeData->surplus = new double[TotalDof];

		for ( int i = 0; i < TotalDof; i++) {
			pNodeData->surplus[i] = surplus[i];
		}

		pData->NodeData.insert(pNodeData);
		SparseGrid.push_front(pData);
	}

	MPI_Barrier(mpiCOMM);

	//Generate New sons for each dimension
	Refine(px1, px2);

	// Compute the integration value
	if ( type == 1)
	{ UpdateError(); }

	if ( rank != 0) {
		AdaptiveCoordAllocator.DeleteItem(px1);
		AdaptiveCoordAllocator.DeleteItem(px2);
	}
}

void AdaptiveSparseGrid::Refine(AdaptiveARRAY<int>* i, AdaptiveARRAY<int>* j) {
	AdaptiveARRAY<int>* px, *px1, *px2;

	//Refine for each direction
	for ( int d = 1; d <= dim; d++) {

		px = AdaptiveCoordAllocator.NewItem(dim);
		px1 = AdaptiveCoordAllocator.NewItem(dim);
		px2 = AdaptiveCoordAllocator.NewItem(dim);

		for ( int k = 1; k <= dim; k++) {
			(*px1)(k) = (*px2)(k)  = (*j)(k);
			(*px)(k) = (*i)(k);
		}

		(*px)(d) += 1;


		//#AE: Build Refinement for CC gridType
		if (gridType == 1) {

			if ((*i)(d) == 1) {
				(*px1)(d) = 1;
				(*px2)(d) = 3;
				FindOrInsertActiveIndex(px, px1, px2);

				//!Special treatment. For the points in level two, there is only one son.

			} else if ( (*i)(d) == 2 ) {

				if ( (*j)(d) == 1 ) {
					(*px1)(d) = 2;

				} else {
					(*px1)(d) = 4;
				}

				(*px2)(d) = (*px1)(d);
				FindOrInsertActiveIndex(px, px1, px2);

			} else {//!Generate neigborhood points
				(*px1)(d) = 2 * ((*j)(d));
				(*px2)(d) = 2 * (*j)(d) - 2;

				FindOrInsertActiveIndex(px, px1, px2);
			}

			//#AE: Build Refinement for each Pfluger gridType
		} else if (gridType == 3) {

			if ((*i)(d) == 1) {
				(*px1)(d) = 2;
				(*px2)(d) = 4;
				FindOrInsertActiveIndex(px, px1, px2);

			} else {
				(*px1)(d) = 2 * ((*j)(d));
				(*px2)(d) = 2 * (*j)(d) - 2;

				FindOrInsertActiveIndex(px, px1, px2);
			}
		} else if (gridType == 4) { // AE Polybase

			if ((*i)(d) == 1) {
				(*px1)(d) = 2;
				(*px2)(d) = 4;
				FindOrInsertActiveIndex(px, px1, px2);

			} else {
				(*px1)(d) = 2 * ((*j)(d));
				(*px2)(d) = 2 * (*j)(d) - 2;

				FindOrInsertActiveIndex(px, px1, px2);
			}
		}
	}
}

void AdaptiveSparseGrid::FindOrInsertActiveIndex(AdaptiveARRAY<int>* i, AdaptiveARRAY<int>* j1, AdaptiveARRAY<int>* j2) {
	pair<set<AdaptiveData*, AdaptiveDataCompare>::iterator, bool> plt;

	AdaptiveData* pData = AdaptiveDataAllocator.NewItem();
	AdaptiveNodeData* pNodeData1 = AdaptiveNodeDataAllocator.NewItem();
	AdaptiveNodeData* pNodeData2 = AdaptiveNodeDataAllocator.NewItem();
	pData->index = i;
	pNodeData1->index = j1;
	pNodeData2->index = j2;

	//! If the multi-index i is already existed, plt.second = false
	plt = ActiveIndex.insert(pData);

	// remove auxiliary key
	if ( plt.second == false ) {

		AdaptiveDataAllocator.DeleteItem(pData);
		AdaptiveCoordAllocator.DeleteItem(i);

		pair<set<AdaptiveNodeData*, AdaptiveNodeDataCompare>::iterator, bool> plt1;

		//! If the multi-index j is already existed, plt.second = false
		plt1 = (*(plt.first))->NodeData.insert(pNodeData1);

		if ( plt1.second == false ) {
			AdaptiveCoordAllocator.DeleteItem(j1);
			AdaptiveNodeDataAllocator.DeleteItem(pNodeData1);
		}

		//! If the multi-index j is already existed, plt.second = false
		plt1 = (*(plt.first))->NodeData.insert(pNodeData2);

		if ( plt1.second == false ) {
			AdaptiveCoordAllocator.DeleteItem(j2);
			AdaptiveNodeDataAllocator.DeleteItem(pNodeData2);
		}
	}

	//! If the multi-index i is not existed, insert the multi-index j directly
	else if (plt.second == true) {
		pData->NodeData.insert(pNodeData1);
		pData->NodeData.insert(pNodeData2);
	}
}

//! Interpolate the function value in the current sparse grid point
void AdaptiveSparseGrid::SpInterpolate(AdaptiveARRAY<double>* px, double* gvalue) {
	double   temp;
	deque<AdaptiveData*>::reverse_iterator it1;
	set<AdaptiveNodeData*, AdaptiveNodeDataCompare>::iterator it2;

	double* local_temp = new double[TotalDof];

	for ( int i = 0 ; i < TotalDof; i++)
	{ local_temp[i] = 0.0; }


	for (it1 = SparseGrid.rbegin(); it1 != SparseGrid.rend(); ++it1) {

		for (it2 = (*it1)->NodeData.begin(); it2 != (*it1)->NodeData.end(); ++it2) {

			temp = 1.0;

			for ( int j = 1; j <= dim; j++) {

				temp *= BasisFunction((*px)(j), (*(*it1)->index)(j), (*(*it2)->index)(j));

				if ( temp == 0.0) { break; }
			}

			for ( int i = 0; i < TotalDof; i++)
			{ local_temp[i] += temp * (((*it2)->surplus)[i]); }
		}
	}


	int  error  = MPI_Allreduce(local_temp, gvalue, TotalDof, MPI_DOUBLE, MPI_SUM, mpiCOMM);
	if (error) { cout << "ERORR 538 : " << error << "|"  << rank <<  endl; }


	delete[] local_temp;
}

//! Interpolate the function value in the current sparse grid point for parallel implementation purpose
void AdaptiveSparseGrid::SpInterpolateLevel(Matrix1<double>& px, double* gvalue) {
	double   temp;
	deque<AdaptiveData*>::reverse_iterator it1;
	set<AdaptiveNodeData*, AdaptiveNodeDataCompare>::iterator it2;

	int number = TotalDof * px.nx;

	double* local_temp = new double[number];

	for ( int i = 0 ; i < number; i++)
	{ local_temp[i] = 0.0; }


	for (it1 = SparseGrid.rbegin(); it1 != SparseGrid.rend(); ++it1) {

		for (it2 = (*it1)->NodeData.begin(); it2 != (*it1)->NodeData.end(); ++it2) {

			for ( int i = 1; i <= px.nx; i++) {
				temp = 1.0;

				for ( int j = 1; j <= dim; j++) {

					temp *= BasisFunction(px(i, j), (*(*it1)->index)(j), (*(*it2)->index)(j));

					if ( temp == 0.0) {
						break;
					}
				}

				for ( int j = 0; j < TotalDof; j++) {
					local_temp[(i - 1)*TotalDof + j] += temp * (((*it2)->surplus)[j]);
				}
			}
		}
	}

	int error = MPI_Allreduce(local_temp, gvalue, number, MPI_DOUBLE, MPI_SUM, mpiCOMM);
	//if (error) { cout << "ERORR 583 : " << error <<  "|" << rank << endl; }

	delete[] local_temp;
}

void AdaptiveSparseGrid::SpInterpolate(double* x, double* value) {
	AdaptiveARRAY<double> px;
	px.redim(dim);

	for ( int i = 1; i <= dim; i++) {
		px(i) = x[i - 1];
	}

	SpInterpolate(&px, value);
}

//! Transform the multi-index(i,j) to the coordinate of interpolation point.
double AdaptiveSparseGrid::IndextoCoordinate(int i, int j ) {
	double m = pow(2.0, (i - 1)) + 1;

	if ( i == 1) {
		return 0.5;
	}

	//#AE: return points based on  grid types
	switch (gridType) {
		case 1:
			return (j - 1.0) / (m - 1);
			break;

		case 2:
			return (-cos(M_PI * (j - 1) / (m - 1)) + 1) / 2.0;
			break;

		case 3:
			m = pow(2.0, i) + 1.0;
			return (j - 1) / (m - 1.0);
			break;
		case 4:
			m = pow(2.0, i) + 1;
			return (j - 1) / (m - 1.0) ;
			break;

		//#AE: throw error when grid type not defined
		default:
			cout << "!! gridtype not valid ( " << gridType << " ) !!" << endl;
			exit(0);

	}
}


//! Transform coordinate of interpolation point to the multi-index i.
int AdaptiveSparseGrid::CoordinateToIndex(double x) {

	if ( x == 0.5 ) {
		return 1;

	} else if ( (x == 0.0) || (x == 1.0) ) {
		return  2;

	} else {
		for ( int k = 3; k <= (Lmax + 1); k++) {
			double m = x * pow(2.0, k - 1);

			if ( fabs((m - floor(m))) <= 1e-20) {
				return k;
				break;
			}
		}
	}
}

double AdaptiveSparseGrid::BasisFunction(double x, int i, int j) {


	//#AE: Select basis function based gridType
	switch (gridType) {
		case 1:
			return LinearBasis(x, i, j);
			break;

		case 2:
			return LagrangePoly(x, i, j);
			break;

		case 3:
			return FlipUpBasis(x, i, j);
			break;

		case 4: // Basis
			return PolyBasis(x, i, j);
			break;


		//#AE: throw error when grid type not defined
		default:
			cout << "!! gridtype not valid ( " << gridType << " ) !!" << endl;
			exit(0);
	}
}

//! Linear basis function
double AdaptiveSparseGrid::LinearBasis(double x, int i, int j) {

	//First Level
	if (  i == 1 ) {
		return 1.0;
	} else {
		double m = pow(2.0, (i - 1));
		double xp = IndextoCoordinate(i, j );

		if ( fabs( x - xp ) >= (1.0 / m)) {
			return 0.0;
		} else {
			return (1 - m * fabs(x - xp));
		}
	}
}


/**
 * #AE
 * Polynomial  Function
 * @param  x [X val]
 * @param  i [Level Depth]
 * @param  j [Basis Function Index]
 * @return   [Value]
 */
double AdaptiveSparseGrid::PolyBasis(double x, int i, int j) {

	if (i < 3) {
		return FlipUpBasis(x, i, j);
	} else {

		double m = pow(2.0, i);
		double xp = IndextoCoordinate(i, j);

		// Wings
		if ( x <= (1.0 / m) && xp == 1. / m)  {
			return (-1.0 * m * x + 2.0);
		} else if ( x >= (1.0 - 1.0 / m )  && xp == (1.0 - 1. / m) )  {
			return (1.0 * m * x + (2.0 - m));
		} else {

			// Body
			double x1 = xp - 1.0 / m;
			double x2 = xp + 1.0 / m;
			double temp = (x - x1) * (x - x2) / ((xp - x1) * (xp - x2));

			if (temp > 0) {
				return temp;
			} else {
				return 0.0;
			}
		}
	}
}


/**
 * #AE
 * FlipUpBasis Function
 * @param  x [X val]
 * @param  i [Level Depth]
 * @param  j [Basis Function Index]
 * @return   [Value]
 */
double AdaptiveSparseGrid::FlipUpBasis(double x, int i, int j) {

	double xp, m;
	if (  i == 1 ) {
		return 1.0;
	} else {
		double m = pow(2.0, i);
		double xp = IndextoCoordinate(i, j);

		// Wing
		if ( x <= (1.0 / m) && xp == 1. / m)  {
			return (-1.0 * m * x + 2.0);
		} else if ( x >= (1.0 - 1.0 / m )  && xp == (1.0 - 1. / m) )  {
			return (1.0 * m * x + (2.0 - m));
		} else {

			// Body
			if ( fabs( x - xp ) >= (1.0 / m)) {
				return 0.0;
			} else {
				return (1 - m * fabs(x - xp));
			}
		}
	}
}


//! Polynomial basis function
double AdaptiveSparseGrid::LagrangePoly(double x, int i, int j) {
	if ( i == 1 ) {
		return 1.0;

	} else {
		int npts = pow(2.0, (i - 1)) + 1;
		Array<double> points; points.redim(npts);
		GetChebyshevPoints(i, points);
		double h = 1.0;
		double temp1, temp2;

		for ( int k = 1; k <= j - 1; k++) {
			h *= ( x - points(k)) / (points(j) - points(k));
		}

		for ( int k = j + 1; k <= npts; k++) {
			h *= ( x - points(k)) / (points(j) - points(k));
		}
		return h;

	}
}

//! Get the extreme of the Chebyshev points
void AdaptiveSparseGrid::GetChebyshevPoints(int i, Array<double>& points) {
	int npts = pow(2.0, (i - 1)) + 1;

	for ( int j = 1; j <= npts; j++) {
		points(j) = (-cos(M_PI * (j - 1) / (npts - 1)) + 1) / 2.0;
	}
}

//! calculate the value: supurlus X integral of the basis function
void  AdaptiveSparseGrid::SpIntegrate(double* value) {
	double* local_stat = new double[TotalDof];

	for ( int i = 0; i < TotalDof; i++)
	{ local_stat[i] = 0.0; }

	double mul;

	deque<AdaptiveData*>::reverse_iterator it1;
	set<AdaptiveNodeData*, AdaptiveNodeDataCompare>::iterator it2;

	//! Linear basis function
	if ( gridType == 1) {
		for (it1 = SparseGrid.rbegin(); it1 != SparseGrid.rend(); ++it1) {
			mul = 1.0;

			for ( int j = 1; j <= dim; j++) {
				mul *= LinearBasisVolumeIntegral((*(*it1)->index)(j));
			}

			for (it2 = (*it1)->NodeData.begin(); it2 != (*it1)->NodeData.end(); ++it2) {
				for ( int i = 0; i < TotalDof; i++) {
					local_stat[i] +=  mul * (((*it2)->surplus)[i]);
				}
			}
		}
	}

	//! Polynomial basis function
	else if ( gridType == 2) {
		//Calulate the Gauss-Legendre quadrature points
		calculate_wts_and_quad_pts(0.0, 1.0);

		for (it1 = SparseGrid.rbegin(); it1 != SparseGrid.rend(); ++it1) {
			for (it2 = (*it1)->NodeData.begin(); it2 != (*it1)->NodeData.end(); ++it2) {
				mul = 1.0;

				for ( int j = 1; j <= dim; j++) {
					mul *= PolynomialVolumeIntegral((*(*it1)->index)(j), (*(*it2)->index)(j));
				}

				for ( int i = 0; i < TotalDof; i++) {
					local_stat[i] +=  mul * (((*it2)->surplus)[i]);
				}
			}
		}

		// #AE :  ToDo Pfluger integral

	} else if ( gridType == 3) {
		for (it1 = SparseGrid.rbegin(); it1 != SparseGrid.rend(); ++it1) {
			mul = 1.0;

			for ( int j = 1; j <= dim; j++) {
				mul *= LinearBasisVolumeIntegral((*(*it1)->index)(j));
			}

			for (it2 = (*it1)->NodeData.begin(); it2 != (*it1)->NodeData.end(); ++it2) {
				for ( int i = 0; i < TotalDof; i++) {
					local_stat[i] +=  mul * (((*it2)->surplus)[i]);
				}
			}
		}
	} else if ( gridType == 4) { // QuadPOly
		for (it1 = SparseGrid.rbegin(); it1 != SparseGrid.rend(); ++it1) {
			mul = 1.0;

			for ( int j = 1; j <= dim; j++) {
				mul *= LinearBasisVolumeIntegral((*(*it1)->index)(j));
			}

			for (it2 = (*it1)->NodeData.begin(); it2 != (*it1)->NodeData.end(); ++it2) {
				for ( int i = 0; i < TotalDof; i++) {
					local_stat[i] +=  mul * (((*it2)->surplus)[i]);
				}
			}
		}
	}

	int error = MPI_Allreduce(local_stat, value, TotalDof, MPI_DOUBLE, MPI_SUM, mpiCOMM);
	//if (error) { cout << "ERORR 806 : " << error <<  "|" << rank << endl; }

	delete[] local_stat;
}

void  AdaptiveSparseGrid::ComputeMeanAndVar(char* filename) {
	double* stat = new double[TotalDof];
	//! Calculate the variance first. It depends on how you store the hierarchical surplus
	//! Here I store it as mean variance mean variance ......alternativly.
	SpIntegrate(stat);

	for ( int i = 1; i < TotalDof; i += 2) {
		stat[i]  -=  pow(stat[i - 1], 2);
	}

	StoreMeanAndVar(stat, filename);

}

//! Compute the integral of plynomial basis function through Gauss-Legendre quadrature
double AdaptiveSparseGrid::PolynomialVolumeIntegral(int i, int j) {
	double mul = 0.0;

	for ( int k = 0; k < no_quad_pts; k++) {
		mul += wts[k] * LagrangePoly(abs[k], i, j);

	}

	return mul;

}

//! integral of the linear basis function
//! if level = 1,  value is 1.0
//! if level = 2,  value is 0.25
//! if others,   value is 2^(1-level).
double AdaptiveSparseGrid::LinearBasisVolumeIntegral(int level) {
	double mul;

	if (  level == 1)
	{ mul = 1.0; }

	else if ( level == 2)
	{ mul = 0.25; }

	else
	{ mul = pow(2.0, (1 - level) / 1.0); }

	return mul;

}

//!Calulate the Gauss-Legendre quadrature points
void   AdaptiveSparseGrid::calculate_wts_and_quad_pts(const double low, const double up) {
	// Local variables
	int i;
	int avg_quad_pts;
	double mean_x, mean_x_diff;
	double guess_root;

	no_quad_pts = (int)((3 * Lmax + 1) / 2) + 1 ;
	abs = new double[no_quad_pts];
	wts = new double[no_quad_pts];

	// NR loop variables
	int j;
	double p1, p2, p3, pp;
	double z1;

	// Error tolerance for NR iterations
#define EPS 1.0e-14

	// self-explanatory
	avg_quad_pts = (no_quad_pts + 1) / 2;
	mean_x = 0.5 * (low + up);
	mean_x_diff = 0.5 * (up - low);

	for (i = 0; i < avg_quad_pts; i++) {

		// Guess roots for Legendre polynomials
		guess_root = cos(3.141592654 * (i + 0.75) / (no_quad_pts + 0.5));

		// Newton-Raphson loop to refine the roots of Legendre polynomials
		do {
			p1 = 1.0;
			p2 = 0.0;

			for (j = 0; j < no_quad_pts; j++) {
				p3 = p2;
				p2 = p1;
				p1 = ((2.0 * j + 1.0) * guess_root * p2 - j * p3) / (j + 1);
			}

			pp = no_quad_pts * (guess_root * p1 - p2) / (guess_root * guess_root - 1.0);
			z1 = guess_root;
			guess_root = z1 - p1 / pp;
		} while (fabs(guess_root - z1) > EPS);

		// Symmetric roots about the mid pt
		abs[i] = mean_x - mean_x_diff * guess_root;
		abs[no_quad_pts - i - 1] = mean_x + mean_x_diff * guess_root;

		// quad weights are symmetric
		wts[i] = 2.0 * mean_x_diff / ((1.0 - guess_root * guess_root) * pp * pp);
		wts[no_quad_pts - i - 1] = wts[i];
	}

#undef EPS

	//// output the abscissa and weights for the Legendre polynomials
	//if(1){
	//	printf("PRINTING THE WEIGHTS AND ABSCISSA FOR GAUSS-LEGENDRE QUADRATURE \n");
	//	printf("   weights    abscissa    \n");
	//	for(i=0;i<no_quad_pts;i++){
	//		printf("%8.3f    %8.3f  \n",wts[i],abs[i]);
	//	}
	//	printf("\n\n\n");
	//}
	//exit(1);
}

int AdaptiveSparseGrid::NumberOfPointsInThisProc() {
	int sum;
	int local_sum = 0.0;

	deque<AdaptiveData*>::reverse_iterator it1;

	for (it1 = SparseGrid.rbegin(); it1 != SparseGrid.rend(); ++it1) {
		local_sum +=  (*it1)->NodeData.size();
	}

	return local_sum;
}

//! Return the number of points in the current sparse grid
int AdaptiveSparseGrid::NumberOfPoints() {
	int local_sum, sum;

	local_sum = NumberOfPointsInThisProc();

	int error = MPI_Allreduce(&local_sum, &sum, 1, MPI_INT, MPI_SUM, mpiCOMM);
	//if (error) { cout << "ERORR 947 : " << error <<  "|" << rank << endl; }

	return sum;
}

//! Return the number of active points in the active index
int   AdaptiveSparseGrid::NumberOfActivePoints() {
	int sum = 0.0;

	set<AdaptiveData*, AdaptiveDataCompare>::iterator it1;

	for (it1 = ActiveIndex.begin(); it1 != ActiveIndex.end(); ++it1)
	{ sum +=  (*it1)->NodeData.size(); }

	return sum;
}

//! Plot the sparse grid in Tecplot format. You may overload this function
//! to write the grid in your own format
void AdaptiveSparseGrid::PlotSparseGrid(char* filename) {
	deque<AdaptiveData*>::reverse_iterator it1;
	set<AdaptiveNodeData*, AdaptiveNodeDataCompare>::iterator it2;

	int number = NumberOfPoints();

	MPI_Status status;

	if (rank == 0) {

		FILE* fp = fopen(filename, "w");

		if ( dim > 3) {
			cout << "ERROR: Dimension can not exceed 3!\n";
			//return;

		} else {
			fprintf(fp, "VARIABLES =");

			for ( int i = 1; i <= dim; i++) {
				fprintf(fp, "\"%d\"", i);
			}

			if ( dim == 1) {
				fprintf(fp, "\"L\"");
			}

			fprintf(fp, "\n");
			fprintf(fp, "ZONE I=%d, F=POINT\n", number);
		}

		int local_number = NumberOfPointsInThisProc();


		int totalnumber = local_number * dim;
		double* local_Coord = new double[totalnumber];

		int n = 0;

		for (it1 = SparseGrid.rbegin(); it1 != SparseGrid.rend(); ++it1) {

			for (it2 = (*it1)->NodeData.begin(); it2 != (*it1)->NodeData.end(); ++it2) {
				for ( int i = 0; i < dim ; i++) {
					local_Coord[n * dim + i] = IndextoCoordinate((*(*it1)->index)(i + 1), (*(*it2)->index)(i + 1));
				}

				n++;
			}
		}

		for ( int i = 0 ; i < local_number; i++) {

			double x;

			for ( int j = 0; j < dim ; j++) {
				x = local_Coord[i * dim + j];
				fprintf(fp, "%e\t ", x );
			}

			if ( dim == 1) {
				fprintf(fp, "%d\t ", CoordinateToIndex(x) );
			}

			fprintf(fp, "\n");
		}


		delete[] local_Coord;

		for ( int p = 1; p < size; p++) {
			MPI_Recv(&local_number, 1 , MPI_INT, p , 0 , mpiCOMM, &status);
			int totalnumber = local_number * dim;
			local_Coord = new double[totalnumber];
			MPI_Recv(local_Coord , totalnumber, MPI_DOUBLE, p , 0 , mpiCOMM, &status);

			for ( int i = 0 ; i < local_number; i++) {
				double x;

				for ( int j = 0; j < dim ; j++) {
					x = local_Coord[i * dim + j];
					fprintf(fp, "%e\t ", x );
				}

				if ( dim == 1)
				{ fprintf(fp, "%d\t ", CoordinateToIndex(x) ); }

				fprintf(fp, "\n");
			}

			delete[] local_Coord;
		}

		fclose(fp);

	} else {
		//! Get number of points from this processor
		int local_number = NumberOfPointsInThisProc();
		//! Send the number of points to the root
		MPI_Send(&local_number, 1, MPI_INT, 0, 0, mpiCOMM);

		int totalnumber = local_number * dim;
		double* local_Coord = new double[totalnumber];

		int n = 0;

		for (it1 = SparseGrid.rbegin(); it1 != SparseGrid.rend(); ++it1) {

			for (it2 = (*it1)->NodeData.begin(); it2 != (*it1)->NodeData.end(); ++it2) {
				for ( int i = 0; i < dim ; i++) {
					local_Coord[n * dim + i] = IndextoCoordinate((*(*it1)->index)(i + 1), (*(*it2)->index)(i + 1));
				}

				n++;
			}
		}

		MPI_Send(local_Coord, totalnumber, MPI_DOUBLE, 0, 0, mpiCOMM);
		delete[] local_Coord;
	}
}


//! Store the hierarchical surplus
//! The basic format of the file is
//! dimension,  number of point,  total degree of freedom, interpolation depth
//! i_1..... i_dim  j_1 .... j_dim surplus_1 .... surplus_TotalDof
//! You need to give the index of the surplus which is going to store.
void AdaptiveSparseGrid::StoreSurplus(AdaptiveARRAY<int>& index, char* filename) {

	deque<AdaptiveData*>::reverse_iterator it1;
	set<AdaptiveNodeData*, AdaptiveNodeDataCompare>::iterator it2;

	int number = NumberOfPoints();
	MPI_Status status;

	int Dof = index.n;

	if (rank == 0) {

		FILE* fp = fopen(filename, "w");
		fprintf(fp, "%d\t %d\t %d\t %d\t %d\n", dim, number, Dof, InterpolationLevel(), gridType);

		//! Get number of points from this processor
		int local_number = NumberOfPointsInThisProc();
		//! Get all of the information from this processors
		int totalnumber = local_number * (2 * dim + Dof);
		double* local_Coord = new double[totalnumber];

		int n = 0;

		for (it1 = SparseGrid.rbegin(); it1 != SparseGrid.rend(); ++it1) {

			for (it2 = (*it1)->NodeData.begin(); it2 != (*it1)->NodeData.end(); ++it2) {
				for ( int i = 0; i < dim ; i++) {
					local_Coord[n * (2 * dim + Dof) + i] = (*(*it1)->index)(i + 1);
				}

				for ( int i = 0; i < dim ; i++) {
					local_Coord[n * (2 * dim + Dof) + dim + i] = (*(*it2)->index)(i + 1);
				}

				for ( int i = 0; i < Dof; i++) {
					local_Coord[ n * (2 * dim + Dof) + 2 * dim + i] = (*it2)->surplus[index(i + 1)];
				}

				n++;
			}

		}

		for ( int i = 0 ; i < local_number; i++) {
			for ( int j = 0; j < dim ; j++) {
				int  index = local_Coord[i * (2 * dim + Dof) + j];
				fprintf(fp, "%d\t ", index );
			}

			for ( int j = 0; j < dim ; j++) {
				int  index = local_Coord[i * (2 * dim + Dof) + dim + j];
				fprintf(fp, "%d\t ", index );
			}

			for ( int j = 0; j < Dof; j++)
				// #AE increase percision
				//{ fprintf(fp, "%e\t ", local_Coord[i * (2 * dim + Dof) + 2 * dim + j]); }
			{ fprintf(fp, "%.15f\t ", local_Coord[i * (2 * dim + Dof) + 2 * dim + j]); }

			{ fprintf(fp, "\n"); }
		}

		delete[] local_Coord;


		//From other processors
		for ( int p = 1; p < size; p++) {

			MPI_Recv(&local_number, 1 , MPI_INT, p , 0 , mpiCOMM, &status);
			totalnumber = local_number * (2 * dim + Dof);
			local_Coord = new double[totalnumber];
			MPI_Recv(local_Coord , totalnumber, MPI_DOUBLE, p , 0 , mpiCOMM, &status);

			for ( int i = 0 ; i < local_number; i++) {
				for ( int j = 0; j < dim ; j++) {
					int  index = local_Coord[i * (2 * dim + Dof) + j];
					fprintf(fp, "%d\t ", index );
				}

				for ( int j = 0; j < dim ; j++) {
					int  index = local_Coord[i * (2 * dim + Dof) + dim + j];
					fprintf(fp, "%d\t ", index );
				}

				for ( int j = 0; j < Dof; j++)
					// #AE increase percision
					//{ fprintf(fp, "%e\t ", local_Coord[i * (2 * dim + Dof) + 2 * dim + j]); }
				{ fprintf(fp, "%.15f\t", local_Coord[i * (2 * dim + Dof) + 2 * dim + j]); }

				{ fprintf(fp, "\n"); }
			}

			delete[] local_Coord;
		}

		fclose(fp);

	} else {

		//! Get number of points from this processor
		int local_number = NumberOfPointsInThisProc();
		//! Send the number of points to the root
		MPI_Send(&local_number, 1, MPI_INT, 0, 0, mpiCOMM);

		//! Get all of the information from this processors
		int totalnumber = local_number * (2 * dim + Dof);
		double* local_Coord = new double[totalnumber];

		int n = 0;

		for (it1 = SparseGrid.rbegin(); it1 != SparseGrid.rend(); ++it1) {

			for (it2 = (*it1)->NodeData.begin(); it2 != (*it1)->NodeData.end(); ++it2) {
				for ( int i = 0; i < dim ; i++) {
					local_Coord[n * (2 * dim + Dof) + i] = (*(*it1)->index)(i + 1);
				}

				for ( int i = 0; i < dim ; i++) {
					local_Coord[n * (2 * dim + Dof) + dim + i] = (*(*it2)->index)(i + 1);
				}

				for ( int i = 0; i < Dof; i++) {
					local_Coord[ n * (2 * dim + Dof) + 2 * dim + i] = (*it2)->surplus[index(i + 1)];
				}

				n++;
			}

		}

		MPI_Send(local_Coord, totalnumber, MPI_DOUBLE, 0, 0, mpiCOMM);
		delete[] local_Coord;

	}

	MPI_Barrier(mpiCOMM);
}


void AdaptiveSparseGrid::Restart(char* filename) {
	Cleanup();

	int i, j;
	ifstream infile;

	if (rank == 0)
	{ infile.open(filename); }

	if (rank == 0) {
		if (!infile.is_open()) {
			printf("Error opening file: %s\n", filename);
			return;
		}
	}



	if (rank == 0) { printf("\nNow Loading the %s \n\n", filename); }

	int number;
	AdaptiveARRAY<int>* px1, *px2;
	AdaptiveData* pData;
	AdaptiveNodeData* pNodeData;

	AdaptiveGrid IndexSet;

	pair<set<AdaptiveData*, AdaptiveDataCompare>::iterator, bool> plt;


	if ( type == 1) {

		if (rank == 0) {
			infile >> dim;
			infile >> number;
			infile >> TotalDof;
			infile >> L;
		}

		MPI_Bcast(&dim, 1, MPI_INT, 0 , mpiCOMM);
		MPI_Bcast(&number, 1, MPI_INT, 0 , mpiCOMM);
		MPI_Bcast(&TotalDof, 1, MPI_INT, 0 , mpiCOMM);
		MPI_Bcast(&L, 1, MPI_INT, 0 , mpiCOMM);

		for (int i = 1; i <= number; i++) {

			//!First load the i index;
			px1 = AdaptiveCoordAllocator.NewItem(dim);
			int sum = 0;

			if ( rank == 0) {
				for ( j = 1 ; j <= dim; j++) {
					infile >> (*px1)(j);
					sum += (*px1)(j);
				}
			}

			MPI_Bcast(px1->pData, dim, MPI_INT, 0 , mpiCOMM);
			MPI_Bcast(&sum, 1, MPI_INT, 0 , mpiCOMM);

			//!Then load the j index;
			px2 = AdaptiveCoordAllocator.NewItem(dim);

			if ( rank == 0) {
				for ( int j = 1 ; j <= dim; j++)
				{ infile >> (*px2)(j); }
			}

			MPI_Bcast(px2->pData, dim, MPI_INT, 0 , mpiCOMM);

			//!Load the surplus
			double* surplus = new double[TotalDof];

			if ( rank == 0) {
				for ( int j = 0; j < TotalDof; j++)
				{ infile >> surplus[j]; }
			}

			MPI_Bcast(surplus, TotalDof, MPI_DOUBLE, 0 , mpiCOMM);


			if ( sum != (L + dim) ) {
				if ( (i % size) == rank) {

					pData = AdaptiveDataAllocator.NewItem();
					pData-> index = px1;

					pNodeData = AdaptiveNodeDataAllocator.NewItem();
					pNodeData->surplus = new double[TotalDof];
					pNodeData->index = px2;

					for ( int k = 0; k < TotalDof; k++)
					{ pNodeData->surplus[k] = surplus[k]; }


					//! If the multi-index i has already existed, plt.second = false
					plt = IndexSet.insert(pData);

					if ( plt.second == false ) {

						AdaptiveDataAllocator.DeleteItem(pData);
						AdaptiveCoordAllocator.DeleteItem(px1);

						(*(plt.first))->NodeData.insert(pNodeData);
					}

					//! If the multi-index i has not existed, insert the multi-index j directly
					else if (plt.second == true)
					{ pData->NodeData.insert(pNodeData); }

				} else {
					AdaptiveCoordAllocator.DeleteItem(px1);
					AdaptiveCoordAllocator.DeleteItem(px2);
				}


				delete[] surplus;
			}
		}

		if (rank == 0)
		{ infile.close(); }

		//! Insert all of the points to the sparse grid.
		SparseGrid.insert(SparseGrid.begin(), IndexSet.begin(), IndexSet.end());
		IndexSet.clear();

		UpdateError();

		ifstream infile1;

		if (rank == 0)
		{ infile1.open(filename); }

		if (rank == 0) {
			infile1 >> dim;
			infile1 >> number;
			infile1 >> TotalDof;
			infile1 >> L;
		}


		for (int i = 1; i <= number; i++) {

			px1 = AdaptiveCoordAllocator.NewItem(dim);
			int sum = 0;

			if ( rank == 0) {
				for ( j = 1 ; j <= dim; j++) {
					infile1 >> (*px1)(j);
					sum += (*px1)(j);
				}
			}

			MPI_Bcast(px1->pData, dim, MPI_INT, 0 , mpiCOMM);
			MPI_Bcast(&sum, 1, MPI_INT, 0 , mpiCOMM);

			//!Then load the j index;
			px2 = AdaptiveCoordAllocator.NewItem(dim);

			if ( rank == 0) {
				for ( int j = 1 ; j <= dim; j++)
				{ infile1 >> (*px2)(j); }
			}

			MPI_Bcast(px2->pData, dim, MPI_INT, 0 , mpiCOMM);

			//!Load the surplus
			double* surplus = new double[TotalDof];

			if ( rank == 0) {
				for ( int j = 0; j < TotalDof; j++)
				{ infile1 >> surplus[j]; }
			}

			MPI_Bcast(surplus, TotalDof, MPI_DOUBLE, 0 , mpiCOMM);


			//! If it is the last level, generate the active index
			if ( sum == (L + dim) ) {

				int refine = IsRefine(surplus, px1, px2);

				if ( refine)
				{ Refine(px1, px2); }

				if ( (i % size) == rank) {

					pData = AdaptiveDataAllocator.NewItem();
					pData-> index = px1;

					pNodeData = AdaptiveNodeDataAllocator.NewItem();
					pNodeData->surplus = new double[TotalDof];
					pNodeData->index = px2;

					for ( int k = 0; k < TotalDof; k++)
					{ pNodeData->surplus[k] = surplus[k]; }


					//! If the multi-index i has already existed, plt.second = false
					plt = IndexSet.insert(pData);

					if ( plt.second == false ) {

						AdaptiveDataAllocator.DeleteItem(pData);
						AdaptiveCoordAllocator.DeleteItem(px1);

						(*(plt.first))->NodeData.insert(pNodeData);
					}

					//! If the multi-index i has not existed, insert the multi-index j directly
					else if (plt.second == true)
					{ pData->NodeData.insert(pNodeData); }

				} else {
					AdaptiveCoordAllocator.DeleteItem(px1);
					AdaptiveCoordAllocator.DeleteItem(px2);
				}


				delete[] surplus;
			}
		}

		if (rank == 0)
		{ infile1.close(); }

		//! Insert all of the points to the sparse grid.
		SparseGrid.insert(SparseGrid.begin(), IndexSet.begin(), IndexSet.end());

		IndexSet.clear();


		L += 1;

	}

	else {

		if (rank == 0) {
			infile >> dim;
			infile >> number;
			infile >> TotalDof;
			infile >> L;
		}

		MPI_Bcast(&dim, 1, MPI_INT, 0 , mpiCOMM);
		MPI_Bcast(&number, 1, MPI_INT, 0 , mpiCOMM);
		MPI_Bcast(&TotalDof, 1, MPI_INT, 0 , mpiCOMM);
		MPI_Bcast(&L, 1, MPI_INT, 0 , mpiCOMM);

		for (int i = 1; i <= number; i++) {

			px1 = AdaptiveCoordAllocator.NewItem(dim);
			int sum = 0;

			if ( rank == 0) {
				for ( j = 1 ; j <= dim; j++) {
					infile >> (*px1)(j);
					sum += (*px1)(j);
				}
			}

			MPI_Bcast(px1->pData, dim, MPI_INT, 0 , mpiCOMM);
			MPI_Bcast(&sum, 1, MPI_INT, 0 , mpiCOMM);

			//!Then load the j index;
			px2 = AdaptiveCoordAllocator.NewItem(dim);

			if ( rank == 0) {
				for ( int j = 1 ; j <= dim; j++)
				{ infile >> (*px2)(j); }
			}

			MPI_Bcast(px2->pData, dim, MPI_INT, 0 , mpiCOMM);

			//!Load the surplus
			double* surplus = new double[TotalDof];

			if ( rank == 0) {
				for ( int j = 0; j < TotalDof; j++)
				{ infile >> surplus[j]; }
			}

			MPI_Bcast(surplus, TotalDof, MPI_DOUBLE, 0 , mpiCOMM);


			//! If it is the last level, generate the active index
			if ( sum == (L + dim) ) {
				int refine = IsRefine(surplus, px1, px2);

				if ( refine)
				{ Refine(px1, px2); }
			}

			if ( (i % size) == rank) {

				pData = AdaptiveDataAllocator.NewItem();
				pData-> index = px1;

				pNodeData = AdaptiveNodeDataAllocator.NewItem();
				pNodeData->surplus = new double[TotalDof];
				pNodeData->index = px2;

				for ( int k = 0; k < TotalDof; k++)
				{ pNodeData->surplus[k] = surplus[k]; }


				//! If the multi-index i has already existed, plt.second = false
				plt = IndexSet.insert(pData);

				if ( plt.second == false ) {

					AdaptiveDataAllocator.DeleteItem(pData);
					AdaptiveCoordAllocator.DeleteItem(px1);

					(*(plt.first))->NodeData.insert(pNodeData);
				}

				//! If the multi-index i has not existed, insert the multi-index j directly
				else if (plt.second == true)
				{ pData->NodeData.insert(pNodeData); }

			} else {
				AdaptiveCoordAllocator.DeleteItem(px1);
				AdaptiveCoordAllocator.DeleteItem(px2);
			}

			delete[] surplus;

		}

		if (rank == 0)
		{ infile.close(); }

		//! Insert all of the points to the sparse grid.
		SparseGrid.insert(SparseGrid.begin(), IndexSet.begin(), IndexSet.end());

		IndexSet.clear();

		L += 1;
	}
}
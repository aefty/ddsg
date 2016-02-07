/**
 * File has been edited by Aryan E. in sections noted by the tag  #AE:
 */


#include "Post.h"

Post::Post() {

	MPI_Comm_rank(mpiCOMM, &rank);
	MPI_Comm_size(mpiCOMM, &size);
}


void Post::LoadData(char* filename) {
	int i, j;

	ifstream infile(filename);

	if (infile.is_open()) {
		infile >> dim;
		infile >> nno;
		infile >> TotalDof;
		infile >> Level;
		infile >> gridType;

		int nsd = 2 * dim;
		index.redim(nno, nsd); index.fill(0);
		surplus.redim(nno, TotalDof); surplus.fill(0.0);
		i = 1;
		while (infile) {
			if (i == nno + 1) { break; }
			for ( j = 1; j <= nsd; j++)
			{    infile >> index(i, j); }
			for ( j = 1; j <= TotalDof; j++)
			{    infile >> surplus(i, j); }

			i++;
		}
		infile.close();
	} else {
		if (rank == 0)
		{ printf("Error opening file: %s\n", filename); }
		exit(1);
	}

}


double Post::IndextoCoordinate(int i, int j ) {
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
			//#AE: Change in the slope of the line
			m = (1 << (i)) + 1;
			return (j - 1.0) / (m - 1) ;
			break;
		//#AE: Quadpoly
		case 4:
			m = (1 << (i)) + 1;
			return (j - 1.0) / (m - 1) ;
			break;
		//#AE: throw error when grid type not defined
		default:
			cout << "!! gridType not valid ( " << gridType << " ) !!" << endl;
			exit(0);
	}
}


//#AE: Select basis function based gridType
double Post::BasisFunction(double x, int i, int j) {

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
			cout << "!! gridType not valid ( " << gridType << " ) !!" << endl;
			exit(0);
	}
};

double Post::LinearBasis(double x, int i, int j) {

	//First Level
	if (  i == 1 ) {
		return 1.0;
	} else {
		double m = pow(2.0, (i - 1));
		double xp = IndextoCoordinate(i, j );
		if ( fabs( x - xp ) >= (1.0 / m))
		{ return 0.0;}
		else {
			return (1 - m * fabs(x - xp));
		}
	}
}


double Post::PolyBasis(double x, int i, int j) {

	if (i < 3) {
		return FlipUpBasis(x, i, j);
	} else {
		double xp, x1, x2, temp;
		double m = pow(2.0, i - 1);

		// Start or end
		if (j == 2 && x <= 0.5 ) {
			x1 = 1. / m;
			x2 = 1.0 - x1;
			xp = 0.5;
			//temp = -i * (x - x1) * (x - x2) / ((xp - x1) * (xp - x2));
			temp = (2.0) / (x1 - x1 * x1) * (x - x1) * (x - x2);
		} else if (j == (1 << i) && x > 0.5) {
			x1 = 1. / m;
			x2 = 1.0 - x1;
			xp = 0.5;
			//temp = -i * (x - x1) * (x - x2) / ((xp - x1) * (xp - x2));
			temp = (2.0) / (x1 - x1 * x1) * (x - x1) * (x - x2);
		} else {
			xp = IndextoCoordinate(i, j); // ith index / m + 0.5*m  (ith index is not 2^i)
			x1 = xp - 0.5 / m;
			x2 = xp + 0.5 / m;
			temp = (x - x1) * (x - x2) / ((xp - x1) * (xp - x2));
		}

		//cout << "(x-" << x1 << ")(x-" << x2 << ")/(" << xp << "-" << x1 << ")(" << xp << "-" << x2 << ")" << " ... " << IndextoCoordinate(i, j) << endl;

		if (temp > 0) {
			return temp;
		} else {
			return 0.0;
		}
	}
}

//! Polynomial basis function
double Post::LagrangePoly(double x, int i, int j) {
	if ( i == 1 )
	{return 1.0;}
	else {
		int npts = pow(2.0, (i - 1)) + 1;
		Array<double> points; points.redim(npts);
		GetChebyshevPoints(i, points);
		double h = 1.0;
		for ( int k = 1; k <= j - 1; k++)
		{ h *= ( x - points(k)) / (points(j) - points(k)); }
		for ( int k = j + 1; k <= npts; k++)
		{ h *= ( x - points(k)) / (points(j) - points(k)); }
		return h;
	}
}


//! Get the extreme of the Chebyshev points
void Post::GetChebyshevPoints(int i, Array<double>& points) {
	int npts = pow(2.0, (i - 1)) + 1;
	for ( int j = 1; j <= npts; j++)
	{ points(j) = (-cos(M_PI * (j - 1) / (npts - 1)) + 1) / 2.0; }
}


/**
 * #AE
 * FlipUpBasis Function
 * @param  x [X val]
 * @param  i [Level Depth]
 * @param  j [Basis Function Index]
 * @return   [Value]
 */
double Post::FlipUpBasis(double x, int i, int j) {

	if (  i == 1 ) {
		return 1.0;
	} else {

		double m = pow(2.0, i);

		//Position marker - 0:middle, 1:start, 2:end;
		int pos = 1 * ( (j == 2) || (i == 2 && j == 1) ) +
		          2 * ( (j == 1 << (i)) || (i == 2 && j == 3) );

		double xp = IndextoCoordinate(i, j);

		if ( fabs( x - xp ) > (1.0 / m)) {
			return 0.0;
		} else {

			if (pos == 1) {               // Start
				return (1 - m * (x - xp) );
			} else if (pos == 2) {        // End
				return (m * (x - xp) + 1 );
			} else {                      // Middle
				return (1 - m * fabs(x - xp));
			}
		}
	}
}


void Post::Interpolate(Array<double>& x, double* value) {
	double  temp;

	for ( int i = 0 ; i < TotalDof; i++) {
		value[i] = 0.0;
	}

	for ( int i = 1; i <= nno; i++) {

		temp = 1.0;
		for ( int j = 1; j <= dim; j++) {
			temp *= BasisFunction(x(j), index(i, j), index(i, j + dim));
			if ( temp == 0.0) { break; }
		}
		for ( int j = 1; j <= TotalDof; j++)
		{ value[j - 1] += temp * surplus(i, j); }
	}
}

void Post::Interpolate(double* x, double* value) {
	Array<double> px;
	px.redim(dim);
	for ( int i = 1; i <= dim; i++) {
		px(i) = x[i - 1];
	}

	Interpolate(px, value);

}

void Post::Integrate(double* value) {
	double temp;

	for ( int i = 0 ; i < TotalDof; i++)
	{ value[i] = 0.0; }

	for ( int i = 1; i <= nno; i++) {

		temp = 1.0;
		for ( int j = 1; j <= dim; j++) {

			temp *= LinearBasisVolumeIntegral(index(i, j));

		}
		for ( int j = 1; j <= TotalDof; j++)
		{ value[j - 1] += temp * surplus(i, j); }
	}

}

//! integral of the linear basis function
//! if level = 1,  value is 1.0
//! if level = 2,  value is 0.25
//! if others,   value is 2^(1-level).
double Post::LinearBasisVolumeIntegral(int level) {
	double mul;
	if (  level == 1)
	{ mul = 1.0; }
	else if ( level == 2)
	{ mul = 0.25; }
	else
	{ mul = pow(2.0, (1 - level) / 1.0); }

	return mul;

}

//! Plot the sparse grid in Tecplot format. You may overload this function
//! to write the grid in your own format
void Post::PlotSparseGrid(char* filename) {
	FILE* fp = fopen(filename, "w");
	if ( dim > 3) {
		cout << "ERROR: Dimension can not exceed 3!" << endl;
		exit(1);
	} else {

		fprintf(fp, "VARIABLES =");
		for ( int i = 1; i <= dim; i++) {
			fprintf(fp, "\"%d\"", i);
		}

		if ( dim == 1) {
			fprintf(fp, "\"L\"");
		}

		fprintf(fp, "\n");
		fprintf(fp, "ZONE I=%d, F=POINT\n", nno);

	}

	for ( int i = 1 ; i <= nno; i++) {
		double x;
		for ( int j = 1; j <= dim ; j++) {
			x = IndextoCoordinate(index(i, j), index(i, j + dim) );
			fprintf(fp, "%e\t ", x );
		}
		if ( dim == 1)
		{ fprintf(fp, "%d\t ", index(i, 1) ); }
		fprintf(fp, "\n");
	}

	fclose(fp);

}
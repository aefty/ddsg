#include <math.h>
#include <time.h>
#include <mpi.h>
#include <iomanip>
#include <stdio.h>
#include "lib/SG.h"


void testFunc(double* x , int dim, double* rtrn) {
    double rtn = 0.0;
    double pi = acos(-1.0L);
    rtrn[0] = fmax(((x[0] * (x[1] * x[1]) - 1.0 / pi) * pi / (pi - 1.0) + 0.4), 0);
    rtrn[0] = sqrt(rtrn[0]);
};


int unitTest();


int main(int argc, char **args) {

    MPI_Init (&argc, &args);
    {
        unitTest();
    }

    MPI_Finalize ();

    return 0;
}


int unitTest() {

    //int fnum = 0; //  function 0 = sqrt(max(..)...) ,  function 99 = y^2 + x^2
    int dim = 2;
    int Lmax = 6; // Refinement levels after level 1 -> total level = Lmax+1;
    double epsilon = 0; //1e-4;//0;//1e-6; //if epsilon = 0, sparse grid with level Lmax
    int type = 3;
    int h = 7; // Discretization ( 2^h)

    SG sp(testFunc, dim, Lmax, epsilon, type);
    sp.solve("results/surplus.data");
    sp.interpolate(h, "results/interpolate.data"); // Create Interpolate from Pre-Proc data. This also create exact values in "orginal.plt"
    sp.exact(h, "results/exact.data");

    /**
     * Post Processing of Surplus File.
     */
    {
        char* surplusFile = "results/surplus.data";
        char* interplateFile = "results/interpolate_post.data";
        char* gridFile = "results/grid.data";

        Post post_interpolate;
        post_interpolate.LoadData(surplusFile);

        {
            ofstream output;
            output.open(interplateFile);

            double *fvalue = new double[post_interpolate.TotalDof];
            AdaptiveARRAY<double> x;
            x.redim(post_interpolate.dim);

            for ( x(1) = 0.0; x(1) <= 1.0; x(1) += 1.0 / (1 << h)) {
                for ( x(2) = 0.0; x(2) <= 1.0; x(2) += 1.0 / (1 << h)) {
                    post_interpolate.Interpolate(x, fvalue);
                    for ( int i = 0; i < post_interpolate.TotalDof; i++) {
                        output << "\t" << fvalue[i];
                    }
                }
                output << "\n";
            }

            output.close();
        }

        // Post Processing Plot of Grid
        post_interpolate.PlotSparseGrid(gridFile);
    }

    /**
     * General Unit Tests
     */
    {
        // Basis Function Test
        if (1) {

            int set [4][8]  = {
                {0, -1, -1, -1, -1,  -1, -1, -1},
                {1,  3, -1, -1, -1,  -1, -1, -1},
                {2,  4,  6,  8, -1,  -1, -1, -1},
                {2,  4,  6,  8,  10, 12, 14, 16}
            };

            // Correct index if basisType 3 is used!
            if (type == 3) {
                set[1][0] = 2;
                set[1][1] = 4;
            }

            double z;

            double step = (1.0 / (1 << h));

            FILE *fp1;

            fp1 = fopen("results/baisFunction.data", "w");

            for (int i = 0; i < 4; ++i) {

                for (int j = 0; j < 8; ++j) {

                    if (set[i][j] < 0) {
                        break;
                    }

                    for (double x = 0; x <= 1.0; x += step ) {
                        z = sp.BasisFunction(x, i + 1, set[i][j]) ;

                        if (z > 0.0) {
                            fprintf(fp1, " %i\t %e\t %e\t \n", i + 1, x, z );
                        }
                    }
                }
            }

            fclose(fp1);
        }

        // Quick Terminal Test
        if (0) {
            int l = 3;
            int j = 4; // Key point
            double val = 0.0;

            for (double x = 0; x < 1.0; x += (1.0 / (1 << h))) {
                val = sp.BasisFunction(x, l, j);
                cout << std::setw(4)  << x << std::setprecision(3) << std::setw(val * 20 + 5) << val << endl;
            }
        }
    }

    return 1;
}